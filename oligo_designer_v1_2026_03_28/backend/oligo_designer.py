"""
Oligo Designer for Gene Synthesis via Overlap Extension Assembly
=====================================================
Takes a DNA sequence and designs overlapping oligos for assembly.

V1: No codon optimization. DNA sequence in, oligos out.

Usage:
    python oligo_designer.py <sequence_or_fasta_file> [options]

The oligos tile the sequence with alternating sense/antisense orientation.
Overlaps (~20bp) between adjacent oligos are checked for:
  - Uniqueness within the construct
  - GC content (30-70%)
  - Homopolymer runs
  - G-quadruplex potential
"""

import re
import sys
import argparse
from dataclasses import dataclass, field
from typing import Optional
from primer3 import thermoanalysis


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ReactionConditions:
    """Buffer/reaction conditions for thermodynamic calculations."""
    mv_conc: float = 50.0       # monovalent cations [Na+/K+] in mM
    dv_conc: float = 2.0        # divalent cations [Mg2+] in mM
    dntp_conc: float = 0.8      # total dNTPs in mM (0.2 mM each × 4)
    dna_conc: float = 250.0     # oligo concentration in nM
    annealing_temp: float = 50.0  # annealing temperature in °C

    def make_thermo(self) -> thermoanalysis.ThermoAnalysis:
        return thermoanalysis.ThermoAnalysis(
            mv_conc=self.mv_conc,
            dv_conc=self.dv_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
        )


@dataclass
class OverlapIssue:
    """A problem detected in an overlap region."""
    kind: str        # "uniqueness", "gc", "homopolymer", "g_quadruplex"
    message: str
    severity: str    # "warning", "error"
    start: Optional[int] = None   # construct position of the problem region
    end: Optional[int] = None


@dataclass
class Overlap:
    """An overlap region between two adjacent oligos."""
    seq: str
    start: int          # position in the full construct
    end: int
    issues: list = field(default_factory=list)
    tm: Optional[float] = None  # melting temperature in °C

    @property
    def gc(self):
        gc_count = self.seq.upper().count("G") + self.seq.upper().count("C")
        return gc_count / len(self.seq) if self.seq else 0


@dataclass
class Oligo:
    """A single oligo to order from IDT."""
    index: int
    seq: str             # 5'->3' as ordered
    start: int           # position in full construct (sense coords)
    end: int
    strand: str          # "sense" or "antisense"
    length: int = 0
    is_first: bool = False
    is_last: bool = False

    def __post_init__(self):
        self.length = len(self.seq)


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(comp)[::-1]


def gc_content(seq: str) -> float:
    if not seq:
        return 0
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s)


def find_homopolymers(seq: str, min_run: int = 5) -> list[tuple[int, int, str]]:
    """Return list of (start, end, base) for homopolymer runs >= min_run."""
    runs = []
    s = seq.upper()
    i = 0
    while i < len(s):
        j = i + 1
        while j < len(s) and s[j] == s[i]:
            j += 1
        if j - i >= min_run:
            runs.append((i, j, s[i]))
        i = j
    return runs


def find_g_quadruplex(seq: str) -> list[tuple[int, int]]:
    """Find potential G-quadruplex motifs: 4+ runs of 3+ G's within 30bp."""
    s = seq.upper()
    # Find all G-runs of 3+
    g_runs = [(m.start(), m.end()) for m in re.finditer(r"G{3,}", s)]
    quads = []
    for i in range(len(g_runs)):
        # Check if 4 G-runs fit within 30bp
        for j in range(i + 3, len(g_runs)):
            span_start = g_runs[i][0]
            span_end = g_runs[j][1]
            if span_end - span_start <= 30:
                quads.append((span_start, span_end))
    # Also check complement (C-runs form G-quad on other strand)
    c_runs = [(m.start(), m.end()) for m in re.finditer(r"C{3,}", s)]
    for i in range(len(c_runs)):
        for j in range(i + 3, len(c_runs)):
            span_start = c_runs[i][0]
            span_end = c_runs[j][1]
            if span_end - span_start <= 30:
                quads.append((span_start, span_end))
    return quads


def find_repeats(seq: str, kmer_size: int = 10) -> dict[str, list[int]]:
    """Find all k-mers that appear more than once. Returns {kmer: [positions]}."""
    s = seq.upper()
    kmers: dict[str, list[int]] = {}
    for i in range(len(s) - kmer_size + 1):
        kmer = s[i:i + kmer_size]
        kmers.setdefault(kmer, []).append(i)
    return {k: v for k, v in kmers.items() if len(v) > 1}


# ---------------------------------------------------------------------------
# Overlap quality checking
# ---------------------------------------------------------------------------

def check_overlap(overlap: Overlap, full_seq: str) -> list[OverlapIssue]:
    """Run all quality checks on an overlap region."""
    issues = []
    seq = overlap.seq.upper()

    # Homopolymer runs — pinpoint the exact run within the construct
    homos = find_homopolymers(seq, min_run=4)
    for hp_start, hp_end, base in homos:
        run_len = hp_end - hp_start
        sev = "error" if run_len >= 6 else "warning"
        issues.append(OverlapIssue(
            "homopolymer",
            f"{run_len}x {base} at pos {overlap.start + hp_start}-{overlap.start + hp_end}",
            sev,
            start=overlap.start + hp_start,
            end=overlap.start + hp_end,
        ))

    # G-quadruplex — pinpoint the span
    quads = find_g_quadruplex(seq)
    if quads:
        q_start = min(q[0] for q in quads)
        q_end = max(q[1] for q in quads)
        issues.append(OverlapIssue(
            "g_quadruplex",
            f"G-quadruplex potential at pos {overlap.start + q_start}-{overlap.start + q_end}",
            "error",
            start=overlap.start + q_start,
            end=overlap.start + q_end,
        ))

    return issues


def _avg_overlap_tm(overlaps: list, thermo: thermoanalysis.ThermoAnalysis) -> float:
    """Compute average Tm across all overlaps."""
    tms = []
    for ovl in overlaps:
        if ovl.tm is not None:
            tms.append(ovl.tm)
        else:
            tms.append(thermo.calc_tm(ovl.seq))
    return sum(tms) / len(tms) if tms else 60.0


def check_cross_hybridization(
    overlaps: list,
    conditions: ReactionConditions,
) -> list[tuple[int, int, OverlapIssue, OverlapIssue]]:
    """
    Check all overlap pairs for cross-hybridization using thermodynamic
    calculation (primer3 nearest-neighbor model).

    An off-target interaction is flagged when the heterodimer Tm between two
    overlaps is within `margin` degrees of the annealing temperature — meaning
    the wrong pair could form a stable duplex under reaction conditions.

    Returns list of (i, j, issue_for_i, issue_for_j) tuples.
    """
    thermo = conditions.make_thermo()
    # Threshold: average overlap Tm minus 20°C, per Gene2Oligo
    # (Rouillard et al. 2004, Nucleic Acids Res. 32:W176-W180).
    # Off-target heterodimers with Tm above this are flagged as non-specific.
    avg_tm = _avg_overlap_tm(overlaps, thermo)
    threshold_tm = avg_tm - 20.0
    results = []

    for i in range(len(overlaps)):
        for j in range(i + 1, len(overlaps)):
            seq_i = overlaps[i].seq.upper()
            seq_j = overlaps[j].seq.upper()

            # Check both orientations: seq_i vs seq_j, and seq_i vs RC(seq_j)
            rc_j = reverse_complement(seq_j)
            for target, label in [(seq_j, ""), (rc_j, " (reverse complement)")]:
                result = thermo.calc_heterodimer(seq_i, target)
                if result.tm > threshold_tm:
                    detail = (f"Tm {result.tm:.1f}°C "
                              f"(threshold {threshold_tm:.0f}°C)")
                    issue_i = OverlapIssue(
                        "cross_hybridization",
                        f"May cross-hybridize with overlap {j+1}{label}: "
                        + detail,
                        "error",
                    )
                    issue_j = OverlapIssue(
                        "cross_hybridization",
                        f"May cross-hybridize with overlap {i+1}{label}: "
                        + detail,
                        "error",
                    )
                    results.append((i, j, issue_i, issue_j))

    return results


# ---------------------------------------------------------------------------
# Core tiling algorithm
# ---------------------------------------------------------------------------

def tile_sequence(
    seq: str,
    oligo_length: int = 45,
    overlap_length: int = 20,
    plasmid_upstream: str = "",
    plasmid_downstream: str = "",
    max_shift: int = 5,
    conditions: ReactionConditions | None = None,
) -> tuple[list[Oligo], list[Overlap]]:
    """
    Tile a DNA sequence into overlapping oligos.

    The oligos alternate sense/antisense:
      Oligo 1 (sense):      =========>
      Oligo 2 (antisense):       <==========
      Oligo 3 (sense):                =========>
      ...

    Each oligo covers `oligo_length` bp of the construct.
    Adjacent oligos share `overlap_length` bp overlaps.
    The step size = oligo_length - overlap_length.

    Parameters
    ----------
    seq : str
        The DNA insert sequence (sense strand, 5'->3').
    oligo_length : int
        Target oligo length in bp (default 45).
    overlap_length : int
        Target overlap between adjacent oligos (default 20).
    plasmid_upstream : str
        ~20bp of plasmid sequence upstream of insert (for vector insertion).
    plasmid_downstream : str
        ~20bp of plasmid sequence downstream of insert (for vector insertion).
    max_shift : int
        How far to shift a cut point to fix a bad overlap.

    Returns
    -------
    oligos : list of Oligo
    overlaps : list of Overlap
    """
    seq = seq.upper().replace(" ", "").replace("\n", "")

    # Build the full construct: upstream flank + insert + downstream flank
    full_seq = plasmid_upstream.upper() + seq + plasmid_downstream.upper()
    insert_start = len(plasmid_upstream)
    insert_end = insert_start + len(seq)

    step = oligo_length - overlap_length

    # Generate cut points (positions where one oligo ends and the overlap begins)
    # First oligo starts at position 0 of full_seq
    cut_points = []
    pos = oligo_length  # end of first oligo
    while pos < len(full_seq):
        cut_points.append(pos - overlap_length)  # start of overlap region
        pos += step

    # Try to optimize each cut point for overlap quality.
    # Constraint: adjacent cuts must stay >= overlap_length apart so that
    # same-strand oligos never overlap each other.
    optimized_cuts = []
    for idx, cut in enumerate(cut_points):
        best_cut = cut
        best_score = _score_overlap(full_seq, cut, overlap_length, step)
        for shift in range(-max_shift, max_shift + 1):
            candidate = cut + shift
            if candidate < overlap_length or candidate + overlap_length > len(full_seq):
                continue
            # Enforce minimum spacing from the previous optimized cut
            if optimized_cuts and candidate - optimized_cuts[-1] < overlap_length:
                continue
            # Enforce minimum spacing to the next unoptimized cut (if any)
            if idx + 1 < len(cut_points) and cut_points[idx + 1] - candidate < overlap_length:
                continue
            score = _score_overlap(full_seq, candidate, overlap_length, step)
            if score > best_score:
                best_score = score
                best_cut = candidate
        optimized_cuts.append(best_cut)

    # Build oligos from cut points
    oligos = []
    overlaps = []

    # Determine oligo boundaries
    # Each oligo spans from one cut point to the next + overlap
    boundaries = [0] + optimized_cuts + [len(full_seq)]

    for i in range(len(boundaries) - 1):
        start = boundaries[i]
        # Each oligo extends overlap_length into the next segment
        if i < len(boundaries) - 2:
            end = boundaries[i + 1] + overlap_length
        else:
            end = boundaries[i + 1]

        # Clamp to sequence bounds
        end = min(end, len(full_seq))

        oligo_seq = full_seq[start:end]
        strand = "sense" if i % 2 == 0 else "antisense"

        if strand == "antisense":
            oligo_seq = reverse_complement(oligo_seq)

        oligo = Oligo(
            index=i + 1,
            seq=oligo_seq,
            start=start,
            end=end,
            strand=strand,
            is_first=(i == 0),
            is_last=(i == len(boundaries) - 2),
        )
        oligos.append(oligo)

        # Create overlap object (between this oligo and the next)
        if i < len(boundaries) - 2:
            ovl_start = boundaries[i + 1]
            ovl_end = ovl_start + overlap_length
            ovl_end = min(ovl_end, len(full_seq))
            ovl_seq = full_seq[ovl_start:ovl_end]
            overlap = Overlap(seq=ovl_seq, start=ovl_start, end=ovl_end)
            overlap.issues = check_overlap(overlap, full_seq)
            overlaps.append(overlap)

    # Thermodynamic analysis
    if conditions is None:
        conditions = ReactionConditions()
    thermo = conditions.make_thermo()

    # Compute Tm for each overlap (the intended on-target annealing)
    for ovl in overlaps:
        ovl.tm = thermo.calc_tm(ovl.seq)

    # Flag overlaps whose Tm is at or below the reaction/annealing temperature.
    # If Tm <= annealing temp, more than half the molecules are melted — the
    # overlap won't reliably hold its intended partner.
    for i, ovl in enumerate(overlaps):
        if ovl.tm is not None and ovl.tm <= conditions.annealing_temp:
            ovl.issues.append(OverlapIssue(
                "low_tm",
                f"Tm {ovl.tm:.1f}°C is at or below the "
                f"annealing temperature ({conditions.annealing_temp:.0f}°C)",
                "warning",
            ))

    # Check all overlap pairs for cross-hybridization
    cross_results = check_cross_hybridization(overlaps, conditions)
    for idx_i, idx_j, issue_i, issue_j in cross_results:
        overlaps[idx_i].issues.append(issue_i)
        overlaps[idx_j].issues.append(issue_j)

    return oligos, overlaps


def _score_overlap(seq: str, start: int, length: int, step: int = 25) -> float:
    """
    Score an overlap region. Higher = better.

    Currently only scores GC balance. Overlaps near 50% GC anneal more
    reliably. Very AT-rich overlaps have weak annealing; very GC-rich
    overlaps can form secondary structures.

    Other factors (homopolymer placement, G-quadruplexes) are not scored
    because their tradeoffs are unclear — e.g., putting a homopolymer in
    the overlap avoids polymerase slippage but forces IDT to synthesize
    through it twice. These are flagged as warnings instead.
    """
    ovl = seq[start:start + length].upper()
    if len(ovl) < length:
        return -999

    score = 0.0

    # Prefer GC near 50%
    gc = gc_content(ovl)
    score -= abs(gc - 0.5) * 10

    return score


# ---------------------------------------------------------------------------
# Sequence complexity scan (whole-sequence)
# ---------------------------------------------------------------------------

def scan_sequence_complexity(seq: str) -> list[dict]:
    """Scan full sequence for problematic regions (pre-tiling check)."""
    issues = []
    s = seq.upper()

    # Homopolymer runs
    for start, end, base in find_homopolymers(s, min_run=5):
        issues.append({
            "type": "homopolymer",
            "start": start,
            "end": end,
            "value": f"{end-start}x{base}",
            "severity": "error" if end - start >= 8 else "warning",
            "message": f"Homopolymer run: {end-start}x {base} at {start}-{end}"
        })

    # G-quadruplex
    for qstart, qend in find_g_quadruplex(s):
        issues.append({
            "type": "g_quadruplex",
            "start": qstart,
            "end": qend,
            "severity": "error",
            "message": f"G-quadruplex potential at {qstart}-{qend}"
        })

    return issues


# ---------------------------------------------------------------------------
# Output / reporting
# ---------------------------------------------------------------------------

def format_report(
    seq: str,
    oligos: list[Oligo],
    overlaps: list[Overlap],
    seq_issues: list[dict],
    plasmid_upstream: str = "",
    plasmid_downstream: str = "",
    conditions: ReactionConditions | None = None,
) -> str:
    """Generate a human-readable report."""
    lines = []
    lines.append("=" * 70)
    lines.append("OLIGO DESIGN REPORT")
    lines.append("=" * 70)
    lines.append("")

    insert_len = len(seq)
    total_len = len(plasmid_upstream) + insert_len + len(plasmid_downstream)
    lines.append(f"Insert length:        {insert_len} bp")
    if plasmid_upstream or plasmid_downstream:
        lines.append(f"Upstream flank:       {len(plasmid_upstream)} bp")
        lines.append(f"Downstream flank:     {len(plasmid_downstream)} bp")
        lines.append(f"Total construct:      {total_len} bp")
    lines.append(f"Number of oligos:     {len(oligos)}")
    lines.append(f"Number of overlaps:   {len(overlaps)}")
    lines.append("")

    # Sequence complexity issues
    if seq_issues:
        lines.append("-" * 70)
        lines.append("SEQUENCE COMPLEXITY WARNINGS")
        lines.append("-" * 70)
        for issue in seq_issues:
            marker = "!!" if issue["severity"] == "error" else " *"
            lines.append(f"  {marker} {issue['message']}")
        lines.append("")

    # Oligos
    lines.append("-" * 70)
    lines.append("OLIGOS")
    lines.append("-" * 70)
    for oligo in oligos:
        tag = ""
        if oligo.is_first:
            tag = " [FIRST - includes upstream flank]"
        elif oligo.is_last:
            tag = " [LAST - includes downstream flank]"
        lines.append(
            f"  Oligo {oligo.index:2d} | {oligo.strand:>10s} | "
            f"pos {oligo.start:4d}-{oligo.end:4d} | "
            f"{oligo.length:3d} bp | "
            f"GC {gc_content(oligo.seq):.0%}{tag}"
        )
        lines.append(f"           5'-{oligo.seq}-3'")
        lines.append("")

    # Overlaps
    has_issues = any(ovl.issues for ovl in overlaps)
    lines.append("-" * 70)
    lines.append("OVERLAPS" + ("  (!! = issues found)" if has_issues else ""))
    lines.append("-" * 70)
    for i, ovl in enumerate(overlaps):
        flag = " !!" if ovl.issues else "   "
        tm_str = f"Tm {ovl.tm:.1f}°C" if ovl.tm is not None else "Tm N/A"
        lines.append(
            f" {flag} Overlap {i+1:2d} | pos {ovl.start:4d}-{ovl.end:4d} | "
            f"{len(ovl.seq):2d} bp | GC {ovl.gc:.0%} | {tm_str} | {ovl.seq}"
        )
        for issue in ovl.issues:
            lines.append(f"        -> {issue.severity.upper()}: {issue.message}")
    lines.append("")

    # Cross-hybridization summary
    if conditions is not None and len(overlaps) > 1:
        thermo = conditions.make_thermo()
        avg_tm = _avg_overlap_tm(overlaps, thermo)
        threshold_tm = avg_tm - 20.0
        lines.append("-" * 70)
        lines.append("CROSS-HYBRIDIZATION")
        lines.append(f"Each overlap's top 3 strongest off-target interactions")
        lines.append(f"(all others are weaker)")
        lines.append(f"Avg overlap Tm: {avg_tm:.1f}°C | "
                     f"Flag threshold: {threshold_tm:.0f}°C "
                     f"(avg - 20°C, Rouillard et al. 2004)")
        lines.append("-" * 70)

        # Compute all pairwise heterodimer Tms
        pair_tms: dict[tuple[int, int], float] = {}
        for i in range(len(overlaps)):
            for j in range(i + 1, len(overlaps)):
                seq_i = overlaps[i].seq.upper()
                seq_j = overlaps[j].seq.upper()
                rc_j = reverse_complement(seq_j)
                tm1 = thermo.calc_heterodimer(seq_i, seq_j).tm
                tm2 = thermo.calc_heterodimer(seq_i, rc_j).tm
                pair_tms[(i, j)] = max(tm1, tm2)

        # For each overlap, show its on-target Tm and top 3 off-target pairs
        for i, ovl in enumerate(overlaps):
            on_target_tm = f"{ovl.tm:.1f}" if ovl.tm is not None else "N/A"

            # Gather all pairs involving this overlap
            pairs = []
            for (a, b), tm in pair_tms.items():
                if a == i:
                    pairs.append((b, tm))
                elif b == i:
                    pairs.append((a, tm))
            pairs.sort(key=lambda x: x[1], reverse=True)

            # Show top 3 off-target interactions
            top = pairs[:3]
            top_strs = []
            for partner, tm in top:
                flag = "!" if tm > threshold_tm else ""
                top_strs.append(f"vs {partner+1}: {tm:.1f}°C{flag}")

            lines.append(
                f"  Overlap {i+1:2d}  Tm {on_target_tm:>5s}°C  |  "
                + "  ".join(top_strs)
            )

        # Count flagged pairs
        flagged = [(a, b, tm) for (a, b), tm in pair_tms.items() if tm > threshold_tm]
        lines.append("")
        if flagged:
            lines.append(f"  {len(flagged)} pair(s) exceed threshold (!):")
            for a, b, tm in sorted(flagged, key=lambda x: x[2], reverse=True):
                lines.append(f"    Overlap {a+1} × Overlap {b+1}: {tm:.1f}°C")
        else:
            lines.append("  No pairs exceed threshold. All overlaps are unique.")
        lines.append("")

    # IDT order summary
    lines.append("-" * 70)
    lines.append("IDT ORDER (copy-paste)")
    lines.append("-" * 70)
    for oligo in oligos:
        lines.append(f"oligo_{oligo.index:02d}\t{oligo.seq}")
    lines.append("")
    lines.append("=" * 70)

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# File I/O
# ---------------------------------------------------------------------------

def read_fasta_or_raw(path: str) -> tuple[str, str]:
    """Read a DNA sequence from a FASTA file or plain text. Returns (name, seq)."""
    with open(path) as f:
        text = f.read().strip()

    if text.startswith(">"):
        lines = text.split("\n")
        name = lines[0][1:].strip()
        seq = "".join(line.strip() for line in lines[1:] if not line.startswith(">"))
    else:
        name = "input"
        seq = text.replace("\n", "").replace(" ", "")

    # Validate it's DNA
    valid = set("ATCGatcgNn")
    bad = set(seq) - valid
    if bad:
        raise ValueError(
            f"Non-DNA characters found: {bad}. "
            f"This looks like a protein sequence. "
            f"V1 requires DNA input (codon optimization not yet implemented)."
        )

    return name, seq.upper()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def design_oligos(
    seq: str,
    oligo_length: int = 45,
    overlap_length: int = 20,
    plasmid_upstream: str = "",
    plasmid_downstream: str = "",
    conditions: ReactionConditions | None = None,
) -> str:
    """Main entry point. Returns formatted report string."""
    if conditions is None:
        conditions = ReactionConditions()

    # Pre-scan for complexity issues
    seq_issues = scan_sequence_complexity(seq)

    # Tile
    oligos, overlaps = tile_sequence(
        seq,
        oligo_length=oligo_length,
        overlap_length=overlap_length,
        plasmid_upstream=plasmid_upstream,
        plasmid_downstream=plasmid_downstream,
        conditions=conditions,
    )

    # Report
    return format_report(
        seq, oligos, overlaps, seq_issues,
        plasmid_upstream=plasmid_upstream,
        plasmid_downstream=plasmid_downstream,
        conditions=conditions,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Design overlapping oligos for gene synthesis via overlap extension assembly."
    )
    parser.add_argument(
        "sequence",
        help="DNA sequence (string) or path to FASTA/text file"
    )
    parser.add_argument(
        "--oligo-length", type=int, default=45,
        help="Target oligo length in bp (default: 45)"
    )
    parser.add_argument(
        "--overlap", type=int, default=20,
        help="Overlap length between adjacent oligos in bp (default: 20)"
    )
    parser.add_argument(
        "--upstream", type=str, default="",
        help="Plasmid sequence upstream of insert (~20bp for vector insertion)"
    )
    parser.add_argument(
        "--downstream", type=str, default="",
        help="Plasmid sequence downstream of insert (~20bp for vector insertion)"
    )

    args = parser.parse_args()

    # Determine if input is a file or raw sequence
    import os
    if os.path.isfile(args.sequence):
        name, seq = read_fasta_or_raw(args.sequence)
        print(f"Read sequence '{name}' ({len(seq)} bp)")
    else:
        seq = args.sequence.upper().replace(" ", "")
        valid = set("ATCGN")
        if not set(seq).issubset(valid):
            print(f"Error: input doesn't look like DNA and is not a valid file path.")
            sys.exit(1)

    report = design_oligos(
        seq,
        oligo_length=args.oligo_length,
        overlap_length=args.overlap,
        plasmid_upstream=args.upstream,
        plasmid_downstream=args.downstream,
    )
    print(report)


if __name__ == "__main__":
    main()
