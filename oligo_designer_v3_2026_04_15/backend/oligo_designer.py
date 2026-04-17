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
            tm_method='santalucia',              # SantaLucia 1998 nearest-neighbor model
            salt_correction_method='owczarzy',    # Owczarzy 2008 Mg2+ correction
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


@dataclass
class GBlockFragment:
    """A double-stranded gBlock fragment (not tiled into oligos)."""
    index: int
    seq: str             # sense strand, 5'->3'
    start: int           # position in the full construct
    end: int
    label: str = ""
    length: int = 0

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


def find_g_quadruplex(seq: str) -> list[tuple[int, int, int, str]]:
    """Find potential G-quadruplex motifs: ≥4 runs of ≥3 Gs within 30 bp.

    C-runs on the sense strand form a G-quadruplex on the antisense strand,
    so both orientations are checked.  Overlapping candidates are
    consolidated into maximal clusters — one reported entry per cluster.

    Returns a list of (start, end, n_runs, strand) tuples where
    strand ∈ {"sense", "antisense"} and n_runs is the number of distinct
    G/C runs inside the reported span.
    """
    s = seq.upper()
    quads: list[tuple[int, int, int, str]] = []
    for pattern, strand in ((r"G{3,}", "sense"), (r"C{3,}", "antisense")):
        runs = [(m.start(), m.end()) for m in re.finditer(pattern, s)]
        i = 0
        while i < len(runs):
            # Greedy extension: grow the cluster as long as the 30-bp window
            # rule (from the first run's start) still holds for the next run.
            j = i
            while j + 1 < len(runs) and runs[j + 1][1] - runs[i][0] <= 30:
                j += 1
            count = j - i + 1
            if count >= 4:
                quads.append((runs[i][0], runs[j][1], count, strand))
                i = j + 1  # consumed this cluster; skip its runs
            else:
                i += 1
    return sorted(quads)


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

    # G-quadruplex — one entry per consolidated cluster, with run count +
    # strand so the user can tell at a glance how severe the motif is.
    for q_start, q_end, n_runs, strand in find_g_quadruplex(seq):
        letter = "G" if strand == "sense" else "C"
        strand_suffix = "" if strand == "sense" else " (antisense)"
        issues.append(OverlapIssue(
            "g_quadruplex",
            f"G-quadruplex: {n_runs} {letter}-runs at "
            f"{overlap.start + q_start}-{overlap.start + q_end}{strand_suffix}",
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
    # Score by Tm: target = annealing_temp + 12°C (overlaps should anneal
    # reliably above the reaction temperature with minimal Tm variance).
    if conditions is None:
        conditions = ReactionConditions()
    _thermo = conditions.make_thermo()
    _target_tm = conditions.annealing_temp + 12.0

    # Greedy left-to-right cut placement with cross-hyb awareness.
    # `placed_kmers` accumulates k-mers from overlaps we've already chosen;
    # each subsequent cut is scored against this set.
    #
    # We also enforce a HARD MAX of `oligo_length` per oligo by bounding
    # each cut against the previous one. The relevant equations:
    #   oligo[0]      length = cut[0] + overlap_length            <= oligo_length
    #   oligo[i]      length = (cut[i] - cut[i-1]) + overlap_length <= oligo_length
    #   oligo[N]      length = len(full_seq) - cut[N-1]           <= oligo_length
    # Translating each to a position bound on the cut being placed:
    #   first cut:  cut[0] in [overlap_length, oligo_length - overlap_length]
    #   middle:     cut[i] in [cut[i-1] + overlap_length, cut[i-1] + step]
    #   last cut:   cut[N-1] in [max(prev+overlap, len_seq - oligo_length),
    #                            prev + step]
    # This may make the search asymmetric (e.g. only LEFT shifts allowed
    # for the first cut), but it guarantees no oligo exceeds the user's
    # length cap — which is the IDT pricing-tier constraint.
    placed_kmers: set[str] = set()
    optimized_cuts: list[int] = []
    n_cuts = len(cut_points)
    for idx, cut in enumerate(cut_points):
        prev_bound = optimized_cuts[-1] if optimized_cuts else -(overlap_length * 2)
        next_bound = cut_points[idx + 1] if idx + 1 < n_cuts else len(full_seq) + overlap_length

        # Per-cut hard bounds enforcing oligo_length on the oligos this cut
        # closes (the one to its LEFT) and, for the last cut, the trailing
        # oligo to its RIGHT.
        if idx == 0:
            cut_min_i = overlap_length
            cut_max_i = oligo_length - overlap_length
        else:
            cut_min_i = optimized_cuts[-1] + overlap_length
            cut_max_i = optimized_cuts[-1] + step
        if idx == n_cuts - 1:
            # Last cut also caps the trailing oligo: len_seq - cut <= oligo_length
            cut_min_i = max(cut_min_i, len(full_seq) - oligo_length)

        best_cut, _best_score, _n_conflicts = _optimize_cut(
            full_seq, cut, overlap_length, step,
            thermo=_thermo, target_tm=_target_tm,
            placed_kmers=placed_kmers,
            prev_cut_bound=prev_bound,
            next_cut_bound=next_bound,
            seq_min=overlap_length,
            seq_max=len(full_seq),
            base_shift=max_shift,
            cut_min=cut_min_i,
            cut_max=cut_max_i,
        )
        optimized_cuts.append(best_cut)
        placed_kmers |= _kmer_set_bidirectional(full_seq[best_cut:best_cut + overlap_length])

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


def tile_sequence_with_gblocks(
    seq: str,
    gblock_regions: list[tuple[int, int, str]],
    oligo_length: int = 60,
    overlap_length: int = 20,
    plasmid_upstream: str = "",
    plasmid_downstream: str = "",
    max_shift: int = 5,
    conditions: ReactionConditions | None = None,
) -> tuple[list[Oligo], list[Overlap], list[GBlockFragment]]:
    """Tile a DNA sequence into oligos + gBlock fragments.

    gBlocks are treated exactly like oligos: they're segments in a unified
    tiling with overlaps at every junction. The only difference is that
    gBlock boundaries are fixed (not placed by the tiler) and gBlocks are
    emitted as GBlockFragment instead of Oligo.

    Algorithm:
      1. Build ordered segment list: oligo-region, gBlock, oligo-region, ...
      2. Tile each oligo-region into individual oligos using step-based cuts
      3. Each segment extends overlap_length into its neighbors
      4. Compute overlaps at every junction (oligo-oligo, oligo-gBlock)
    """
    seq = seq.upper().replace(" ", "").replace("\n", "")
    if conditions is None:
        conditions = ReactionConditions()

    full_seq = plasmid_upstream.upper() + seq + plasmid_downstream.upper()
    insert_start = len(plasmid_upstream)
    step = oligo_length - overlap_length
    thermo = conditions.make_thermo()
    target_tm = conditions.annealing_temp + 12.0

    # Running k-mer set used throughout this assembly to avoid cross-hyb
    # between any two placed overlaps (gBlock boundaries and oligo cuts
    # share the same pool).
    placed_kmers: set[str] = set()

    # Convert gBlock positions to full_seq coordinates, optimize boundaries.
    # Each gBlock contributes two boundary overlaps; both go into placed_kmers
    # as soon as they're selected so later boundaries / oligo cuts avoid them.
    gblocks_full = []
    for i, (gs, ge, label) in enumerate(sorted(gblock_regions, key=lambda g: g[0])):
        gs_f, ge_f = gs + insert_start, ge + insert_start
        # Optimize start boundary: score the overlap that starts at the cut.
        if gs_f > 0:
            prev_bound = gblocks_full[-1][1] if gblocks_full else -(overlap_length * 2)
            gs_f, _s, _c = _optimize_cut(
                full_seq, gs_f, overlap_length, step,
                thermo=thermo, target_tm=target_tm,
                placed_kmers=placed_kmers,
                prev_cut_bound=prev_bound,
                next_cut_bound=ge_f,
                seq_min=overlap_length,
                seq_max=ge_f - overlap_length,
                base_shift=max_shift,
            )
            placed_kmers |= _kmer_set_bidirectional(
                full_seq[gs_f:gs_f + overlap_length]
            )
        # Optimize end boundary: the overlap ends AT the cut, so its start is
        # (cut - overlap_length).  _optimize_cut scores overlap starting at
        # the cut, so we translate: treat ge_f - overlap_length as the
        # nominal cut, then translate back.
        if ge_f < len(full_seq):
            end_start_nominal = ge_f - overlap_length
            end_start, _s, _c = _optimize_cut(
                full_seq, end_start_nominal, overlap_length, step,
                thermo=thermo, target_tm=target_tm,
                placed_kmers=placed_kmers,
                prev_cut_bound=gs_f,
                next_cut_bound=len(full_seq) + overlap_length,
                seq_min=gs_f + overlap_length,
                seq_max=len(full_seq),
                base_shift=max_shift,
            )
            ge_f = end_start + overlap_length
            placed_kmers |= _kmer_set_bidirectional(
                full_seq[end_start:ge_f]
            )
        gblocks_full.append((gs_f, ge_f, label))

    # Build ordered segment list in full_seq coordinates:
    #   [(type, start, end, label?), ...]
    segments = []
    cursor = 0
    for gs, ge, label in gblocks_full:
        if cursor < gs:
            segments.append(("oligo_region", cursor, gs))
        segments.append(("gblock", gs, ge, label))
        cursor = ge
    if cursor < len(full_seq):
        segments.append(("oligo_region", cursor, len(full_seq)))

    # Expand each segment into individual fragments. Each fragment is one
    # "piece" in the assembly: either a single oligo or a gBlock.
    # Every fragment extends overlap_length past the next fragment's start.
    fragments = []  # [(type, start, end, label)]  — in full_seq coords

    for seg_idx, seg in enumerate(segments):
        if seg[0] == "gblock":
            _, gs, ge, label = seg
            fragments.append(("gblock", gs, ge, label))
            continue

        _, rs, re = seg
        region_len = re - rs
        if region_len <= 0:
            continue

        # If this oligo region abuts a gBlock, its boundary oligo will
        # later be extended by overlap_length INTO the gBlock to create
        # the Gibson overlap. That extension has to be counted here when
        # we decide how many cuts to place — otherwise the boundary
        # oligo can blow past oligo_length.
        bord_L = seg_idx > 0 and segments[seg_idx - 1][0] == "gblock"
        bord_R = seg_idx < len(segments) - 1 and segments[seg_idx + 1][0] == "gblock"

        # Effective endpoints: where the first and last oligos actually
        # start / end after gBlock extensions.
        rs_eff = rs - overlap_length if bord_L else rs
        re_eff = re + overlap_length if bord_R else re
        eff_span = re_eff - rs_eff

        # Minimum number of oligos to tile eff_span with each oligo <= oligo_length:
        #   N*oligo_length - (N-1)*overlap >= eff_span
        #   N >= (eff_span - overlap) / step
        if eff_span <= oligo_length:
            n_oligos = 1
        else:
            n_oligos = (eff_span - overlap_length + step - 1) // step  # ceil
            n_oligos = max(1, n_oligos)
        n_region_cuts = n_oligos - 1

        # Nominal cut positions: same tiling rule as tile_sequence, just
        # anchored to the effective endpoints. The while-loop would
        # over-count on short regions where we manually set n_region_cuts,
        # so we just generate exactly n_region_cuts evenly-spaced cuts.
        cuts: list[int] = []
        if n_region_cuts > 0:
            pos = rs_eff + oligo_length
            for _ in range(n_region_cuts):
                cuts.append(pos - overlap_length)
                pos += step

        # Cut placement. Bounds enforce per-oligo length cap, properly
        # accounting for the gBlock extensions at region boundaries.
        #   oligo_0 length = cut[0] + overlap - rs_eff         <= oligo_length
        #       => cut[0] <= rs_eff + step
        #   oligo_i length = cut[i] - cut[i-1] + overlap        <= oligo_length
        #       => cut[i] <= cut[i-1] + step
        #   oligo_N length = re_eff - cut[N-1]                  <= oligo_length
        #       => cut[N-1] >= re_eff - oligo_length
        # Physical constraint: cuts must stay inside [rs, re - overlap]
        # so the outgoing overlap [cut, cut+overlap] fits inside the
        # physical region. cut == rs is allowed — it corresponds to a
        # valid "bridge oligo" pattern where the boundary oligo is just
        # the gBlock-extension plus the outgoing overlap (length =
        # 2*overlap_length, no unique middle content). This is the
        # minimum boundary oligo size and it's what lets small
        # oligo_length values (e.g. 50bp with 20bp overlaps) work.
        opt_cuts: list[int] = []
        for idx, cut in enumerate(cuts):
            prev_bound = opt_cuts[-1] if opt_cuts else rs - overlap_length
            next_bound = cuts[idx + 1] if idx + 1 < n_region_cuts else re + overlap_length

            if idx == 0:
                cut_min_i = rs
                cut_max_i = rs_eff + step
            else:
                cut_min_i = opt_cuts[-1] + overlap_length
                cut_max_i = opt_cuts[-1] + step
            if idx == n_region_cuts - 1:
                cut_min_i = max(cut_min_i, re_eff - oligo_length)

            best_cut, _s, _c = _optimize_cut(
                full_seq, cut, overlap_length, step,
                thermo=thermo, target_tm=target_tm,
                placed_kmers=placed_kmers,
                prev_cut_bound=prev_bound,
                next_cut_bound=next_bound,
                seq_min=rs,
                seq_max=re,
                base_shift=max_shift,
                cut_min=cut_min_i,
                cut_max=cut_max_i,
            )
            opt_cuts.append(best_cut)
            placed_kmers |= _kmer_set_bidirectional(
                full_seq[best_cut:best_cut + overlap_length]
            )

        # Build oligo fragments from boundaries within this region
        bounds = [rs] + opt_cuts + [re]
        for i in range(len(bounds) - 1):
            frag_start = bounds[i]
            frag_end = bounds[i + 1] + overlap_length if i < len(bounds) - 2 else bounds[i + 1]
            frag_end = min(frag_end, len(full_seq))
            fragments.append(("oligo", frag_start, frag_end, ""))

    # Sort fragments by start position. gBlocks sort before oligos at the
    # same position so that the extension step correctly identifies neighbors.
    type_order = {"gblock": 0, "oligo": 1}
    fragments.sort(key=lambda f: (f[1], type_order.get(f[0], 2)))

    # Extend boundary oligos into adjacent gBlocks so they overlap for
    # Gibson assembly. This is the key step that treats gBlocks like oligos.
    for i in range(len(fragments)):
        typ, start, end, lbl = fragments[i]
        if typ != "oligo":
            continue
        # Extend forward into next gBlock
        if i + 1 < len(fragments) and fragments[i + 1][0] == "gblock":
            gb_start = fragments[i + 1][1]
            new_end = min(gb_start + overlap_length, fragments[i + 1][2])
            end = max(end, new_end)
        # Extend backward into previous gBlock
        if i > 0 and fragments[i - 1][0] == "gblock":
            gb_end = fragments[i - 1][2]
            new_start = max(gb_end - overlap_length, fragments[i - 1][1])
            start = min(start, new_start)
        fragments[i] = (typ, start, end, lbl)

    # Build final objects
    oligos = []
    overlaps_list = []
    gblock_fragments = []
    oligo_idx = 0

    # All emitted positions (oligo/gBlock/overlap start/end) are in
    # FULL_SEQ coordinates — consistent with tile_sequence() and with
    # what the frontend's fullSequence (upstream + insert + downstream)
    # expects. When plasmid flanks are empty, insert_start == 0 and
    # full_seq coords == insert coords, so legacy behavior is preserved.
    for i, (frag_type, frag_start, frag_end, label) in enumerate(fragments):
        if frag_type == "gblock":
            gblock_fragments.append(GBlockFragment(
                index=len(gblock_fragments),
                seq=full_seq[frag_start:frag_end],
                start=frag_start, end=frag_end, label=label,
            ))
        else:
            frag_seq = full_seq[frag_start:frag_end]
            strand = "sense" if oligo_idx % 2 == 0 else "antisense"
            if strand == "antisense":
                frag_seq = reverse_complement(frag_seq)
            oligos.append(Oligo(
                index=oligo_idx, seq=frag_seq,
                start=frag_start, end=frag_end, strand=strand,
                is_first=(frag_start == 0),
                is_last=(frag_end == len(full_seq)),
            ))
            oligo_idx += 1

        # Overlap with next fragment: the physical intersection of this
        # fragment's range and the next fragment's range.
        if i < len(fragments) - 1:
            next_frag = fragments[i + 1]
            ovl_start = max(frag_start, next_frag[1])
            ovl_end = min(frag_end, next_frag[2])

            if ovl_end > ovl_start and ovl_end - ovl_start >= 6:
                ovl_seq = full_seq[ovl_start:ovl_end]
                ovl_obj = Overlap(
                    seq=ovl_seq,
                    start=ovl_start,
                    end=ovl_end,
                    tm=thermo.calc_tm(ovl_seq),
                )
                ovl_obj.issues = check_overlap(ovl_obj, full_seq)
                overlaps_list.append(ovl_obj)

    # Sort and re-index
    oligos.sort(key=lambda o: o.start)
    for i, o in enumerate(oligos):
        o.index = i
    overlaps_list.sort(key=lambda o: o.start)

    # Cross-hybridization check
    cross_results = check_cross_hybridization(overlaps_list, conditions)
    for idx_i, idx_j, issue_i, issue_j in cross_results:
        overlaps_list[idx_i].issues.append(issue_i)
        overlaps_list[idx_j].issues.append(issue_j)

    return oligos, overlaps_list, gblock_fragments


# ---------------------------------------------------------------------------
# Cross-hybridization avoidance during cut placement
#
# Rationale: the existing tiler scored each cut purely by its own Tm match,
# and shifted only ±max_shift (default 5 bp) from each nominal cut.  In
# linker-heavy constructs (GGGGS repeats) this is catastrophic — every
# overlap that lands inside a linker shares its 10-mer spelling with every
# other such overlap, producing cross-hybridization flags that the shifter
# can't escape.
#
# We fix this by:
#   (a) maintaining a running set of k-mers (from both strands) contributed
#       by already-placed overlaps, and
#   (b) penalizing any candidate overlap that shares k-mers with that set.
# If ±max_shift produces no conflict-free candidate, we widen the search to
# the largest span the spacing constraints allow and try again.  The k-mer
# penalty guides the scorer toward a unique window without adding any
# thermodynamic cost in the inner loop.
# ---------------------------------------------------------------------------

KMER_SIZE_FOR_XHYB = 10
KMER_PENALTY = 5.0  # per unique matching k-mer — dominates ~2 °C of Tm drift


def _kmer_set_bidirectional(seq: str, k: int = KMER_SIZE_FOR_XHYB) -> set[str]:
    """k-mers of seq plus k-mers of its reverse complement.

    A placed overlap can cross-hybridize with another overlap (sense-sense)
    or with another overlap's reverse complement (sense-antisense), so we
    store both orientations and later check one strand of each candidate.
    """
    s = seq.upper()
    rc = reverse_complement(s)
    kmers: set[str] = set()
    if len(s) >= k:
        for i in range(len(s) - k + 1):
            kmers.add(s[i:i + k])
    if len(rc) >= k:
        for i in range(len(rc) - k + 1):
            kmers.add(rc[i:i + k])
    return kmers


def _count_kmer_conflicts(
    seq: str, placed: set[str], k: int = KMER_SIZE_FOR_XHYB
) -> int:
    """Number of distinct k-mers in seq that already appear in `placed`."""
    if not placed or len(seq) < k:
        return 0
    s = seq.upper()
    matches = set()
    for i in range(len(s) - k + 1):
        kmer = s[i:i + k]
        if kmer in placed:
            matches.add(kmer)
    return len(matches)


def _score_overlap(
    seq: str, start: int, length: int, step: int = 25,
    thermo: thermoanalysis.ThermoAnalysis | None = None,
    target_tm: float | None = None,
    placed_kmers: set[str] | None = None,
) -> float:
    """
    Score an overlap region. Higher = better.

    When a thermo object and target_tm are provided, scores by how close
    the overlap's Tm is to the target (minimizes Tm variance across the
    construct). Falls back to GC balance scoring if thermo is not available.

    If `placed_kmers` is supplied, subtracts a penalty for each k-mer this
    candidate shares with previously-placed overlaps (cross-hyb avoidance).
    """
    ovl = seq[start:start + length].upper()
    if len(ovl) < length:
        return -999

    if thermo is not None and target_tm is not None:
        tm = thermo.calc_tm(ovl)
        score = -abs(tm - target_tm)
    else:
        gc = gc_content(ovl)
        score = -abs(gc - 0.5) * 10

    if placed_kmers:
        score -= _count_kmer_conflicts(ovl, placed_kmers) * KMER_PENALTY

    return score


def _optimize_cut(
    full_seq: str,
    nominal_cut: int,
    overlap_length: int,
    step: int,
    *,
    thermo: thermoanalysis.ThermoAnalysis,
    target_tm: float,
    placed_kmers: set[str],
    prev_cut_bound: int,     # no candidate may land closer than overlap_length
    next_cut_bound: int,     # no candidate may land closer than overlap_length
    seq_min: int,            # hard lower bound (usually = overlap_length)
    seq_max: int,            # hard upper bound (usually = len(full_seq))
    base_shift: int = 5,
    widen_cap: int | None = None,
    cut_min: int | None = None,  # hard lower bound on the cut position itself
    cut_max: int | None = None,  # hard upper bound on the cut position itself
) -> tuple[int, float, int]:
    """Pick the best cut position around `nominal_cut`.

    First searches ±base_shift; if the winner still has k-mer conflicts with
    previously-placed overlaps, widens the search up to `widen_cap` (or the
    largest span the spacing constraints permit) to try to find a clean
    window.

    `cut_min` and `cut_max`, if provided, are HARD bounds on the cut
    position. They are used by the caller to enforce per-oligo length
    limits (e.g. cut[i] <= cut[i-1] + step keeps oligo i within
    `oligo_length`). Any candidate outside [cut_min, cut_max] is rejected.

    Returns (best_cut, best_score, n_kmer_conflicts_at_best).
    """
    # Effective hard bounds combining all sources.
    hard_min = max(seq_min, prev_cut_bound + overlap_length)
    if cut_min is not None:
        hard_min = max(hard_min, cut_min)
    hard_max = min(seq_max - overlap_length, next_cut_bound - overlap_length)
    if cut_max is not None:
        hard_max = min(hard_max, cut_max)

    # If bounds are infeasible, fall back to the closest feasible position.
    # This can happen when the caller's nominal placement is too tight; we
    # prefer "shorter than ideal" over "infeasible" so the assembly still
    # tiles the whole sequence.
    if hard_min > hard_max:
        clipped = max(seq_min, min(seq_max - overlap_length, nominal_cut))
        return clipped, -999.0, 0

    # Max symmetric shift that still satisfies spacing constraints.
    max_feasible_shift = max(
        0,
        min(nominal_cut - hard_min, hard_max - nominal_cut),
    )
    if widen_cap is None:
        widen_cap = max_feasible_shift
    else:
        widen_cap = min(widen_cap, max_feasible_shift)

    def _search(shift: int) -> tuple[int, float, int]:
        # Start with no candidate. We'll accept the first valid one we see
        # and improve from there. (Initializing with `nominal_cut` directly
        # was buggy — if nominal violated cut_min/cut_max we'd return an
        # out-of-bounds cut.)
        best_c: int | None = None
        best_s: float = float("-inf")
        best_n: int = 0
        for s in range(-shift, shift + 1):
            c = nominal_cut + s
            if c < hard_min or c > hard_max:
                continue
            score = _score_overlap(
                full_seq, c, overlap_length, step,
                thermo=thermo, target_tm=target_tm, placed_kmers=placed_kmers,
            )
            if score > best_s:
                best_s = score
                best_c = c
                best_n = _count_kmer_conflicts(
                    full_seq[c:c + overlap_length], placed_kmers,
                )
        if best_c is None:
            # No candidate in the search window — fall back to the closest
            # feasible position to nominal so we still produce a tiling.
            best_c = max(hard_min, min(hard_max, nominal_cut))
            best_s = _score_overlap(
                full_seq, best_c, overlap_length, step,
                thermo=thermo, target_tm=target_tm, placed_kmers=placed_kmers,
            )
            best_n = _count_kmer_conflicts(
                full_seq[best_c:best_c + overlap_length], placed_kmers,
            )
        return best_c, best_s, best_n

    best_cut, best_score, conflicts = _search(min(base_shift, widen_cap))
    # Widen only if we still have conflicts and there's room to grow.
    if conflicts > 0 and widen_cap > base_shift:
        wider_cut, wider_score, wider_conflicts = _search(widen_cap)
        if wider_score > best_score:
            best_cut, best_score, conflicts = wider_cut, wider_score, wider_conflicts

    return best_cut, best_score, conflicts


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
    for qstart, qend, n_runs, strand in find_g_quadruplex(s):
        letter = "G" if strand == "sense" else "C"
        strand_suffix = "" if strand == "sense" else " (antisense)"
        issues.append({
            "type": "g_quadruplex",
            "start": qstart,
            "end": qend,
            "severity": "error",
            "message": f"G-quadruplex: {n_runs} {letter}-runs at "
                       f"{qstart}-{qend}{strand_suffix}"
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
