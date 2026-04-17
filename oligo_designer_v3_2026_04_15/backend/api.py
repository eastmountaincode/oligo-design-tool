"""
FastAPI backend for oligo designer.
Can run standalone or as a Vercel serverless function.
"""

import asyncio
import json
import queue
import threading
import time
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from Bio.Seq import Seq
import re
import sys
import os
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

# Built-in custom codon tables shipped with the project
_BACKEND_DIR = Path(__file__).resolve().parent
CUSTOM_TABLES_DIR = _BACKEND_DIR.parent.parent.parent / "data"
CUSTOM_TABLES: dict[str, str] = {}
if CUSTOM_TABLES_DIR.exists():
    for f in CUSTOM_TABLES_DIR.glob("*CodonTable*.txt"):
        if f.name.startswith("."):
            continue
        CUSTOM_TABLES[f.stem] = f.read_text()
from oligo_designer import (
    tile_sequence, tile_sequence_with_gblocks, scan_sequence_complexity,
    gc_content, design_oligos, ReactionConditions, GBlockFragment,
)

app = FastAPI(title="Oligo Designer", version="1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # tighten in production
    allow_methods=["*"],
    allow_headers=["*"],
)


# ---------------------------------------------------------------------------
# Oligo design endpoint
# ---------------------------------------------------------------------------

class DesignRequest(BaseModel):
    sequence: str = Field(..., description="DNA sequence (ATCGs only)")
    oligo_length: int = Field(60, ge=30, le=200)
    overlap_length: int = Field(20, ge=15, le=40)
    plasmid_upstream: str = Field("", description="~20bp upstream flank for vector insertion")
    plasmid_downstream: str = Field("", description="~20bp downstream flank for vector insertion")
    mv_conc: float = Field(50.0, description="Monovalent cation concentration [Na+/K+] in mM")
    dv_conc: float = Field(2.0, description="Divalent cation concentration [Mg2+] in mM")
    dntp_conc: float = Field(0.8, description="Total dNTP concentration in mM")
    dna_conc: float = Field(250.0, description="Oligo concentration in nM")
    annealing_temp: float = Field(50.0, description="Annealing temperature in °C")
    gblock_regions: list[dict] = Field(default_factory=list, description="gBlock regions [{start_bp, end_bp, label}]")


class OligoOut(BaseModel):
    index: int
    seq: str
    start: int
    end: int
    strand: str
    length: int
    is_first: bool
    is_last: bool
    gc: float


class OverlapOut(BaseModel):
    index: int
    seq: str
    start: int
    end: int
    gc: float
    tm: float | None
    issues: list[dict]


class GBlockOut(BaseModel):
    index: int
    seq: str
    start: int
    end: int
    length: int
    label: str
    gc: float
    fragment_type: str = "gblock"


class DesignResponse(BaseModel):
    insert_length: int
    total_length: int
    num_oligos: int
    num_overlaps: int
    avg_overlap_tm: float | None
    oligos: list[OligoOut]
    overlaps: list[OverlapOut]
    gblocks: list[GBlockOut] = []
    sequence_issues: list[dict]
    report_text: str


@app.post("/api/design", response_model=DesignResponse)
def design(req: DesignRequest):
    seq = req.sequence.upper().replace(" ", "").replace("\n", "")

    bad = set(seq) - set("ATCGN")
    if bad:
        raise HTTPException(400, f"Non-DNA characters: {bad}")

    conditions = ReactionConditions(
        mv_conc=req.mv_conc,
        dv_conc=req.dv_conc,
        dntp_conc=req.dntp_conc,
        dna_conc=req.dna_conc,
        annealing_temp=req.annealing_temp,
    )

    gbs = sorted(req.gblock_regions, key=lambda g: g.get("start_bp", 0)) if req.gblock_regions else []

    if not gbs:
        oligos, overlaps = tile_sequence(
            seq, oligo_length=req.oligo_length, overlap_length=req.overlap_length,
            plasmid_upstream=req.plasmid_upstream, plasmid_downstream=req.plasmid_downstream,
            conditions=conditions,
        )
        report = design_oligos(
            seq, oligo_length=req.oligo_length, overlap_length=req.overlap_length,
            plasmid_upstream=req.plasmid_upstream, plasmid_downstream=req.plasmid_downstream,
            conditions=conditions,
        )
        gblock_frags: list[GBlockFragment] = []
    else:
        gb_tuples = [(g.get("start_bp", 0), g.get("end_bp", 0), g.get("label", ""))
                     for g in gbs]
        oligos, overlaps, gblock_frags = tile_sequence_with_gblocks(
            seq,
            gblock_regions=gb_tuples,
            oligo_length=req.oligo_length,
            overlap_length=req.overlap_length,
            plasmid_upstream=req.plasmid_upstream,
            plasmid_downstream=req.plasmid_downstream,
            conditions=conditions,
        )
        report = f"Designed {len(oligos)} oligos + {len(gblock_frags)} gBlock(s)"

    seq_issues = scan_sequence_complexity(seq)

    final_oligos = [
        OligoOut(
            index=o.index, seq=o.seq, start=o.start, end=o.end,
            strand=o.strand, length=o.length,
            is_first=o.is_first, is_last=o.is_last,
            gc=round(gc_content(o.seq), 3),
        )
        for o in oligos
    ]
    final_overlaps = [
        OverlapOut(
            index=i + 1, seq=ovl.seq, start=ovl.start, end=ovl.end,
            gc=round(ovl.gc, 3),
            tm=round(ovl.tm, 1) if ovl.tm is not None else None,
            issues=[{"kind": iss.kind, "message": iss.message, "severity": iss.severity,
                     "start": iss.start, "end": iss.end}
                    for iss in ovl.issues],
        )
        for i, ovl in enumerate(overlaps)
    ]
    gblock_outs = [
        GBlockOut(
            index=gf.index, seq=gf.seq, start=gf.start, end=gf.end,
            length=gf.length, label=gf.label,
            gc=round(gc_content(gf.seq), 3),
        )
        for gf in gblock_frags
    ]

    overlap_tms = [ovl.tm for ovl in final_overlaps if ovl.tm is not None]
    avg_tm = sum(overlap_tms) / len(overlap_tms) if overlap_tms else None

    return DesignResponse(
        insert_length=len(seq),
        total_length=len(req.plasmid_upstream) + len(seq) + len(req.plasmid_downstream),
        num_oligos=len(final_oligos),
        num_overlaps=len(final_overlaps),
        avg_overlap_tm=round(avg_tm, 1) if avg_tm is not None else None,
        oligos=final_oligos,
        overlaps=final_overlaps,
        gblocks=gblock_outs,
        sequence_issues=seq_issues,
        report_text=report,
    )


# ---------------------------------------------------------------------------
# Codon optimization endpoint
# ---------------------------------------------------------------------------

THREE_TO_ONE = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'End': '*',
}


def parse_codon_table(text: str) -> dict:
    """Parse Kazusa/GCG format codon table into {aa_1letter: {codon: frequency}}.

    Frequencies are normalized to sum to 1.0 per amino acid.
    Uses the /1000 column for relative weights.
    """
    raw = parse_codon_table_raw(text)
    table: dict[str, dict[str, float]] = {}
    for aa, codons in raw.items():
        total = sum(codons.values())
        if total > 0:
            table[aa] = {c: v / total for c, v in codons.items()}
        else:
            table[aa] = codons
    return table


def parse_codon_table_raw(text: str) -> dict[str, dict[str, float]]:
    """Parse Kazusa/GCG format codon table, preserving raw /1000 values."""
    table: dict[str, dict[str, float]] = {}
    for line in text.strip().splitlines():
        parts = line.strip().split()
        if len(parts) < 4:
            continue
        aa_3, codon = parts[0], parts[1].upper()
        if len(codon) != 3 or not all(c in 'ATCG' for c in codon):
            continue
        aa_1 = THREE_TO_ONE.get(aa_3, None)
        if aa_1 is None or aa_1 == '*':
            continue
        per_thousand = float(parts[3])
        if aa_1 not in table:
            table[aa_1] = {}
        table[aa_1][codon] = per_thousand
    return table


class GBlockRegionInput(BaseModel):
    start_aa: int = Field(..., ge=0, description="Start residue (0-indexed, inclusive)")
    end_aa: int = Field(..., ge=0, description="End residue (0-indexed, inclusive)")
    label: str = Field("", description="Optional label for this gBlock region")


# IDT gBlock synthesis constraints (from idtdna.com FAQ, April 2026)
GBLOCK_MIN_BP = 125
GBLOCK_MAX_BP = 3000
GBLOCK_HOMO_GC = 6   # IDT: avoid runs of 6+ G or C
GBLOCK_HOMO_AT = 10   # IDT: avoid runs of 10+ A or T
GBLOCK_GC_MIN = 0.25
GBLOCK_GC_MAX = 0.75


class CodonOptRequest(BaseModel):
    protein_sequence: str = Field(..., description="Amino acid sequence (single-letter codes)")
    species: str = Field("h_sapiens_9606", description="Species for codon optimization")
    custom_codon_table: str | None = Field(None, description="Custom codon table in Kazusa/GCG format")
    avoid_homopolymers_gc: int = Field(4, ge=3, le=8, description="Avoid runs of N+ identical G or C bases")
    avoid_homopolymers_at: int = Field(6, ge=3, le=8, description="Avoid runs of N+ identical A or T bases")
    gc_min: float = Field(0.25, ge=0.0, le=1.0)
    gc_max: float = Field(0.75, ge=0.0, le=1.0)
    gc_window: int = Field(50, ge=20, le=200)
    uniquify_kmers: int = Field(10, ge=6, le=20, description="Max repeated k-mer length to eliminate")
    avoid_patterns: list[str] = Field(default_factory=list, description="Additional patterns to avoid (regex)")
    min_codon_frequency: float = Field(10.0, ge=0.0, le=100.0, description="Minimum codon frequency (/1000) allowed")
    beam_k: int = Field(250, ge=1, le=500, description="Beam search width per bucket. Higher = better quality, slower.")
    gblock_regions: list[GBlockRegionInput] = Field(default_factory=list, description="Regions to be ordered as IDT gBlocks (relaxed synthesis constraints)")


class CodonDetail(BaseModel):
    amino_acid: str
    codon: str
    per_thousand: float


class CodonTableEntry(BaseModel):
    amino_acid: str
    codon: str
    per_thousand: float


class SequenceWarning(BaseModel):
    kind: str
    message: str
    start: int
    end: int
    group: int | None = None  # For repeat k-mers: matching warnings share the same group number


class GBlockSegmentInfo(BaseModel):
    start_aa: int
    end_aa: int
    start_bp: int
    end_bp: int
    label: str = ""
    len_bp: int
    valid_size: bool
    gc: float
    max_g_run: int
    max_c_run: int


class CodonOptResponse(BaseModel):
    input_protein: str
    optimized_dna: str
    back_translated_protein: str
    length_dna: int
    length_protein: int
    gc_content: float
    gc_window_min: float | None = None
    gc_window_max: float | None = None
    gc_window_size: int = 50
    constraints_pass: bool
    constraints_summary: list[dict]
    warnings: list[SequenceWarning]
    codons: list[CodonDetail]
    codon_table: list[CodonTableEntry]
    codon_table_name: str
    optimization_report: dict = {}
    gblock_segments: list[GBlockSegmentInfo] = []
    suggested_gblocks: list[GBlockSegmentInfo] = []   # auto-detected when beam fails
    gblocks_auto_applied: bool = False                 # True if suggestions were used to produce this result


def _find_repetitive_aa_regions(
    protein: str, k: int = 5, min_hits: int = 2,
) -> list[tuple[int, int]]:
    """Find contiguous residue ranges whose k-aa subsequences appear more
    than once in the protein.

    This is the general mechanism behind the "linker problem": any k-aa
    motif that repeats forces the naive reverse translation to repeat at
    the DNA level too, which (a) breaks dnachisel's k-mer uniqueness
    constraint and (b) creates cross-hybridization between distant oligo
    overlaps that land inside copies of the motif.  Detecting repeats at
    the AA level is sequence-agnostic — it catches canonical GGGGS
    linkers, custom linkers, tandem repeats, and homologous framework
    regions alike.

    k = 5 aa ≈ 15 bp, which is comfortably above dnachisel's default
    10-bp uniqueness window — so any flagged region really is one the
    optimizer cannot untangle with codon swaps alone.

    Returns a list of inclusive [start_aa, end_aa] ranges, merged.
    """
    n = len(protein)
    if n < k:
        return []
    from collections import defaultdict
    positions: dict[str, list[int]] = defaultdict(list)
    for i in range(n - k + 1):
        positions[protein[i:i + k]].append(i)
    mask = [False] * n
    for kmer, pos_list in positions.items():
        if len(pos_list) >= min_hits:
            for p in pos_list:
                for j in range(p, min(p + k, n)):
                    mask[j] = True
    ranges: list[tuple[int, int]] = []
    i = 0
    while i < n:
        if mask[i]:
            j = i
            while j < n and mask[j]:
                j += 1
            ranges.append((i, j - 1))
            i = j
        else:
            i += 1
    return ranges


def _auto_detect_gblocks(
    protein: str,
    raw_table: dict,
    *,
    gc_min: float, gc_max: float, gc_window: int,
    homo_gc: int, homo_at: int,
    avoid_patterns: list[str],
    min_codon_frequency: float,
    beam_k: int,
    max_iterations: int = 20,
    grow_step: int = 10,
    initial_radius: int = 5,
    min_aa: int | None = None,
    progress_callback=None,
) -> list[GBlockRegionInput]:
    """Iteratively grow gBlock regions around beam failure points until the
    per-segment beam succeeds.

    Flow:
      1. Seed regions from GGGGS linkers (≥2 repeats) with a small flanking
         buffer — linkers are the most common beam-failure cause AND they
         cross-hybridize badly when tiled, so they belong in gBlocks anyway.
      2. Run the per-segment beam; on failure, grow the nearest region (or
         create a new one) and retry.
      3. After convergence, SHRINK each region as much as possible while
         still letting the beam pass — users pay per gBlock base pair.
      4. Pad each surviving region to the IDT minimum length.

    Returns a list of GBlockRegionInput (padded to IDT minimum size) that,
    when used as gBlock regions, allow the remaining oligo-tiled segments to
    pass the beam under strict constraints.
    """
    from gc_refinement import refine_gc_windows
    from dnachisel import reverse_translate

    if min_aa is None:
        min_aa = -(-GBLOCK_MIN_BP // 3)  # ceil(125/3) = 42

    n = len(protein)

    # Phase 1: seed from repetitive AA regions.  Any k-aa motif that occurs
    # more than once will (a) defeat the beam's uniqueness constraint and
    # (b) cross-hybridize at the oligo layer — both exactly the failures a
    # gBlock is meant to cover.  See _find_repetitive_aa_regions for why
    # this generalizes the "linker problem" to arbitrary sequence.
    repeat_regions = _find_repetitive_aa_regions(protein, k=5, min_hits=2)
    regions: list[list[int]] = []
    for (la, lb) in repeat_regions:
        regions.append([max(0, la - 2), min(n - 1, lb + 2)])
    regions.sort()
    if regions:
        merged: list[list[int]] = []
        for r in regions:
            if merged and r[0] <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], r[1])
            else:
                merged.append(r)
        regions = merged

    def _run_segmented() -> tuple[bool, int | None, str | None, dict | None, list]:
        """Run per-segment beam with current regions.
        Returns (success, fail_codon, optimized_dna, opt_report, freq_warnings)."""
        segs = []
        cursor = 0
        for r in sorted(regions):
            if cursor < r[0]:
                segs.append(("oligo", cursor, r[0] - 1))
            segs.append(("gblock", r[0], r[1]))
            cursor = r[1] + 1
        if cursor < n:
            segs.append(("oligo", cursor, n - 1))

        dna_parts = []
        all_warnings = []
        total_loss = 0.0
        ideal_total = 0.0
        all_ok = True
        seg_reports = []

        for kind, start, end in segs:
            sub = protein[start:end + 1]
            sub_dna = reverse_translate(sub)
            if kind == "gblock":
                kwargs = dict(
                    gc_min=GBLOCK_GC_MIN, gc_max=GBLOCK_GC_MAX,
                    gc_window=gc_window,
                    homo_gc=GBLOCK_HOMO_GC, homo_at=GBLOCK_HOMO_AT,
                    avoid_patterns=avoid_patterns,
                    min_codon_frequency=min_codon_frequency,
                    beam_k=beam_k,
                )
            else:
                kwargs = dict(
                    gc_min=gc_min, gc_max=gc_max, gc_window=gc_window,
                    homo_gc=homo_gc, homo_at=homo_at,
                    avoid_patterns=avoid_patterns,
                    min_codon_frequency=min_codon_frequency,
                    beam_k=beam_k,
                )
            refined, warnings, report = refine_gc_windows(sub_dna, raw_table, **kwargs)
            if report.get("source") != "beam":
                local_fail = report.get("failed_at_codon", 0)
                return False, start + local_fail, None, None, []
            dna_parts.append(refined)
            for w in warnings:
                w["start"] += start * 3
                w["end"] += start * 3
                w["codon_index"] += start
            all_warnings.extend(warnings)
            total_loss += report.get("total_loss", 0.0)
            ideal_total += report.get("ideal_score", 0.0)
            seg_reports.append({"kind": kind, "start_aa": start, "end_aa": end,
                                "source": report.get("source"), "loss": round(report.get("total_loss", 0.0), 1)})

        opt_dna = "".join(dna_parts)
        opt_report = {
            "source": "beam",
            "beam_k": beam_k,
            "ideal_score": round(ideal_total, 1),
            "chosen_score": round(ideal_total - total_loss, 1),
            "total_loss": round(total_loss, 1),
            "segments": seg_reports,
        }
        return True, None, opt_dna, opt_report, all_warnings

    success_dna = None
    success_report = None
    success_warnings: list = []

    for _it in range(max_iterations):
        if progress_callback:
            progress_callback("gblock_detect", _it, max_iterations)
            time.sleep(0)  # yield GIL for SSE
        ok, fail_at, opt_dna_result, opt_report_result, warn_result = _run_segmented()
        if ok:
            success_dna = opt_dna_result
            success_report = opt_report_result
            success_warnings = warn_result
            break
        if fail_at is None:
            break

        # Find or grow a region around the failure point
        target = None
        for r in regions:
            if r[0] - grow_step <= fail_at <= r[1] + grow_step:
                target = r
                break
        if target is not None:
            target[0] = max(0, target[0] - grow_step)
            target[1] = min(n - 1, target[1] + grow_step)
        else:
            regions.append([
                max(0, fail_at - initial_radius),
                min(n - 1, fail_at + initial_radius),
            ])

        # Merge overlapping regions
        regions.sort()
        merged: list[list[int]] = []
        for r in regions:
            if merged and r[0] <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], r[1])
            else:
                merged.append(r)
        regions = merged

        # Safety: bail if we've covered > 50% of the protein
        covered = sum(b - a + 1 for a, b in regions)
        if covered > n * 0.5:
            break

    if success_dna is None:
        # Couldn't make the beam pass even by growing — give up on auto-detect.
        return [], None, None, []

    # Phase 3: shrink.  The growth loop uses grow_step-sized chunks, so the
    # final region is typically 5–10 aa larger than strictly necessary on
    # each side.  We try larger shrink amounts first and fall back to half
    # that if the bigger shrink breaks the beam.  Boundaries are never
    # allowed to land inside a repetitive region (would leave a cross-hyb-
    # prone fragment in the tiled region).
    def _inside_repeat(pos: int) -> bool:
        return any(la < pos <= lb for (la, lb) in repeat_regions)

    def _left_crosses_repeat(lo: int, hi: int) -> bool:
        """Shrinking the left edge from lo to hi would exclude a repeat
        whose left edge sat between them."""
        return any(lo < la <= hi for (la, lb) in repeat_regions)

    def _right_crosses_repeat(lo: int, hi: int) -> bool:
        """Shrinking the right edge from hi down to lo would exclude a
        repeat whose right edge sat between them."""
        return any(lo < lb <= hi for (la, lb) in repeat_regions)

    shrink_attempts = (grow_step, max(1, grow_step // 2))
    for reg_idx, r in enumerate(regions):
        if progress_callback:
            progress_callback("gblock_shrink", reg_idx, len(regions))
            time.sleep(0)

        # Shrink from the LEFT
        for amt in shrink_attempts:
            new_a = r[0] + amt
            if new_a >= r[1]:
                continue
            if _inside_repeat(new_a):
                continue
            if _left_crosses_repeat(r[0], new_a):
                continue
            original_a = r[0]
            r[0] = new_a
            ok, _fa, opt_dna_r, opt_rep_r, warn_r = _run_segmented()
            if ok:
                success_dna = opt_dna_r
                success_report = opt_rep_r
                success_warnings = warn_r
                break
            r[0] = original_a

        # Shrink from the RIGHT
        for amt in shrink_attempts:
            new_b = r[1] - amt
            if new_b <= r[0]:
                continue
            if _inside_repeat(new_b + 1):
                continue
            if _right_crosses_repeat(new_b, r[1]):
                continue
            original_b = r[1]
            r[1] = new_b
            ok, _fa, opt_dna_r, opt_rep_r, warn_r = _run_segmented()
            if ok:
                success_dna = opt_dna_r
                success_report = opt_rep_r
                success_warnings = warn_r
                break
            r[1] = original_b

    # Phase 4: pad each region to IDT minimum (42 aa = 126 bp ≥ 125 bp).
    # Prefer padding toward a nearby repetitive region (pull it fully
    # inside the gBlock so the tiled flank stays non-repetitive) before
    # defaulting to symmetric padding.
    padded = []
    for a, b in regions:
        length = b - a + 1
        if length < min_aa:
            deficit = min_aa - length
            for (la, lb) in repeat_regions:
                if deficit <= 0:
                    break
                if lb + 1 == a and la >= a - deficit:
                    absorb = min(deficit, a - la)
                    a -= absorb
                    deficit -= absorb
                elif la == b + 1 and lb <= b + deficit:
                    absorb = min(deficit, lb - b)
                    b += absorb
                    deficit -= absorb
            if deficit > 0:
                expand_left = deficit // 2
                expand_right = deficit - expand_left
                a = max(0, a - expand_left)
                b = min(n - 1, b + expand_right)
                if (b - a + 1) < min_aa:
                    if a == 0:
                        b = min(n - 1, a + min_aa - 1)
                    else:
                        a = max(0, b - min_aa + 1)
        padded.append((a, b))

    # Post-pad merge: symmetric padding can push two adjacent regions into
    # overlap (common on short proteins with closely-spaced linkers).  Fuse
    # them into a single larger gBlock rather than emit invalid overlapping
    # segments.
    padded.sort()
    fused: list[tuple[int, int]] = []
    for a, b in padded:
        if fused and a <= fused[-1][1] + 1:
            fa, fb = fused[-1]
            fused[-1] = (fa, max(fb, b))
        else:
            fused.append((a, b))

    result = [
        GBlockRegionInput(start_aa=a, end_aa=b, label=f"auto ({a}-{b})")
        for (a, b) in fused
    ]
    return result, success_dna, success_report, success_warnings


def _build_suggested_gblock_info(
    suggestions: list[GBlockRegionInput],
    opt_dna: str,
) -> list[GBlockSegmentInfo]:
    """Convert auto-detect suggestions into GBlockSegmentInfo for the response."""
    def _max_run(seq, bases):
        best = cur = 0
        for ch in seq:
            if ch in bases:
                cur += 1
                best = max(best, cur)
            else:
                cur = 0
        return best

    result = []
    for gb in suggestions:
        bp_start = gb.start_aa * 3
        bp_end = (gb.end_aa + 1) * 3
        region_dna = opt_dna[bp_start:bp_end] if bp_end <= len(opt_dna) else ""
        result.append(GBlockSegmentInfo(
            start_aa=gb.start_aa, end_aa=gb.end_aa,
            start_bp=bp_start, end_bp=bp_end,
            label=gb.label,
            len_bp=len(region_dna),
            valid_size=GBLOCK_MIN_BP <= len(region_dna) <= GBLOCK_MAX_BP,
            gc=round(gc_content(region_dna), 3) if region_dna else 0,
            max_g_run=_max_run(region_dna, "G"),
            max_c_run=_max_run(region_dna, "C"),
        ))
    return result


def _run_codon_optimize(req: CodonOptRequest, progress_callback=None) -> CodonOptResponse:
    """Core optimization pipeline. progress_callback(stage: str, current: int, total: int)
    is called periodically so streaming endpoints can report progress."""
    from dnachisel import (
        DnaOptimizationProblem, AvoidPattern,
        EnforceGCContent, EnforceTranslation, UniquifyAllKmers,
        reverse_translate,
    )

    protein = req.protein_sequence.upper().replace(" ", "").replace("\n", "")

    bad = set(protein) - set("ACDEFGHIKLMNPQRSTVWY*")
    if bad:
        raise HTTPException(400, f"Non-amino-acid characters: {bad}")

    if protein.endswith("*"):
        protein = protein[:-1]

    # Resolve codon table (used for frequency lookup by the beam search)
    if req.custom_codon_table:
        usage_table = parse_codon_table(req.custom_codon_table)
        if not usage_table:
            raise HTTPException(400, "Could not parse custom codon table")
        table_name = "Custom"
    elif req.species in CUSTOM_TABLES:
        usage_table = parse_codon_table(CUSTOM_TABLES[req.species])
        if not usage_table:
            raise HTTPException(400, f"Could not parse bundled table: {req.species}")
        table_name = req.species
    else:
        import python_codon_tables as pct
        usage_table = pct.get_codons_table(req.species)
        table_name = req.species

    # Build raw /1000 table (needed for both DP refinement and per-codon detail)
    VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
    if req.custom_codon_table:
        raw_table = parse_codon_table_raw(req.custom_codon_table)
    elif req.species in CUSTOM_TABLES:
        raw_table = parse_codon_table_raw(CUSTOM_TABLES[req.species])
    else:
        # Built-in tables have frequencies summing to ~1 per aa; scale to /1000
        raw_table = {}
        for aa in usage_table:
            if aa not in VALID_AA:
                continue
            raw_table[aa] = {c: round(v * 1000, 1) for c, v in usage_table[aa].items()}

    # Build constraints
    constraints = _build_constraints(
        req.avoid_homopolymers_gc, req.avoid_homopolymers_at,
        req.gc_min, req.gc_max, req.gc_window, req.avoid_patterns,
    )

    # Upfront feasibility check: if any amino acid in the protein has no codon
    # clearing min_codon_frequency, the user's floor is impossible for this
    # protein and table. Fail with a specific per-AA error before wasting time
    # on the beam.
    if req.min_codon_frequency > 0:
        impossible: dict[str, float] = {}
        for aa in set(protein):
            if aa not in VALID_AA:
                continue
            codons_for_aa = raw_table.get(aa, {})
            if not codons_for_aa:
                continue
            best = max(codons_for_aa.values())
            if best < req.min_codon_frequency:
                impossible[aa] = best
        if impossible:
            parts = []
            for aa, best in sorted(impossible.items(), key=lambda kv: kv[1]):
                positions = [i for i, a in enumerate(protein) if a == aa]
                pos_preview = ", ".join(str(p) for p in positions[:3])
                if len(positions) > 3:
                    pos_preview += f", … ({len(positions)} total)"
                parts.append(
                    f"{aa} (best {best:.1f}/1000, positions {pos_preview})"
                )
            raise HTTPException(
                400,
                "Min codon frequency floor is impossible for this protein in "
                f"the {table_name} table. No synonymous codon above "
                f"{req.min_codon_frequency:.0f}/1000 exists for: "
                + "; ".join(parts)
                + f". Lower the floor to ≤ {min(impossible.values()):.1f}/1000 or use a different codon table."
            )

    # Beam search is the sole optimizer. Start from a naive reverse translation
    # — the beam only reads the amino acid sequence from this, so the starting
    # codon choices are discarded.
    opt_dna = reverse_translate(protein)

    from gc_refinement import refine_gc_windows

    # Progress throttle helper
    def _make_beam_progress(n_codons: int, offset: int = 0, total: int | None = None):
        if n_codons <= 0:
            return None
        total = total or n_codons
        step = max(1, total // 50)
        def cb(ci: int, _n: int) -> None:
            if progress_callback is None:
                return
            pos = offset + ci
            if pos % step == 0 or ci == _n - 1:
                progress_callback("beam", pos, total)
                time.sleep(0)
        return cb

    gblocks = sorted(req.gblock_regions, key=lambda g: g.start_aa) if req.gblock_regions else []
    auto_suggested: list[GBlockRegionInput] = []
    gblocks_auto_applied = False

    if not gblocks:
        # No gBlock regions — original single-beam path
        _beam_progress = _make_beam_progress(len(opt_dna) // 3)
        opt_dna, freq_warnings, opt_report = refine_gc_windows(
            opt_dna, raw_table,
            gc_min=req.gc_min, gc_max=req.gc_max, gc_window=req.gc_window,
            homo_gc=req.avoid_homopolymers_gc, homo_at=req.avoid_homopolymers_at,
            avoid_patterns=req.avoid_patterns,
            min_codon_frequency=req.min_codon_frequency,
            beam_k=req.beam_k,
            progress_callback=_beam_progress,
        )

        # AUTO-DETECT: if beam failed, iteratively find minimal gBlock regions
        # and re-run with per-segment optimization.
        if opt_report.get("source") in ("beam_pruned_empty", "beam_no_solution"):
            if progress_callback:
                progress_callback("gblock_detect", 0, 1)

            auto_suggested, cached_dna, cached_report, cached_warnings = _auto_detect_gblocks(
                protein, raw_table,
                gc_min=req.gc_min, gc_max=req.gc_max, gc_window=req.gc_window,
                homo_gc=req.avoid_homopolymers_gc, homo_at=req.avoid_homopolymers_at,
                avoid_patterns=req.avoid_patterns,
                min_codon_frequency=req.min_codon_frequency,
                beam_k=req.beam_k,
                progress_callback=progress_callback,
            )

            if auto_suggested:
                gblocks = sorted(auto_suggested, key=lambda g: g.start_aa)
                gblocks_auto_applied = True
                # Use the cached result from the successful auto-detect iteration
                # instead of re-running the per-segment beam.
                if cached_dna:
                    opt_dna = cached_dna
                    opt_report = cached_report
                    freq_warnings = cached_warnings
                if progress_callback:
                    progress_callback("gblock_detect", 1, 1)

    if gblocks and not gblocks_auto_applied:
        # gBlock-aware per-segment optimization: run the beam independently on
        # each segment (gBlock vs oligo-tiled), using appropriate constraints.
        n_codons = len(protein)
        segments = []
        cursor = 0
        for gb in gblocks:
            if cursor < gb.start_aa:
                segments.append(("oligo", cursor, gb.start_aa - 1))
            segments.append(("gblock", gb.start_aa, gb.end_aa))
            cursor = gb.end_aa + 1
        if cursor < n_codons:
            segments.append(("oligo", cursor, n_codons - 1))

        dna_parts = []
        all_freq_warnings = []
        total_loss = 0.0
        ideal_total = 0.0
        all_succeeded = True
        seg_reports = []
        codon_offset = 0

        for kind, start, end in segments:
            sub_protein = protein[start:end + 1]
            sub_dna = reverse_translate(sub_protein)

            if kind == "gblock":
                kwargs = dict(
                    gc_min=GBLOCK_GC_MIN, gc_max=GBLOCK_GC_MAX,
                    gc_window=req.gc_window,
                    homo_gc=GBLOCK_HOMO_GC, homo_at=GBLOCK_HOMO_AT,
                    avoid_patterns=req.avoid_patterns,
                    min_codon_frequency=req.min_codon_frequency,
                    beam_k=req.beam_k,
                    progress_callback=_make_beam_progress(len(sub_protein), codon_offset, n_codons),
                )
            else:
                kwargs = dict(
                    gc_min=req.gc_min, gc_max=req.gc_max,
                    gc_window=req.gc_window,
                    homo_gc=req.avoid_homopolymers_gc, homo_at=req.avoid_homopolymers_at,
                    avoid_patterns=req.avoid_patterns,
                    min_codon_frequency=req.min_codon_frequency,
                    beam_k=req.beam_k,
                    progress_callback=_make_beam_progress(len(sub_protein), codon_offset, n_codons),
                )

            refined, warnings, report = refine_gc_windows(sub_dna, raw_table, **kwargs)
            dna_parts.append(refined)
            for w in warnings:
                w["start"] += start * 3
                w["end"] += start * 3
                w["codon_index"] += start
            all_freq_warnings.extend(warnings)
            total_loss += report.get("total_loss", 0.0)
            ideal_total += report.get("ideal_score", 0.0)
            source = report.get("source", "?")
            if source != "beam":
                all_succeeded = False
            seg_reports.append({
                "kind": kind,
                "start_aa": start, "end_aa": end,
                "source": source,
                "loss": round(report.get("total_loss", 0.0), 1),
            })
            codon_offset += len(sub_protein)

        opt_dna = "".join(dna_parts)
        freq_warnings = all_freq_warnings
        opt_report = {
            "source": "beam" if all_succeeded else "beam_partial",
            "beam_k": req.beam_k,
            "ideal_score": round(ideal_total, 1),
            "chosen_score": round(ideal_total - total_loss, 1),
            "total_loss": round(total_loss, 1),
            "segments": seg_reports,
        }

    back_protein = str(Seq(opt_dna).translate())

    # Per-codon detail using raw /1000 values
    codons_detail = []
    for i in range(0, len(opt_dna) - 2, 3):
        codon = opt_dna[i:i+3]
        aa = str(Seq(codon).translate())
        aa_raw = raw_table.get(aa, {})
        per_k = aa_raw.get(codon, 0.0)
        codons_detail.append(CodonDetail(
            amino_acid=aa, codon=codon, per_thousand=round(per_k, 1)
        ))

    # Actual minimum codon frequency in the output sequence (diagnostic).
    # The user can compare this against their requested floor to decide
    # whether to manually raise it.
    if codons_detail:
        opt_report["min_frequency_achieved"] = round(min(c.per_thousand for c in codons_detail), 1)

    table_entries = []
    for aa in sorted(raw_table.keys()):
        if aa not in VALID_AA:
            continue
        for codon, per_k in sorted(raw_table[aa].items(), key=lambda x: -x[1]):
            table_entries.append(CodonTableEntry(
                amino_acid=aa, codon=codon, per_thousand=round(per_k, 1)
            ))

    # Re-evaluate constraints against the (possibly refined) sequence
    eval_problem = DnaOptimizationProblem(
        sequence=opt_dna,
        constraints=constraints,
    )
    kmer_problem = DnaOptimizationProblem(
        sequence=opt_dna,
        constraints=[UniquifyAllKmers(req.uniquify_kmers)],
    )

    # Constraint results — tag each as enforced or detection-only
    summary = []
    for ev in eval_problem.constraints_evaluations():
        summary.append({
            "constraint": str(ev.specification),
            "passing": bool(ev.passes),
            "message": str(ev.message) if hasattr(ev, "message") else "",
            "enforced": True,
        })
    for ev in kmer_problem.constraints_evaluations():
        summary.append({
            "constraint": str(ev.specification),
            "passing": bool(ev.passes),
            "message": str(ev.message) if hasattr(ev, "message") else "",
            "enforced": False,
        })

    # Min codon frequency constraint
    if req.min_codon_frequency > 0:
        floor_pass = len(freq_warnings) == 0
        summary.append({
            "constraint": f"MinCodonFrequency(min:{req.min_codon_frequency}/1000)",
            "passing": floor_pass,
            "message": "" if floor_pass else f"{len(freq_warnings)} codon(s) below {req.min_codon_frequency}/1000",
            "enforced": True,
        })

    # Warnings: scan + k-mer failing locations (grouped by matching sequence)
    # Generate two sets of warnings: oligo-strict for non-gBlock regions,
    # IDT-grade for gBlock regions.
    opt_warnings = _scan_warnings(
        opt_dna, gc_min=req.gc_min, gc_max=req.gc_max, gc_window=req.gc_window,
        homo_gc=req.avoid_homopolymers_gc, homo_at=req.avoid_homopolymers_at,
    )
    opt_warnings.extend(_kmer_warnings(opt_dna, kmer_problem, req.uniquify_kmers))
    opt_warnings.extend(freq_warnings)

    # Filter warnings inside gBlock regions: suppress homopolymer/gc_window
    # warnings that are within IDT's looser limits.
    if gblocks:
        gb_ranges = [(gb.start_aa * 3, (gb.end_aa + 1) * 3) for gb in gblocks]

        def _inside_gblock(w_start: int, w_end: int) -> bool:
            for gs, ge in gb_ranges:
                if w_start >= gs and w_end <= ge:
                    return True
            return False

        def _overlaps_gblock(w_start: int, w_end: int) -> bool:
            for gs, ge in gb_ranges:
                if w_start < ge and w_end > gs:
                    return True
            return False

        filtered = []
        for w in opt_warnings:
            ws, we = w.get("start", 0), w.get("end", 0)
            kind = w.get("kind", "")

            if kind == "homopolymer" and _inside_gblock(ws, we):
                # Re-check with IDT limits: suppress if run < IDT threshold
                region = opt_dna[ws:we]
                max_gc_run = max((len(m.group()) for m in re.finditer(r"[GC]+", region)), default=0)
                max_at_run = max((len(m.group()) for m in re.finditer(r"[AT]+", region)), default=0)
                if max_gc_run < GBLOCK_HOMO_GC and max_at_run < GBLOCK_HOMO_AT:
                    continue  # within IDT limits, suppress

            if kind == "gc_window" and _overlaps_gblock(ws, we):
                # GC window that spans into a gBlock boundary — suppress, since
                # the GC content in the gBlock portion is IDT's concern, not ours.
                continue

            filtered.append(w)
        opt_warnings = filtered

    all_pass = eval_problem.all_constraints_pass()
    # Override all_pass if gblocks are present — the DNA Chisel constraint
    # evaluator doesn't know about gBlock-aware relaxation, so its failures
    # inside gBlock regions are not real failures.
    if gblocks and not all_pass:
        # Re-check: are all failures inside gBlock regions?
        gb_ranges_set = [(gb.start_aa * 3, (gb.end_aa + 1) * 3) for gb in gblocks]
        real_failures = []
        for ev in eval_problem.constraints_evaluations():
            if not ev.passes:
                # Check if all failing locations are inside gBlock regions
                locs = getattr(ev, 'locations', []) if hasattr(ev, 'locations') else []
                if not locs:
                    real_failures.append(ev)
                else:
                    for loc in locs:
                        loc_start = loc.start if hasattr(loc, 'start') else 0
                        loc_end = loc.end if hasattr(loc, 'end') else len(opt_dna)
                        inside = any(gs <= loc_start and loc_end <= ge for gs, ge in gb_ranges_set)
                        if not inside:
                            real_failures.append(ev)
                            break
        all_pass = len(real_failures) == 0

    # Build gBlock segment info for the response
    def _max_run(seq, bases):
        best = cur = 0
        for ch in seq:
            if ch in bases:
                cur += 1
                best = max(best, cur)
            else:
                cur = 0
        return best

    gb_segments = []
    for gb in gblocks:
        bp_start = gb.start_aa * 3
        bp_end = (gb.end_aa + 1) * 3
        region_dna = opt_dna[bp_start:bp_end]
        gb_segments.append(GBlockSegmentInfo(
            start_aa=gb.start_aa, end_aa=gb.end_aa,
            start_bp=bp_start, end_bp=bp_end,
            label=gb.label,
            len_bp=len(region_dna),
            valid_size=GBLOCK_MIN_BP <= len(region_dna) <= GBLOCK_MAX_BP,
            gc=round(gc_content(region_dna), 3) if region_dna else 0,
            max_g_run=_max_run(region_dna, "G"),
            max_c_run=_max_run(region_dna, "C"),
        ))

    gc_lo, gc_hi = _gc_window_range(opt_dna, req.gc_window)
    return CodonOptResponse(
        input_protein=protein,
        optimized_dna=opt_dna,
        back_translated_protein=back_protein,
        length_dna=len(opt_dna),
        length_protein=len(protein),
        gc_content=round(gc_content(opt_dna), 3),
        gc_window_min=round(gc_lo, 3) if gc_lo is not None else None,
        gc_window_max=round(gc_hi, 3) if gc_hi is not None else None,
        gc_window_size=req.gc_window,
        constraints_pass=all_pass,
        constraints_summary=summary,
        warnings=[SequenceWarning(**w) for w in opt_warnings],
        codons=codons_detail,
        codon_table=table_entries,
        codon_table_name=table_name,
        optimization_report=opt_report,
        gblock_segments=gb_segments,
        suggested_gblocks=_build_suggested_gblock_info(auto_suggested, opt_dna),
        gblocks_auto_applied=gblocks_auto_applied,
    )


@app.post("/api/codon-optimize", response_model=CodonOptResponse)
def codon_optimize(req: CodonOptRequest):
    return _run_codon_optimize(req)


@app.post("/api/codon-optimize-stream")
async def codon_optimize_stream(req: CodonOptRequest):
    """Streams Server-Sent Events with progress updates and the final result.

    Event format (newline-delimited JSON, one event per line):
      {"type":"progress","stage":"beam","current":42,"total":300}
      {"type":"result","data":{...CodonOptResponse...}}
      {"type":"error","message":"..."}
    """
    q: queue.Queue = queue.Queue()
    SENTINEL = object()

    def emit_progress(stage: str, current: int, total: int) -> None:
        q.put({"type": "progress", "stage": stage, "current": current, "total": total})

    def worker() -> None:
        try:
            response = _run_codon_optimize(req, progress_callback=emit_progress)
            q.put({"type": "result", "data": response.model_dump()})
        except HTTPException as e:
            q.put({"type": "error", "message": e.detail, "status": e.status_code})
        except Exception as e:
            q.put({"type": "error", "message": str(e), "status": 500})
        finally:
            q.put(SENTINEL)

    threading.Thread(target=worker, daemon=True).start()

    async def event_stream():
        loop = asyncio.get_event_loop()
        last_progress: tuple[str, int] | None = None
        while True:
            item = await loop.run_in_executor(None, q.get)
            if item is SENTINEL:
                break
            # Throttle progress events: only emit when stage or current value changes
            if item.get("type") == "progress":
                key = (item["stage"], item["current"])
                if key == last_progress:
                    continue
                last_progress = key
            yield f"data: {json.dumps(item)}\n\n"

    return StreamingResponse(event_stream(), media_type="text/event-stream")


# ---------------------------------------------------------------------------
# Constraint checking (for manual edits)
# ---------------------------------------------------------------------------

class CheckConstraintsRequest(BaseModel):
    dna_sequence: str = Field(..., description="DNA sequence to check")
    avoid_homopolymers_gc: int = Field(4, ge=3, le=8)
    avoid_homopolymers_at: int = Field(6, ge=3, le=8)
    gc_min: float = Field(0.25, ge=0.0, le=1.0)
    gc_max: float = Field(0.75, ge=0.0, le=1.0)
    gc_window: int = Field(50, ge=20, le=200)
    uniquify_kmers: int = Field(10, ge=6, le=20)
    avoid_patterns: list[str] = Field(default_factory=list)
    min_codon_frequency: float = Field(10.0, ge=0.0, le=100.0)


class CheckConstraintsResponse(BaseModel):
    constraints_pass: bool
    constraints_summary: list[dict]
    gc_content: float
    warnings: list[SequenceWarning]


def _build_constraints(homo_gc: int, homo_at: int, gc_min: float, gc_max: float,
                       gc_window: int, avoid_patterns: list[str]):
    """Build the standard DNA Chisel constraint list used by both endpoints."""
    from dnachisel import AvoidPattern, EnforceGCContent, EnforceTranslation
    constraints = [EnforceTranslation()]
    for base in "GC":
        constraints.append(AvoidPattern(base * homo_gc))
    for base in "AT":
        constraints.append(AvoidPattern(base * homo_at))
    constraints.append(EnforceGCContent(mini=gc_min, maxi=gc_max, window=gc_window))
    for pat in avoid_patterns:
        constraints.append(AvoidPattern(pat))
    return constraints


def _gc_window_range(dna: str, window: int) -> tuple[float | None, float | None]:
    """Return (min, max) GC fraction across all sliding windows of given size.
    Returns (None, None) if the sequence is shorter than the window."""
    s = dna.upper()
    if len(s) < window:
        return (None, None)
    lo = 1.0
    hi = 0.0
    for i in range(0, len(s) - window + 1):
        w = s[i:i + window]
        gc = (w.count('G') + w.count('C')) / window
        if gc < lo:
            lo = gc
        if gc > hi:
            hi = gc
    return (lo, hi)


def _scan_warnings(dna: str, gc_min: float = 0.25, gc_max: float = 0.75,
                    gc_window: int = 50,
                    homo_gc: int = 4, homo_at: int = 6) -> list[dict]:
    """Scan for soft warnings — things worth knowing about but not constraint violations."""
    warnings = []
    s = dna.upper()

    # Homopolymer runs
    for base, threshold in [('G', homo_gc), ('C', homo_gc), ('A', homo_at), ('T', homo_at)]:
        for m in re.finditer(f'{base}{{{threshold},}}', s):
            warnings.append({
                "kind": "homopolymer",
                "message": f"{base}x{m.end() - m.start()} at {m.start()}-{m.end()}",
                "start": m.start(),
                "end": m.end(),
            })

    # GC content in sliding windows — use same thresholds as constraints
    window = gc_window
    for i in range(0, len(s) - window + 1, 10):
        w = s[i:i + window]
        gc = (w.count('G') + w.count('C')) / len(w)
        if gc < gc_min or gc > gc_max:
            warnings.append({
                "kind": "gc_window",
                "message": f"GC {gc:.0%} in {window}bp window at {i}-{i + window}",
                "start": i,
                "end": i + window,
            })

    return warnings


def _kmer_warnings(dna: str, kmer_problem, kmer_size: int) -> list[dict]:
    """Extract k-mer repeat warnings from DNA Chisel, grouped by shared sequence.

    Uses DNA Chisel's merged failing regions, then for each pair of regions
    that share an exact k-mer, assigns them the same group letter.
    """
    s = dna.upper()

    # Collect all failing locations from DNA Chisel
    locations = []
    for ev in kmer_problem.constraints_evaluations():
        if not ev.passes:
            for loc in ev.locations:
                locations.append((loc.start, loc.end))

    if not locations:
        return []

    # Merge overlapping locations into contiguous regions
    locations.sort()
    merged = []
    for start, end in locations:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    # For each region, find all k-mers it contains
    region_kmers: list[set[str]] = []
    for start, end in merged:
        kmers = set()
        for i in range(start, min(end - kmer_size + 1, len(s) - kmer_size + 1)):
            kmers.add(s[i:i + kmer_size])
        region_kmers.append(kmers)

    # Find which k-mers are actually repeated (appear in 2+ regions)
    kmer_to_regions: dict[str, list[int]] = {}
    for idx, kmers in enumerate(region_kmers):
        for kmer in kmers:
            kmer_to_regions.setdefault(kmer, []).append(idx)
    repeated_kmers = {k: v for k, v in kmer_to_regions.items() if len(v) > 1}

    # Group regions by shared repeated k-mers using union-find
    parent = list(range(len(merged)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for kmer, region_ids in repeated_kmers.items():
        for i in range(1, len(region_ids)):
            union(region_ids[0], region_ids[i])

    # For each group, find the shared k-mer(s) that caused the grouping
    root_to_group: dict[int, int] = {}
    root_to_shared_kmer: dict[int, str] = {}
    next_group = 1

    # First pass: find shared k-mers per group
    group_members: dict[int, list[int]] = {}
    for idx in range(len(merged)):
        root = find(idx)
        group_members.setdefault(root, []).append(idx)

    for root, members in group_members.items():
        if len(members) < 2:
            continue
        # Find k-mers shared across ALL members of this group
        shared = region_kmers[members[0]].copy()
        for m in members[1:]:
            shared &= region_kmers[m]
        if not shared:
            # If no single k-mer is shared by ALL, find the most common one
            all_kmers: dict[str, int] = {}
            for m in members:
                for k in region_kmers[m]:
                    if k in repeated_kmers:
                        all_kmers[k] = all_kmers.get(k, 0) + 1
            shared = {max(all_kmers, key=all_kmers.get)} if all_kmers else set()
        if shared:
            root_to_shared_kmer[root] = sorted(shared)[0]

    result = []
    for idx, (start, end) in enumerate(merged):
        root = find(idx)
        if root not in root_to_group:
            root_to_group[root] = next_group
            next_group += 1
        group = root_to_group[root]

        # Use the shared k-mer for the group, not the region-start k-mer
        rep_kmer = root_to_shared_kmer.get(root, s[start:start + kmer_size])

        result.append({
            "kind": "repeat_kmer",
            "message": f"Repeated {kmer_size}-mer at {start}-{end} ({rep_kmer})",
            "start": start,
            "end": end,
            "group": group,
        })

    # Remove groups with only one member — a solo "repeat" is meaningless
    from collections import Counter
    group_counts = Counter(w["group"] for w in result)
    result = [w for w in result if group_counts[w["group"]] > 1]

    # Re-number groups sequentially
    old_to_new: dict[int, int] = {}
    new_num = 1
    for w in result:
        if w["group"] not in old_to_new:
            old_to_new[w["group"]] = new_num
            new_num += 1
        w["group"] = old_to_new[w["group"]]

    return result


@app.post("/api/check-constraints", response_model=CheckConstraintsResponse)
def check_constraints(req: CheckConstraintsRequest):
    from dnachisel import (
        DnaOptimizationProblem, AvoidPattern,
        EnforceGCContent, EnforceTranslation, UniquifyAllKmers,
    )

    dna = req.dna_sequence.upper().replace(" ", "").replace("\n", "")

    bad = set(dna) - set("ATCGN")
    if bad:
        raise HTTPException(400, f"Non-DNA characters: {bad}")

    constraints = _build_constraints(
        req.avoid_homopolymers_gc, req.avoid_homopolymers_at,
        req.gc_min, req.gc_max, req.gc_window, req.avoid_patterns,
    )

    problem = DnaOptimizationProblem(
        sequence=dna,
        constraints=constraints,
    )

    summary = []
    constraint_warnings = []
    for ev in problem.constraints_evaluations():
        summary.append({
            "constraint": str(ev.specification),
            "passing": bool(ev.passes),
            "message": str(ev.message) if hasattr(ev, "message") else "",
            "enforced": True,
        })
        # Extract failing locations from constraints to show on the track
        if not ev.passes and "EnforceGCContent" in str(ev.specification):
            for loc in ev.locations:
                # Find the worst 50bp window within the failing region
                worst_gc = None
                worst_start = loc.start
                scan_start = max(0, loc.start - req.gc_window)
                scan_end = min(len(dna), loc.end + req.gc_window)
                for pos in range(scan_start, scan_end - req.gc_window + 1):
                    w = dna[pos:pos + req.gc_window]
                    gc = (w.count('G') + w.count('C')) / len(w)
                    if gc < req.gc_min or gc > req.gc_max:
                        if worst_gc is None or abs(gc - 0.5) > abs(worst_gc - 0.5):
                            worst_gc = gc
                            worst_start = pos
                if worst_gc is not None:
                    direction = "high" if worst_gc > req.gc_max else "low"
                    constraint_warnings.append({
                        "kind": "gc_window",
                        "message": f"GC {worst_gc:.0%} ({direction}) at {loc.start}-{loc.end}",
                        "start": loc.start,
                        "end": loc.end,
                    })
    # Run k-mer check separately for grouped warnings
    kmer_problem = DnaOptimizationProblem(
        sequence=dna,
        constraints=[UniquifyAllKmers(req.uniquify_kmers)],
    )
    for ev in kmer_problem.constraints_evaluations():
        summary.append({
            "constraint": str(ev.specification),
            "passing": bool(ev.passes),
            "message": str(ev.message) if hasattr(ev, "message") else "",
            "enforced": False,
        })

    warnings = _scan_warnings(dna, gc_min=req.gc_min, gc_max=req.gc_max, gc_window=req.gc_window,
                              homo_gc=req.avoid_homopolymers_gc, homo_at=req.avoid_homopolymers_at)
    # Add GC constraint warnings
    existing = {(w["kind"], w["start"], w["end"]) for w in warnings}
    for cw in constraint_warnings:
        if (cw["kind"], cw["start"], cw["end"]) not in existing:
            warnings.append(cw)
    # Add grouped k-mer warnings
    warnings.extend(_kmer_warnings(dna, kmer_problem, req.uniquify_kmers))

    return CheckConstraintsResponse(
        constraints_pass=problem.all_constraints_pass(),
        constraints_summary=summary,
        gc_content=round(gc_content(dna), 3),
        warnings=[SequenceWarning(**w) for w in warnings],
    )


# ---------------------------------------------------------------------------
# Utility endpoints
# ---------------------------------------------------------------------------

@app.get("/api/codon-tables")
def codon_tables():
    import python_codon_tables as pct
    return {
        "tables": sorted(pct.available_codon_tables_names),
        "custom_tables": sorted(CUSTOM_TABLES.keys()),
    }


@app.get("/api/health")
def health():
    return {"status": "ok"}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
