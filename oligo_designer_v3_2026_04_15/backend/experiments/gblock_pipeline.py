"""
gBlock-aware codon optimization + oligo design pipeline prototype.

End-to-end demo: takes a protein, gBlock region markings, and constraints,
and produces (a) codon-optimized DNA with per-region constraint handling,
(b) an oligo plan where gBlock regions are emitted as whole fragments and
everything else is tiled into overlapping oligos.

Run:
    cd backend && python experiments/gblock_pipeline.py

IDT gBlock specs (from idtdna.com, April 2026):
  - Size:  125–3000 bp  (42–1000 aa)
  - Homo G/C: < 6  (max run of 5)
  - Homo A/T: < 10 (max run of 9)
  - GC content: avoid < 25% and > 75%
  - Pricing: ~$0.07/bp
"""

from __future__ import annotations
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path

_THIS = Path(__file__).resolve()
_BACKEND = _THIS.parent.parent
sys.path.insert(0, str(_BACKEND))

from Bio.Seq import Seq
from api import parse_codon_table_raw
from gc_refinement import refine_gc_windows

PROJECT_ROOT = _BACKEND.parent.parent.parent


# ---------------------------------------------------------------------------
# IDT gBlock constraint constants
# ---------------------------------------------------------------------------

GBLOCK_MIN_BP = 125
GBLOCK_MAX_BP = 3000
GBLOCK_HOMO_GC = 6   # max G/C run < 6
GBLOCK_HOMO_AT = 10   # max A/T run < 10
GBLOCK_GC_MIN = 0.25
GBLOCK_GC_MAX = 0.75


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class GBlockRegion:
    """A region marked as gBlock, in amino acid coordinates (0-indexed, inclusive)."""
    start_aa: int
    end_aa: int
    label: str = ""

    @property
    def start_bp(self) -> int:
        return self.start_aa * 3

    @property
    def end_bp(self) -> int:
        return (self.end_aa + 1) * 3

    @property
    def len_aa(self) -> int:
        return self.end_aa - self.start_aa + 1

    @property
    def len_bp(self) -> int:
        return self.len_aa * 3

    def valid_size(self) -> bool:
        return GBLOCK_MIN_BP <= self.len_bp <= GBLOCK_MAX_BP


@dataclass
class PipelineResult:
    protein: str
    dna: str
    gblock_regions: list[GBlockRegion]
    beam_report: dict
    segments: list[dict]  # [{type: "oligo_tiled"|"gblock", start_bp, end_bp, dna, ...}]
    success: bool


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_protein(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "".join(l.strip() for l in lines if l and not l.startswith(">"))


def load_raw_table(path: Path) -> dict[str, dict[str, float]]:
    return parse_codon_table_raw(path.read_text())


def reverse_translate(protein: str, raw_table: dict) -> str:
    return "".join(
        max(raw_table[aa].items(), key=lambda kv: kv[1])[0]
        for aa in protein
    )


def gc_frac(seq: str) -> float:
    if not seq:
        return 0
    return (seq.count("G") + seq.count("C")) / len(seq)


def max_run(seq: str, bases: str) -> int:
    best = cur = 0
    for ch in seq:
        if ch in bases:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


# ---------------------------------------------------------------------------
# Step 1: gBlock-aware codon optimization
#
# Strategy: run the beam independently on each segment (gBlock vs non-gBlock),
# using appropriate constraints for each. This avoids the complexity of a
# single beam with a per-codon mask and is more robust — each segment gets
# the full beam width.
# ---------------------------------------------------------------------------

def codon_optimize_with_gblocks(
    protein: str,
    raw_table: dict,
    gblock_regions: list[GBlockRegion],
    oligo_constraints: dict,
) -> tuple[str, dict]:
    """Run per-segment codon optimization.

    Returns (optimized_dna, report).
    """
    n = len(protein)

    # Build sorted, non-overlapping segment list
    gblocks = sorted(gblock_regions, key=lambda g: g.start_aa)
    segments = []
    cursor = 0
    for gb in gblocks:
        if cursor < gb.start_aa:
            segments.append(("oligo", cursor, gb.start_aa - 1))
        segments.append(("gblock", gb.start_aa, gb.end_aa))
        cursor = gb.end_aa + 1
    if cursor < n:
        segments.append(("oligo", cursor, n - 1))

    # Optimize each segment
    dna_parts = []
    total_loss = 0.0
    ideal_total = 0.0
    all_succeeded = True
    seg_reports = []

    for kind, start, end in segments:
        sub_protein = protein[start:end + 1]
        sub_dna_naive = reverse_translate(sub_protein, raw_table)

        if kind == "gblock":
            kwargs = dict(
                gc_min=GBLOCK_GC_MIN,
                gc_max=GBLOCK_GC_MAX,
                gc_window=oligo_constraints.get("gc_window", 50),
                homo_gc=GBLOCK_HOMO_GC,
                homo_at=GBLOCK_HOMO_AT,
                avoid_patterns=oligo_constraints.get("avoid_patterns", []),
                min_codon_frequency=oligo_constraints.get("min_codon_frequency", 10.0),
                beam_k=oligo_constraints.get("beam_k", 250),
            )
        else:
            kwargs = dict(
                gc_min=oligo_constraints.get("gc_min", 0.25),
                gc_max=oligo_constraints.get("gc_max", 0.75),
                gc_window=oligo_constraints.get("gc_window", 50),
                homo_gc=oligo_constraints.get("homo_gc", 4),
                homo_at=oligo_constraints.get("homo_at", 6),
                avoid_patterns=oligo_constraints.get("avoid_patterns", []),
                min_codon_frequency=oligo_constraints.get("min_codon_frequency", 10.0),
                beam_k=oligo_constraints.get("beam_k", 250),
            )

        refined, _warnings, report = refine_gc_windows(
            sub_dna_naive, raw_table, **kwargs
        )
        dna_parts.append(refined)
        loss = report.get("total_loss", 0.0)
        ideal = report.get("ideal_score", 0.0)
        total_loss += loss
        ideal_total += ideal
        source = report.get("source", "?")
        if source != "beam":
            all_succeeded = False
        seg_reports.append({
            "kind": kind,
            "range_aa": f"{start}-{end}",
            "len_aa": end - start + 1,
            "source": source,
            "loss": round(loss, 1),
        })

    full_dna = "".join(dna_parts)
    return full_dna, {
        "success": all_succeeded,
        "total_loss": round(total_loss, 1),
        "ideal_score": round(ideal_total, 1),
        "segments": seg_reports,
    }


# ---------------------------------------------------------------------------
# Step 2: gBlock-aware oligo design
#
# Split the DNA at gBlock boundaries. Tile the non-gBlock regions as oligos.
# Emit gBlock regions as whole fragments. At boundaries, extend each piece
# by `overlap_length` bp into the adjacent piece so the assembly has proper
# overlaps.
# ---------------------------------------------------------------------------

def design_oligos_with_gblocks(
    dna: str,
    gblock_regions: list[GBlockRegion],
    oligo_length: int = 45,
    overlap_length: int = 20,
) -> list[dict]:
    """Design oligos + gBlock fragments with proper boundary overlaps.

    Returns a list of segments, each either:
      {"type": "gblock", "label": ..., "dna": ..., "start": ..., "end": ...,
       "len_bp": ..., "gc": ..., "cost_estimate": ...}
    or
      {"type": "oligo_region", "start": ..., "end": ..., "n_oligos": ...,
       "oligos": [...]}
    """
    n = len(dna)
    gblocks = sorted(gblock_regions, key=lambda g: g.start_bp)

    # Build intervals: list of (type, start_bp, end_bp) — non-overlapping, covering [0, n)
    intervals = []
    cursor = 0
    for gb in gblocks:
        if cursor < gb.start_bp:
            intervals.append(("oligo", cursor, gb.start_bp))
        intervals.append(("gblock", gb.start_bp, gb.end_bp, gb))
        cursor = gb.end_bp
    if cursor < n:
        intervals.append(("oligo", cursor, n))

    output = []
    for iv in intervals:
        if iv[0] == "gblock":
            _, start, end, gb = iv
            # Extend the gBlock by overlap_length on each side for assembly overlap,
            # clamped to sequence boundaries.
            ext_start = max(0, start - overlap_length)
            ext_end = min(n, end + overlap_length)
            frag_dna = dna[ext_start:ext_end]
            output.append({
                "type": "gblock",
                "label": gb.label or f"gBlock {start}-{end}",
                "core_start": start,
                "core_end": end,
                "ext_start": ext_start,
                "ext_end": ext_end,
                "dna": frag_dna,
                "len_bp": len(frag_dna),
                "gc": round(gc_frac(frag_dna), 3),
                "max_g_run": max_run(frag_dna, "G"),
                "max_c_run": max_run(frag_dna, "C"),
                "valid_size": GBLOCK_MIN_BP <= len(frag_dna) <= GBLOCK_MAX_BP,
                "cost_estimate_usd": round(len(frag_dna) * 0.07, 2),
            })
        else:
            _, start, end = iv
            region_dna = dna[start:end]
            if len(region_dna) == 0:
                continue
            # Simple tiling (no shift optimization — just demonstrate the concept)
            step = oligo_length - overlap_length
            oligos = []
            pos = 0
            idx = 0
            while pos < len(region_dna):
                olig_end = min(pos + oligo_length, len(region_dna))
                olig_seq = region_dna[pos:olig_end]
                strand = "sense" if idx % 2 == 0 else "antisense"
                oligos.append({
                    "index": idx,
                    "seq": olig_seq,
                    "start": start + pos,
                    "end": start + olig_end,
                    "length": len(olig_seq),
                    "strand": strand,
                    "gc": round(gc_frac(olig_seq), 3),
                })
                pos += step
                idx += 1
            output.append({
                "type": "oligo_region",
                "start": start,
                "end": end,
                "len_bp": end - start,
                "n_oligos": len(oligos),
                "oligos": oligos,
            })

    return output


# ---------------------------------------------------------------------------
# Full pipeline
# ---------------------------------------------------------------------------

def run_pipeline(
    protein: str,
    raw_table: dict,
    gblock_regions: list[GBlockRegion],
    oligo_constraints: dict | None = None,
    oligo_length: int = 45,
    overlap_length: int = 20,
) -> PipelineResult:
    if oligo_constraints is None:
        oligo_constraints = {
            "gc_min": 0.25, "gc_max": 0.75, "gc_window": 50,
            "homo_gc": 4, "homo_at": 6,
            "avoid_patterns": [],
            "min_codon_frequency": 10.0,
            "beam_k": 250,
        }

    # Validate gBlock sizes
    for gb in gblock_regions:
        if not gb.valid_size():
            print(f"  WARNING: {gb.label or 'gBlock'} is {gb.len_bp} bp "
                  f"(IDT requires {GBLOCK_MIN_BP}-{GBLOCK_MAX_BP} bp)")

    # Step 1: codon optimize
    dna, opt_report = codon_optimize_with_gblocks(
        protein, raw_table, gblock_regions, oligo_constraints
    )

    # Verify back-translation
    back = str(Seq(dna).translate()).rstrip("*")
    prot_clean = protein.rstrip("*")
    if back != prot_clean:
        print("  ERROR: back-translation mismatch!")

    # Step 2: design oligos + gBlock fragments
    segments = design_oligos_with_gblocks(
        dna, gblock_regions, oligo_length, overlap_length
    )

    return PipelineResult(
        protein=protein,
        dna=dna,
        gblock_regions=gblock_regions,
        beam_report=opt_report,
        segments=segments,
        success=opt_report["success"],
    )


def print_result(result: PipelineResult):
    print(f"\n{'=' * 72}")
    print(f"PIPELINE RESULT  (success={result.success})")
    print(f"{'=' * 72}")
    print(f"  protein: {len(result.protein)} aa")
    print(f"  DNA:     {len(result.dna)} bp")
    print(f"  back-translation OK: {str(Seq(result.dna).translate()).rstrip('*') == result.protein.rstrip('*')}")
    print(f"\n  Codon optimization:")
    print(f"    total_loss: {result.beam_report['total_loss']}")
    print(f"    ideal_score: {result.beam_report['ideal_score']}")
    for seg in result.beam_report["segments"]:
        status = "OK" if seg["source"] == "beam" else f"FAIL ({seg['source']})"
        print(f"    [{seg['kind']:>6}] aa {seg['range_aa']:<10} "
              f"({seg['len_aa']:>3} aa)  {status}  loss={seg['loss']}")

    print(f"\n  Assembly plan:")
    total_oligos = 0
    total_gblocks = 0
    total_cost = 0.0
    for seg in result.segments:
        if seg["type"] == "gblock":
            total_gblocks += 1
            total_cost += seg["cost_estimate_usd"]
            size_ok = "OK" if seg["valid_size"] else "TOO SMALL/LARGE"
            print(f"    [gBlock]  {seg['label']:<30}  {seg['len_bp']:>4} bp  "
                  f"GC={seg['gc']:.0%}  maxG={seg['max_g_run']}  maxC={seg['max_c_run']}  "
                  f"~${seg['cost_estimate_usd']:.2f}  {size_ok}")
        else:
            total_oligos += seg["n_oligos"]
            print(f"    [oligos]  bp {seg['start']}-{seg['end']:<5}  "
                  f"{seg['len_bp']:>4} bp  {seg['n_oligos']} oligos")

    print(f"\n  Totals: {total_oligos} oligos + {total_gblocks} gBlocks  "
          f"(gBlock cost ~${total_cost:.2f})")


# ---------------------------------------------------------------------------
# Demo: S1297 with various gBlock markings
# ---------------------------------------------------------------------------

def main():
    protein_path = (
        PROJECT_ROOT / "data" / "example_input_data_2026_03_27"
        / "glp_1_constructs_v2_S1297.protein.fasta"
    )
    table_path = PROJECT_ROOT / "data" / "HumColi_CodonTable.txt"

    protein = load_protein(protein_path)
    raw_table = load_raw_table(table_path)

    print(f"S1297 protein: {len(protein)} aa, {len(protein) * 3} bp")

    # Domain annotations for S1297
    # CD5(0-23) | linker(24-33) | Furin(34-37) | exenatide(38-77)
    # | linker(78-97) | AgeI(98-99) | Fc(100-331) | linker(332-341)
    # | TEV(342-348) | linker(349-358) | BamH1(359-361) | scFv(362-636)

    configs = [
        (
            "CONFIG A: Fc as gBlock only (auto-detect showed 2 hotspots in Fc)",
            [GBlockRegion(100, 331, "Fc")],
        ),
        (
            "CONFIG B: Fc + scFv as gBlocks",
            [GBlockRegion(100, 331, "Fc"), GBlockRegion(362, 636, "scFv")],
        ),
        (
            "CONFIG C: just the three auto-detected hotspots (may be too small for IDT)",
            [
                GBlockRegion(130, 155, "Fc hotspot 1"),    # padded to ~78 bp
                GBlockRegion(265, 295, "Fc hotspot 2"),    # padded to ~93 bp
                GBlockRegion(455, 500, "scFv hotspot+linker"),  # padded to ~138 bp
            ],
        ),
        (
            "CONFIG D: hotspots padded to IDT minimum (42+ aa = 126+ bp)",
            [
                GBlockRegion(120, 165, "Fc region 1"),     # 46 aa = 138 bp
                GBlockRegion(255, 300, "Fc region 2"),     # 46 aa = 138 bp
                GBlockRegion(450, 510, "scFv region"),     # 61 aa = 183 bp
            ],
        ),
    ]

    for label, gblocks in configs:
        print(f"\n{'#' * 72}")
        print(f"# {label}")
        print(f"{'#' * 72}")
        for gb in gblocks:
            print(f"  {gb.label}: aa {gb.start_aa}-{gb.end_aa} "
                  f"({gb.len_aa} aa, {gb.len_bp} bp, "
                  f"valid IDT size: {gb.valid_size()})")

        t0 = time.time()
        result = run_pipeline(protein, raw_table, gblocks)
        elapsed = time.time() - t0
        print(f"  elapsed: {elapsed:.1f}s")
        print_result(result)


if __name__ == "__main__":
    main()
