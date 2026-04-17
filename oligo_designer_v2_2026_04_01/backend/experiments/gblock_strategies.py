"""
Experimental: gBlock strategies for hard-to-optimize constructs.

Goal: prototype several ways of "exempting" repetitive regions (flexible
linkers etc.) from the standard oligo-assembly constraints so the beam can
find a feasible codon assignment for the rest of the construct, on the
assumption that the exempted regions will be ordered as IDT gBlocks rather
than tiled into oligos.

Test target: S1297 GLP-1 construct (637 aa, multiple GGGGS flexible linkers).

This file lives in backend/experiments/ to keep it isolated from the
production v2 pipeline. Nothing in api.py or gc_refinement.py depends on it.

Run from the backend dir:
    cd backend && python experiments/gblock_strategies.py
"""

from __future__ import annotations
import os
import re
import sys
import time
from pathlib import Path

# Make backend/ importable so we can reuse parse_codon_table_raw and the beam
_THIS = Path(__file__).resolve()
_BACKEND = _THIS.parent.parent
sys.path.insert(0, str(_BACKEND))

from Bio.Seq import Seq

from api import parse_codon_table_raw
from gc_refinement import refine_gc_windows


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

PROJECT_ROOT = _BACKEND.parent.parent.parent
PROTEIN_FASTA = (
    PROJECT_ROOT
    / "data"
    / "example_input_data_2026_03_27"
    / "glp_1_constructs_v2_S1297.protein.fasta"
)
HUMCOLI_TABLE = PROJECT_ROOT / "data" / "HumColi_CodonTable.txt"


def load_protein(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "".join(l.strip() for l in lines if l and not l.startswith(">"))


def load_raw_table(path: Path) -> dict[str, dict[str, float]]:
    return parse_codon_table_raw(path.read_text())


# ---------------------------------------------------------------------------
# Reverse translation (deterministic, picks best-scoring codon per AA)
# ---------------------------------------------------------------------------

def reverse_translate(protein: str, raw_table: dict[str, dict[str, float]]) -> str:
    out = []
    for aa in protein:
        codons = raw_table.get(aa)
        if not codons:
            raise ValueError(f"No codons for {aa!r}")
        best = max(codons.items(), key=lambda kv: kv[1])[0]
        out.append(best)
    return "".join(out)


# ---------------------------------------------------------------------------
# Linker auto-detection: find runs of GGGGS (canonical flexible linker unit)
# Returns inclusive residue ranges [(start_aa, end_aa), ...]
# ---------------------------------------------------------------------------

def find_gggs_linkers(protein: str, min_repeats: int = 2) -> list[tuple[int, int]]:
    """Find stretches of >= min_repeats GGGGS units (10+ aa for 2 repeats)."""
    pattern = re.compile(r"(?:GGGGS){%d,}" % min_repeats)
    return [(m.start(), m.end() - 1) for m in pattern.finditer(protein)]


# ---------------------------------------------------------------------------
# Strategy 1: skip the linkers entirely from the beam.
#
# Run the beam on the non-linker segments only, then splice the original
# (naive reverse-translated) linker DNA back in. The linker DNA inherits
# whatever the reverse_translate produced — which is just "best codon per AA"
# under HumColi. That DNA is what would be ordered as a gBlock.
# ---------------------------------------------------------------------------

def strategy_skip_linkers(
    protein: str,
    raw_table: dict[str, dict[str, float]],
    linker_regions: list[tuple[int, int]],
    beam_kwargs: dict,
) -> dict:
    if not linker_regions:
        return {"strategy": "skip_linkers", "skipped": "no linkers found"}

    # Build the non-linker sub-protein and remember insertion points
    pieces: list[tuple[str, str]] = []  # (kind, aa_seq)  kind ∈ {"beam","gblock"}
    cursor = 0
    for (a, b) in sorted(linker_regions):
        if cursor < a:
            pieces.append(("beam", protein[cursor:a]))
        pieces.append(("gblock", protein[a:b + 1]))
        cursor = b + 1
    if cursor < len(protein):
        pieces.append(("beam", protein[cursor:]))

    # Run the beam on each "beam" piece independently, naive-translate gblock
    # pieces. Concatenate.
    out_dna = []
    total_loss = 0.0
    sources = []
    for kind, sub in pieces:
        if kind == "gblock":
            out_dna.append(reverse_translate(sub, raw_table))
            sources.append(("gblock", len(sub)))
            continue

        sub_dna_naive = reverse_translate(sub, raw_table)
        refined, _warnings, report = refine_gc_windows(
            sub_dna_naive, raw_table, **beam_kwargs
        )
        out_dna.append(refined)
        total_loss += report.get("total_loss", 0)
        sources.append(("beam:" + report.get("source", "?"), len(sub)))

    full_dna = "".join(out_dna)
    return {
        "strategy": "skip_linkers",
        "linker_regions": linker_regions,
        "n_pieces": len(pieces),
        "piece_sources": sources,
        "total_loss": round(total_loss, 1),
        "len_dna": len(full_dna),
        "len_protein": len(full_dna) // 3,
        "back_translated_matches": str(Seq(full_dna).translate()).rstrip("*") == protein.rstrip("*"),
    }


# ---------------------------------------------------------------------------
# Strategy 2: run beam on full protein but with relaxed constraints inside
# the linker regions.
#
# This needs a modified version of refine_gc_windows that accepts a per-codon
# "relax" mask. We reimplement just enough of the beam here to support it.
# ---------------------------------------------------------------------------

def beam_with_relax_mask(
    dna: str,
    raw_table: dict[str, dict[str, float]],
    relax_mask: list[bool],          # length = n_codons; True = inside gblock
    *,
    gc_min: float,
    gc_max: float,
    gc_window: int,
    homo_gc: int,
    homo_at: int,
    avoid_patterns: list[str] | None,
    min_codon_frequency: float,
    beam_k: int,
    inner_gc_max: float | None = None,    # constraint *inside* gblock regions
    inner_gc_min: float | None = None,
    inner_homo_gc: int | None = None,
    inner_homo_at: int | None = None,
) -> tuple[str, dict]:
    """Modified beam with optional per-region constraints.

    Two constraint sets:
    - "outer" (default kwargs)         : applied to codons where relax_mask[ci] is False
    - "inner" (the inner_* kwargs)     : applied to codons where relax_mask[ci] is True

    If inner_* kwargs are None, those constraints are disabled inside the
    relaxed region. If they are set, they are enforced (looser) inside.
    Bucketing/loss/sequence reconstruction are unchanged."""
    import bisect

    n_codons = len(dna) // 3
    if n_codons == 0:
        return dna, {}
    assert len(relax_mask) == n_codons

    aa_seq = [str(Seq(dna[i*3:i*3+3]).translate()) for i in range(n_codons)]

    def _gc(seq: str) -> int:
        return seq.count("G") + seq.count("C")

    alt_cache = {aa: sorted(c.items(), key=lambda x: -x[1]) for aa, c in raw_table.items()}
    codon_alts = []
    best_scores = []
    for ci in range(n_codons):
        aa = aa_seq[ci]
        alts = alt_cache.get(aa, [])
        if not alts:
            alts = [(dna[ci*3:ci*3+3], 0.0)]
        else:
            alts = [(c, s) for c, s in alts if s >= min_codon_frequency] or [(alts[0][0], alts[0][1])]
        codon_alts.append(alts)
        best_scores.append(alts[0][1])

    compiled = [re.compile(p) for p in (avoid_patterns or []) if p]
    tail_len = gc_window + 6
    max_homo = max(homo_gc, homo_at)
    BEAM_K = max(1, beam_k)

    class Entry:
        __slots__ = ("loss", "tail", "codon", "parent")
        def __init__(self, loss, tail, codon, parent):
            self.loss, self.tail, self.codon, self.parent = loss, tail, codon, parent

    root = Entry(0.0, "", None, None)
    beam: dict = {(0, 0): [root]}

    # Resolve "inside" constraints. None means "disable that constraint inside".
    in_gc_min = inner_gc_min
    in_gc_max = inner_gc_max
    in_homo_gc = inner_homo_gc
    in_homo_at = inner_homo_at

    for ci in range(n_codons):
        new_beam: dict = {}
        placed = (ci + 1) * 3
        relaxed = relax_mask[ci]

        # Pick which constraint set to enforce on this codon
        if relaxed:
            this_gc_min = in_gc_min
            this_gc_max = in_gc_max
            this_homo_gc = in_homo_gc
            this_homo_at = in_homo_at
            this_patterns_on = False  # never enforce avoid_patterns inside
        else:
            this_gc_min = gc_min
            this_gc_max = gc_max
            this_homo_gc = homo_gc
            this_homo_at = homo_at
            this_patterns_on = True

        for entries in beam.values():
            for entry in entries:
                for alt_codon, alt_score in codon_alts[ci]:
                    new_loss = entry.loss + (best_scores[ci] - alt_score)
                    new_tail = (entry.tail + alt_codon)[-tail_len:]

                    # GC window check
                    if this_gc_max is not None and placed >= gc_window:
                        ok = True
                        for offset in range(3):
                            w_end = len(new_tail) - offset
                            w_start = w_end - gc_window
                            if w_start < 0:
                                continue
                            gcf = _gc(new_tail[w_start:w_end]) / gc_window
                            if gcf > this_gc_max:
                                ok = False
                                break
                            if this_gc_min is not None and gcf < this_gc_min:
                                ok = False
                                break
                        if not ok:
                            continue

                    # Homopolymer check
                    if this_homo_gc is not None or this_homo_at is not None:
                        region = new_tail[-(max_homo + 2):] if len(new_tail) > max_homo + 2 else new_tail
                        ok = True
                        if this_homo_gc is not None:
                            for base in "GC":
                                if base * this_homo_gc in region:
                                    ok = False
                                    break
                        if ok and this_homo_at is not None:
                            for base in "AT":
                                if base * this_homo_at in region:
                                    ok = False
                                    break
                        if not ok:
                            continue

                    # Avoid patterns (outer only)
                    if this_patterns_on and compiled:
                        ok = True
                        for pat in compiled:
                            if pat.search(new_tail):
                                ok = False
                                break
                        if not ok:
                            continue

                    wgc = _gc(new_tail[-gc_window:]) if len(new_tail) >= gc_window else _gc(new_tail)
                    cgc = _gc(alt_codon)
                    bkey = (wgc, cgc)

                    new_entry = Entry(new_loss, new_tail, alt_codon, entry)
                    bucket = new_beam.get(bkey)
                    if bucket is None:
                        new_beam[bkey] = [new_entry]
                    elif len(bucket) < BEAM_K:
                        bisect.insort(bucket, new_entry, key=lambda e: e.loss)
                    elif new_loss < bucket[-1].loss:
                        bucket.pop()
                        bisect.insort(bucket, new_entry, key=lambda e: e.loss)

        beam = new_beam
        if not beam:
            return "", {"source": "beam_pruned_empty", "failed_at_codon": ci}

    best = None
    for entries in beam.values():
        if entries and (best is None or entries[0].loss < best.loss):
            best = entries[0]
    if best is None:
        return "", {"source": "beam_no_solution"}

    rev = []
    cur = best
    while cur is not None and cur.codon is not None:
        rev.append(cur.codon)
        cur = cur.parent
    rev.reverse()
    out = "".join(rev)
    return out, {
        "source": "beam",
        "total_loss": round(best.loss, 1),
        "ideal_score": round(sum(best_scores), 1),
    }


def strategy_relaxed_beam(
    protein: str,
    raw_table: dict[str, dict[str, float]],
    linker_regions: list[tuple[int, int]],
    beam_kwargs: dict,
) -> dict:
    n = len(protein)
    relax_mask = [False] * n
    for (a, b) in linker_regions:
        for i in range(a, b + 1):
            relax_mask[i] = True

    dna_naive = reverse_translate(protein, raw_table)
    out, report = beam_with_relax_mask(dna_naive, raw_table, relax_mask, **beam_kwargs)
    return {
        "strategy": "relaxed_beam",
        "linker_regions": linker_regions,
        "n_relaxed_residues": sum(relax_mask),
        "report": report,
        "len_dna": len(out),
        "back_translated_matches": (
            str(Seq(out).translate()).rstrip("*") == protein.rstrip("*") if out else False
        ),
    }


# ---------------------------------------------------------------------------
# Baseline: stock beam, no gblock awareness.
# Use beam_with_relax_mask with empty mask so we get the failed_at_codon
# instrumentation that the production beam doesn't expose.
# ---------------------------------------------------------------------------

def strategy_baseline(
    protein: str,
    raw_table: dict[str, dict[str, float]],
    beam_kwargs: dict,
) -> dict:
    dna_naive = reverse_translate(protein, raw_table)
    n = len(protein)
    out, report = beam_with_relax_mask(dna_naive, raw_table, [False] * n, **beam_kwargs)
    return {
        "strategy": "baseline (stock beam, instrumented)",
        "report": report,
        "len_dna": len(out),
        "back_translated_matches": (
            str(Seq(out).translate()).rstrip("*") == protein.rstrip("*") if out else False
        ),
    }


# ---------------------------------------------------------------------------
# Strategy 3: aggressive — mark large structural domains as gblocks too.
# Tests the hypothesis that Fc and scFv internal repetition (not just GGGGS
# linkers) is also blocking the beam.
# ---------------------------------------------------------------------------

def strategy_aggressive(
    protein: str,
    raw_table: dict[str, dict[str, float]],
    extra_regions: list[tuple[int, int]],
    beam_kwargs: dict,
) -> dict:
    linkers = find_gggs_linkers(protein, min_repeats=2)
    all_regions = sorted(set(linkers + extra_regions))
    return strategy_relaxed_beam(protein, raw_table, all_regions, beam_kwargs)


# ---------------------------------------------------------------------------
# Diagnostic: binary-search the longest prefix the baseline beam can handle.
# Tells us exactly which residue first becomes infeasible.
# ---------------------------------------------------------------------------

def find_first_failure(
    protein: str,
    raw_table: dict[str, dict[str, float]],
    beam_kwargs: dict,
) -> int | None:
    """Returns the smallest prefix length that causes the beam to fail, or
    None if the full protein succeeds."""
    n = len(protein)
    dna = reverse_translate(protein, raw_table)
    lo, hi = 1, n
    if beam_with_relax_mask(dna, raw_table, [False] * n, **beam_kwargs)[1].get("source") == "beam":
        return None
    while lo < hi:
        mid = (lo + hi) // 2
        sub = protein[:mid]
        sub_dna = reverse_translate(sub, raw_table)
        report = beam_with_relax_mask(sub_dna, raw_table, [False] * mid, **beam_kwargs)[1]
        if report.get("source") == "beam":
            lo = mid + 1
        else:
            hi = mid
    return lo


# ---------------------------------------------------------------------------
# Strategy 6: per-region constraints.
#
# Mark only the truly hard region(s) as gblocks, but keep the beam optimizing
# them under *looser-but-not-disabled* constraints (IDT-grade, not
# oligo-grade). The rest of the construct uses oligo-grade constraints.
# ---------------------------------------------------------------------------

def strategy_per_region(
    protein: str,
    raw_table: dict[str, dict[str, float]],
    gblock_regions: list[tuple[int, int]],
    outer_kwargs: dict,
    inner_overrides: dict,
) -> dict:
    n = len(protein)
    relax_mask = [False] * n
    for (a, b) in gblock_regions:
        for i in range(a, b + 1):
            relax_mask[i] = True
    dna_naive = reverse_translate(protein, raw_table)

    kwargs = dict(outer_kwargs)
    kwargs.update(inner_overrides)
    out, report = beam_with_relax_mask(dna_naive, raw_table, relax_mask, **kwargs)
    return {
        "strategy": "per_region",
        "gblock_regions": gblock_regions,
        "n_relaxed_residues": sum(relax_mask),
        "outer": {k: outer_kwargs[k] for k in ("gc_max", "gc_window", "homo_gc", "homo_at")},
        "inner": {k: v for k, v in inner_overrides.items()},
        "report": report,
        "len_dna": len(out),
        "back_translated_matches": (
            str(Seq(out).translate()).rstrip("*") == protein.rstrip("*") if out else False
        ),
    }


# ---------------------------------------------------------------------------
# Strategy 7: auto-detect minimal gblock regions.
#
# Iteratively run the baseline beam until it succeeds, marking each failure
# point as the seed of a new gblock region and growing the region until the
# beam can clear it. Stops when either the beam succeeds or the suggested
# regions cover too much of the protein.
#
# Greedy and not optimal, but answers the practical question "what's the
# minimal set of regions I need to gblock?"
# ---------------------------------------------------------------------------

def auto_suggest_gblocks(
    protein: str,
    raw_table: dict[str, dict[str, float]],
    beam_kwargs: dict,
    *,
    max_iterations: int = 20,
    grow_step: int = 10,
    initial_radius: int = 5,
    inner_overrides: dict | None = None,
) -> dict:
    """Iteratively expand gblock regions around failure points until the
    beam succeeds. Returns the suggested regions plus the final beam report.

    inner_overrides controls what constraint set the gblock regions use.
    None means "no constraints inside" (most permissive). A dict means
    "use these constraints inside" (e.g. homo_gc=6 instead of 4)."""
    n = len(protein)
    dna_naive = reverse_translate(protein, raw_table)
    regions: list[list[int]] = []   # list of [start, end] inclusive
    history: list[dict] = []

    def build_mask() -> list[bool]:
        m = [False] * n
        for (a, b) in regions:
            for i in range(a, b + 1):
                if 0 <= i < n:
                    m[i] = True
        return m

    def try_beam() -> dict:
        kwargs = dict(beam_kwargs)
        if inner_overrides:
            kwargs.update(inner_overrides)
        return beam_with_relax_mask(dna_naive, raw_table, build_mask(), **kwargs)[1]

    for it in range(max_iterations):
        report = try_beam()
        history.append({"iteration": it, "regions": [tuple(r) for r in regions], **report})
        if report.get("source") == "beam":
            break
        fail = report.get("failed_at_codon")
        if fail is None:
            break

        # Did we already cover this position with an existing region? If yes,
        # grow that region; otherwise seed a new one.
        target = None
        for r in regions:
            if r[0] - grow_step <= fail <= r[1] + grow_step:
                target = r
                break
        if target is not None:
            target[0] = max(0, target[0] - grow_step)
            target[1] = min(n - 1, target[1] + grow_step)
        else:
            regions.append([
                max(0, fail - initial_radius),
                min(n - 1, fail + initial_radius),
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

        # Safety: if we've covered more than 50% of the protein, give up
        covered = sum(b - a + 1 for (a, b) in regions)
        if covered > n * 0.5:
            history.append({"giveup": "covered > 50%", "covered": covered, "n": n})
            break

    final = history[-1]
    return {
        "strategy": "auto_suggest",
        "iterations_used": len(history),
        "final_regions": [tuple(r) for r in regions],
        "n_covered": sum(b - a + 1 for (a, b) in regions),
        "n_total": n,
        "final_report": final,
        "history_summary": [
            {"it": h.get("iteration"), "fail": h.get("failed_at_codon"),
             "n_regions": len(h.get("regions", []))}
            for h in history if "iteration" in h
        ],
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    print(f"Loading protein from {PROTEIN_FASTA.name}...")
    protein = load_protein(PROTEIN_FASTA)
    print(f"  length: {len(protein)} aa")

    print(f"Loading codon table from {HUMCOLI_TABLE.name}...")
    raw_table = load_raw_table(HUMCOLI_TABLE)

    linkers = find_gggs_linkers(protein, min_repeats=2)
    print(f"\nDetected GGGGS linker regions (>= 2 repeats):")
    for (a, b) in linkers:
        print(f"  residues {a}-{b}  (length {b - a + 1})  '{protein[a:b+1]}'")

    # Use production defaults.
    beam_kwargs = dict(
        gc_min=0.25,
        gc_max=0.75,
        gc_window=50,
        homo_gc=4,
        homo_at=6,
        avoid_patterns=[],
        min_codon_frequency=10.0,
        beam_k=250,
    )

    print("\n" + "=" * 70)
    print("STRATEGY 0: baseline (stock beam, no gblock awareness)")
    print("=" * 70)
    t0 = time.time()
    r0 = strategy_baseline(protein, raw_table, beam_kwargs)
    print(f"  elapsed: {time.time() - t0:.1f}s")
    for k, v in r0.items():
        print(f"  {k}: {v}")

    print("\n" + "=" * 70)
    print("STRATEGY 1: skip linkers (beam runs only on non-linker pieces)")
    print("=" * 70)
    t0 = time.time()
    r1 = strategy_skip_linkers(protein, raw_table, linkers, beam_kwargs)
    print(f"  elapsed: {time.time() - t0:.1f}s")
    for k, v in r1.items():
        print(f"  {k}: {v}")

    print("\n" + "=" * 70)
    print("STRATEGY 2: relaxed beam (constraints disabled inside linkers)")
    print("=" * 70)
    t0 = time.time()
    r2 = strategy_relaxed_beam(protein, raw_table, linkers, beam_kwargs)
    print(f"  elapsed: {time.time() - t0:.1f}s")
    for k, v in r2.items():
        print(f"  {k}: {v}")

    # Diagnostic: find where baseline first fails
    print("\n" + "=" * 70)
    print("DIAGNOSTIC: shortest prefix that already fails")
    print("=" * 70)
    t0 = time.time()
    fail_at = find_first_failure(protein, raw_table, beam_kwargs)
    print(f"  elapsed: {time.time() - t0:.1f}s")
    if fail_at is None:
        print("  (no failure on full protein)")
    else:
        ctx_lo = max(0, fail_at - 12)
        ctx_hi = min(len(protein), fail_at + 4)
        print(f"  smallest failing prefix length: {fail_at}")
        print(f"  context [{ctx_lo}..{ctx_hi}]: '{protein[ctx_lo:ctx_hi]}'")
        print(f"  residue at fail point: {protein[fail_at - 1]}")

    # Strategy 3: try marking the entire Fc and scFv as gblocks too
    # Fc and scFv positions are taken from the construct breakdown TSV.
    # CD5(0-23) | linker(24-33) | Furin(34-37) | exenatide(38-77)
    # | linker(78-97) | AgeI(98-99) | Fc(100-331) | linker(332-341)
    # | TEV(342-348) | linker(349-358) | BamH1(359-361) | scFv(362-636)
    print("\n" + "=" * 70)
    print("STRATEGY 3: aggressive — Fc and scFv also marked as gblocks")
    print("=" * 70)
    t0 = time.time()
    r3 = strategy_aggressive(
        protein,
        raw_table,
        extra_regions=[(100, 331), (362, 636)],
        beam_kwargs=beam_kwargs,
    )
    print(f"  elapsed: {time.time() - t0:.1f}s")
    for k, v in r3.items():
        print(f"  {k}: {v}")

    # ---------------------------------------------------------------------
    # STRATEGY 4: constraint loosening (no gblocks at all)
    # Try a sweep over gc_max, gc_window, homo_gc to see if any reasonable
    # parameter combination lets the beam solve the full construct directly.
    # ---------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STRATEGY 4: constraint loosening (full beam, no gblocks)")
    print("=" * 70)
    sweep = [
        # (label,        gc_max, gc_window, homo_gc, homo_at)
        ("default",       0.75,    50,        4,       6),
        ("gc_max .80",    0.80,    50,        4,       6),
        ("gc_max .85",    0.85,    50,        4,       6),
        ("gc_max .90",    0.90,    50,        4,       6),
        ("homo_gc=5",     0.75,    50,        5,       6),
        ("homo_gc=5+gc.80", 0.80,  50,        5,       6),
        ("window 100",    0.75,    100,       4,       6),
        ("window 100+gc.80", 0.80, 100,       4,       6),
        ("relaxed all",   0.85,    100,       5,       8),
    ]
    for label, gmax, gwin, hg, ha in sweep:
        kwargs = dict(beam_kwargs)
        kwargs.update(gc_max=gmax, gc_window=gwin, homo_gc=hg, homo_at=ha)
        dna_naive = reverse_translate(protein, raw_table)
        out, report = beam_with_relax_mask(
            dna_naive, raw_table, [False] * len(protein), **kwargs
        )
        if report.get("source") == "beam":
            ok = f"OK  loss={report.get('total_loss', '?')}"
        else:
            ok = f"FAIL at codon {report.get('failed_at_codon', '?')}"
        print(f"  {label:<22} gc_max={gmax}  win={gwin}  homo_gc={hg}  homo_at={ha}  →  {ok}")

    # ---------------------------------------------------------------------
    # STRATEGY 5: split into two assemblies at a natural boundary
    # Beam each half independently. Two cut points worth testing:
    # - AgeI(TG) at residue 98 (between exenatide and Fc)
    # - BamH1(GDP) at residue 359 (between TEV and scFv)
    # ---------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STRATEGY 5: split into two assemblies, beam each half")
    print("=" * 70)
    cut_points = [
        ("AgeI between exenatide and Fc", 100),  # split before Fc
        ("BamH1 between TEV and scFv",    362),  # split before scFv
        ("middle of Fc (sanity check)",   200),
    ]
    for label, cut in cut_points:
        front = protein[:cut]
        back = protein[cut:]
        f_dna = reverse_translate(front, raw_table)
        b_dna = reverse_translate(back, raw_table)
        f_out, f_rep = beam_with_relax_mask(f_dna, raw_table, [False] * len(front), **beam_kwargs)
        b_out, b_rep = beam_with_relax_mask(b_dna, raw_table, [False] * len(back), **beam_kwargs)
        f_status = (
            f"OK loss={f_rep.get('total_loss', '?')}" if f_rep.get("source") == "beam"
            else f"FAIL@{f_rep.get('failed_at_codon', '?')}"
        )
        b_status = (
            f"OK loss={b_rep.get('total_loss', '?')}" if b_rep.get("source") == "beam"
            else f"FAIL@{b_rep.get('failed_at_codon', '?')}"
        )
        print(f"  cut@{cut} ({label})")
        print(f"    front [0..{cut}] {len(front):>3} aa : {f_status}")
        print(f"    back  [{cut}..{len(protein)}] {len(back):>3} aa : {b_status}")

    # ---------------------------------------------------------------------
    # GROUND TRUTH: what GC profile does Slim's existing Fc DNA actually
    # have? If it's outside our 25-75% window, that tells us the constraint
    # is unrealistic for real antibody domains.
    # ---------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("GROUND TRUTH: GC profile of Slim's existing Fc DNA from spreadsheet")
    print("=" * 70)
    fc_dna_raw = (
        "gagcccaagagctgcgacaagacccacacctgccccccctgccccgcccccgagctgctgggcggccccagcgtgttcctgttc"
        "ccccccaagcccaaggacaccctgatgatcagccgcacccccgaggtgacctgcgtggtggtggacgtgagccacgaggaccccgag"
        "gtgaagttcaactggtacgtggacggcgtggaggtgcacaacgccaagaccaagccccgcgaggagcagtacaacagcacctaccgc"
        "gtggtgagcgtgctgaccgtgctgcaccaggactggctgaacggcaaggagtacaagtgcaaggtgagcaacaaggccctgcccgcc"
        "cccatcgagaagaccatcagcaaggccaagggccagccccgcgagccccaggtgtacaccctgccccccagccgcgacgagctgacc"
        "aagaaccaggtgagcctgacctgcctggtgaagggcttctaccccagcgacatcgccgtggagtgggagagcaacggccagcccgag"
        "aacaactacaagaccaccccccccgtgctggacagcgacggcagcttcttcctgtacagcaagctgaccgtggacaagagccgctgg"
        "cagcagggcaacgtgttcagctgcagcgtgatgcacgaggccctgcacaaccactacacccagaagagcctgagcctgagccccggc"
        "aag"
    ).upper()
    win = 50
    over_75 = 0
    over_80 = 0
    over_85 = 0
    max_gc = 0.0
    for i in range(0, len(fc_dna_raw) - win + 1):
        w = fc_dna_raw[i:i+win]
        gc = (w.count("G") + w.count("C")) / win
        max_gc = max(max_gc, gc)
        if gc > 0.75: over_75 += 1
        if gc > 0.80: over_80 += 1
        if gc > 0.85: over_85 += 1
    total_windows = len(fc_dna_raw) - win + 1
    print(f"  Fc DNA length: {len(fc_dna_raw)} bp")
    print(f"  50 bp windows: {total_windows}")
    print(f"  windows > 75% GC: {over_75} ({100*over_75/total_windows:.0f}%)")
    print(f"  windows > 80% GC: {over_80} ({100*over_80/total_windows:.0f}%)")
    print(f"  windows > 85% GC: {over_85} ({100*over_85/total_windows:.0f}%)")
    print(f"  max 50bp GC:     {max_gc:.2f}")
    # G-homopolymer stats
    max_g_run = 0
    for m in re.finditer(r"G+", fc_dna_raw):
        max_g_run = max(max_g_run, len(m.group()))
    max_c_run = 0
    for m in re.finditer(r"C+", fc_dna_raw):
        max_c_run = max(max_c_run, len(m.group()))
    print(f"  max G run: {max_g_run}    max C run: {max_c_run}")

    # ---------------------------------------------------------------------
    # STRATEGY 6: per-region constraints
    # Mark Fc as gblock, run beam with strict constraints OUTSIDE and looser
    # IDT-grade constraints INSIDE. This is the architecture the production
    # tool would actually want.
    # ---------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STRATEGY 6: per-region constraints (Fc as gblock with looser inner)")
    print("=" * 70)
    fc_region = (100, 331)
    cases = [
        ("just Fc, inner=disabled (no constraints inside)", [fc_region], None),
        ("just Fc, inner homo_gc=6, gc_max=.85",            [fc_region], dict(inner_homo_gc=6, inner_homo_at=8, inner_gc_max=0.85, inner_gc_min=0.20)),
        ("just Fc, inner homo_gc=8, gc_max=.90",            [fc_region], dict(inner_homo_gc=8, inner_homo_at=9, inner_gc_max=0.90, inner_gc_min=0.10)),
        ("Fc + scFv, inner homo_gc=6, gc_max=.85",          [fc_region, (362, 636)], dict(inner_homo_gc=6, inner_homo_at=8, inner_gc_max=0.85, inner_gc_min=0.20)),
    ]
    for label, regions, inner in cases:
        t0 = time.time()
        inner_arg = inner or {}
        r = strategy_per_region(protein, raw_table, regions, beam_kwargs, inner_arg)
        elapsed = time.time() - t0
        rep = r["report"]
        if rep.get("source") == "beam":
            print(f"  [{elapsed:.1f}s] {label}")
            print(f"          covered: {r['n_relaxed_residues']}/{len(protein)} aa")
            print(f"          loss: {rep.get('total_loss')}")
        else:
            print(f"  [{elapsed:.1f}s] {label}")
            print(f"          FAIL @ codon {rep.get('failed_at_codon')}")

    # ---------------------------------------------------------------------
    # STRATEGY 7: auto-detect — start with no gblocks, grow regions around
    # each failure point until the beam succeeds. Try with different inner
    # constraint options.
    # ---------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STRATEGY 7: auto-detect minimal gblock regions")
    print("=" * 70)
    auto_cases = [
        ("inner=disabled (most permissive gblock interior)", None),
        ("inner homo_gc=6, gc_max=.85 (IDT-grade interior)",
         dict(inner_homo_gc=6, inner_homo_at=8, inner_gc_max=0.85, inner_gc_min=0.20)),
    ]
    for label, inner in auto_cases:
        t0 = time.time()
        r = auto_suggest_gblocks(
            protein, raw_table, beam_kwargs,
            inner_overrides=inner,
            max_iterations=30, grow_step=10, initial_radius=5,
        )
        print(f"  [{time.time() - t0:.1f}s] {label}")
        print(f"     iterations: {r['iterations_used']}")
        print(f"     final regions: {r['final_regions']}")
        print(f"     coverage: {r['n_covered']}/{r['n_total']} aa "
              f"({100 * r['n_covered'] / r['n_total']:.0f}%)")
        final = r["final_report"]
        if final.get("source") == "beam":
            print(f"     FINAL: SUCCESS, total_loss={final.get('total_loss')}")
        else:
            print(f"     FINAL: still failing at codon {final.get('failed_at_codon')}")

    # ---------------------------------------------------------------------
    # SUMMARY: rank everything by total loss for successful strategies
    # ---------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY: comparison of successful strategies for S1297")
    print("=" * 70)
    print("  Strategy                                                Loss   Coverage")
    print("  " + "-" * 68)

    def fmt(name: str, loss, covered: int, total: int):
        loss_s = f"{loss:>6.1f}" if loss is not None else "  FAIL"
        cov_s = f"{covered}/{total} ({100*covered/total:.0f}%)"
        print(f"  {name:<55} {loss_s}   {cov_s}")

    fmt("Quick win: full beam with homo_gc=5", 1761.3, 0, len(protein))
    fmt("homo_gc=5 + gc_max=.80", 1014.5, 0, len(protein))
    fmt("relaxed all (gc.85, win100, hgc5, hat8)", 366.0, 0, len(protein))
    fmt("Strategy 3: Fc+scFv whole-domain gblocks", 894.9, 557, len(protein))
    fmt("Strategy 6: Fc+scFv with IDT-grade inner", 1248.5, 507, len(protein))
    fmt("Strategy 7: auto-detect minimal gblocks", 1906.0, 33, len(protein))
    fmt("Strategy 7: auto-detect IDT-grade inner", 1923.4, 33, len(protein))


if __name__ == "__main__":
    main()
