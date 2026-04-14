"""
Global codon optimization via beam search with incremental constraint enforcement.

Given a starting DNA encoding the target protein, this module searches the full
synonymous codon space for a sequence that minimizes total frequency loss while
enforcing GC windows, homopolymers, avoid patterns, and minimum codon frequency
as hard constraints during the search — not as post-hoc validation.

The beam is bucketed by (window_gc, last_codon_gc) to preserve spatial diversity,
so paths whose GC totals match but whose tail distributions differ survive in
different buckets.
"""

from __future__ import annotations
import bisect
import re
from Bio.Seq import Seq


def _gc(seq: str) -> int:
    return seq.count("G") + seq.count("C")


def refine_gc_windows(
    dna: str,
    raw_table: dict[str, dict[str, float]],
    gc_min: float,
    gc_max: float,
    gc_window: int,
    homo_gc: int,
    homo_at: int,
    avoid_patterns: list[str] | None = None,
    min_codon_frequency: float = 10.0,
    beam_k: int = 15,
    progress_callback=None,
) -> tuple[str, list[dict], dict]:
    """Refine codon choices globally to minimize total frequency loss.

    All constraints are enforced during the search:
    - GC content within [gc_min, gc_max] for every sliding window
    - No G/C homopolymer runs >= homo_gc
    - No A/T homopolymer runs >= homo_at
    - No matches for avoid_patterns
    - No codon with frequency < min_codon_frequency (when alternatives exist)

    Returns (refined_dna, freq_warnings, optimization_report).
    """
    n_codons = len(dna) // 3
    if n_codons == 0:
        return dna, [], {}

    aa_seq = [str(Seq(dna[i*3:i*3+3]).translate()) for i in range(n_codons)]

    # Build alternatives per amino acid, sorted by score descending
    alt_cache: dict[str, list[tuple[str, float]]] = {}
    for aa, codons in raw_table.items():
        alt_cache[aa] = sorted(codons.items(), key=lambda x: -x[1])

    # Per-position alternatives, filtered by min frequency. The API-layer
    # feasibility check has already verified that every AA has at least one
    # codon above the floor, so the filtered list is guaranteed non-empty.
    codon_alts: list[list[tuple[str, float]]] = []
    best_scores: list[float] = []
    for ci in range(n_codons):
        aa = aa_seq[ci]
        alts = alt_cache.get(aa, [])
        if not alts:
            alts = [(dna[ci*3:ci*3+3], 0.0)]
        else:
            alts = [(c, s) for c, s in alts if s >= min_codon_frequency]
        codon_alts.append(alts)
        best_scores.append(alts[0][1])

    compiled_patterns = []
    if avoid_patterns:
        for pat in avoid_patterns:
            if pat:
                compiled_patterns.append(re.compile(pat))

    # Beam search
    # Each entry: (loss, tail, codon, parent)
    #   tail: last tail_len bases for incremental constraint checking
    #   parent: pointer to previous entry for sequence reconstruction
    # Beam bucketed by window GC count for spatial diversity.

    tail_len = gc_window + 6
    max_homo = max(homo_gc, homo_at)
    BEAM_K = max(1, beam_k)  # entries per bucket

    class Entry:
        __slots__ = ("loss", "tail", "codon", "parent")
        def __init__(self, loss: float, tail: str, codon: str | None, parent: "Entry | None"):
            self.loss = loss
            self.tail = tail
            self.codon = codon
            self.parent = parent

    root = Entry(0.0, "", None, None)
    # Bucket key: (window_gc_count, last_codon_gc_count). The second dimension
    # diversifies paths whose GC totals match but whose tails differ — those
    # have different futures because the next sliding window will see different
    # bases falling out.
    beam: dict[tuple[int, int], list[Entry]] = {(0, 0): [root]}

    for ci in range(n_codons):
        if progress_callback is not None:
            progress_callback(ci, n_codons)
        new_beam: dict[tuple[int, int], list[Entry]] = {}
        placed = (ci + 1) * 3

        for entries in beam.values():
            for entry in entries:
                for alt_codon, alt_score in codon_alts[ci]:
                    new_loss = entry.loss + (best_scores[ci] - alt_score)
                    new_tail = (entry.tail + alt_codon)[-tail_len:]

                    # GC window check: up to 3 newly complete windows
                    if placed >= gc_window:
                        gc_ok = True
                        for offset in range(3):
                            w_end = len(new_tail) - offset
                            w_start = w_end - gc_window
                            if w_start < 0:
                                continue
                            gc_frac = _gc(new_tail[w_start:w_end]) / gc_window
                            if gc_frac < gc_min or gc_frac > gc_max:
                                gc_ok = False
                                break
                        if not gc_ok:
                            continue

                    # Homopolymer check at boundary
                    region = new_tail[-(max_homo + 2):] if len(new_tail) > max_homo + 2 else new_tail
                    homo_ok = True
                    for base in "GC":
                        if base * homo_gc in region:
                            homo_ok = False
                            break
                    if homo_ok:
                        for base in "AT":
                            if base * homo_at in region:
                                homo_ok = False
                                break
                    if not homo_ok:
                        continue

                    # Avoid-pattern check at boundary
                    if compiled_patterns:
                        pat_ok = True
                        for pat in compiled_patterns:
                            if pat.search(new_tail):
                                pat_ok = False
                                break
                        if not pat_ok:
                            continue

                    # Bucket by (window GC, last-codon GC) for diversity
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
            return dna, _freq_warnings(dna, n_codons, aa_seq, raw_table, alt_cache, min_codon_frequency, beam_failed=True), _early_report(
                "beam_pruned_empty", BEAM_K, dna, n_codons, aa_seq, raw_table, best_scores,
            )

    # Best solution across all buckets
    best: Entry | None = None
    for entries in beam.values():
        if entries and (best is None or entries[0].loss < best.loss):
            best = entries[0]

    if best is None:
        return dna, _freq_warnings(dna, n_codons, aa_seq, raw_table, alt_cache, min_codon_frequency, beam_failed=True), _early_report(
            "beam_no_solution", BEAM_K, dna, n_codons, aa_seq, raw_table, best_scores,
        )

    # Reconstruct sequence from parent chain
    codons_rev: list[str] = []
    curr: Entry | None = best
    while curr is not None and curr.codon is not None:
        codons_rev.append(curr.codon)
        curr = curr.parent
    codons_rev.reverse()
    result_dna = "".join(codons_rev)

    ideal_score = sum(best_scores)
    chosen_score = sum(
        raw_table.get(aa_seq[ci], {}).get(result_dna[ci*3:ci*3+3], 0)
        for ci in range(n_codons)
    )
    report = {
        "ideal_score": round(ideal_score, 1),
        "chosen_score": round(chosen_score, 1),
        "total_loss": round(ideal_score - chosen_score, 1),
        "source": "beam",
        "beam_k": BEAM_K,
    }

    return result_dna, _freq_warnings(result_dna, n_codons, aa_seq, raw_table, alt_cache, min_codon_frequency), report


def _early_report(source: str, beam_k: int, dna: str, n_codons: int,
                  aa_seq: list, raw_table: dict, best_scores: list) -> dict:
    """Partial report for early-return paths (beam failed). Scores reflect the
    naive reverse-translation we fall back to, not a real optimized output."""
    ideal_score = sum(best_scores)
    chosen_score = sum(
        raw_table.get(aa_seq[ci], {}).get(dna[ci*3:ci*3+3], 0)
        for ci in range(n_codons)
    )
    return {
        "source": source,
        "beam_k": beam_k,
        "ideal_score": round(ideal_score, 1),
        "chosen_score": round(chosen_score, 1),
        "total_loss": round(ideal_score - chosen_score, 1),
    }


def _freq_warnings(
    dna: str, n_codons: int, aa_seq: list[str],
    raw_table: dict, alt_cache: dict, min_codon_frequency: float,
    beam_failed: bool = False,
) -> list[dict]:
    """Warnings for codons still below the frequency floor.

    When beam_failed is True, the reason text reflects that no global
    solution exists rather than claiming a per-codon constraint conflict.
    """
    warnings: list[dict] = []
    if min_codon_frequency <= 0:
        return warnings
    for ci in range(n_codons):
        codon = dna[ci*3:ci*3+3]
        aa = aa_seq[ci]
        score = raw_table.get(aa, {}).get(codon, 0)
        if score < min_codon_frequency:
            has_better = any(s >= min_codon_frequency for _, s in alt_cache.get(aa, []))
            if beam_failed:
                reason = f"below floor ({min_codon_frequency:.0f}/1000) — no global solution exists under current constraints"
            elif has_better:
                reason = "other constraints (GC%, homopolymer, avoid pattern) prevent using a higher-frequency codon"
            else:
                reason = f"no {aa} codon in this table has frequency >= {min_codon_frequency}/1000"
            warnings.append({
                "kind": "low_codon_frequency",
                "message": f"{codon} ({aa}) at pos {ci*3}: {score:.1f}/1000 — {reason}",
                "start": ci * 3,
                "end": ci * 3 + 3,
                "codon_index": ci,
                "frequency": score,
            })
    return warnings
