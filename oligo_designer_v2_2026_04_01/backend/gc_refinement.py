"""
Post-optimization GC refinement pass using dynamic programming.

DNA Chisel's greedy left-to-right solver can make suboptimal codon choices
when resolving GC-content window constraints. It discovers violations at
the end of a window and is forced to make large sacrifices on the last 1-2
codons instead of spreading the cost across many codons.

This module re-examines each GC-constrained window and uses DP to find the
minimum-cost set of codon swaps that satisfies the constraint. "Cost" = total
loss in codon frequency (/1000 score) relative to the best codon for each
amino acid.
"""

from __future__ import annotations
from Bio.Seq import Seq


def _gc_count(seq: str) -> int:
    return sum(1 for b in seq if b in "GC")


def _find_suboptimal_windows(dna: str, raw_table: dict, gc_max: float,
                              gc_window: int) -> list[tuple[int, int]]:
    """Find windows containing downgraded codons where GC is near the limit.

    Returns (window_start, window_end) for windows where:
    - At least one codon was downgraded from its best, AND
    - GC content is within 10% of the limit (meaning GC pressure caused it)
    """
    n_codons = len(dna) // 3
    regions = []

    # Find codon positions where the current codon is not the best
    downgraded = set()
    for ci in range(n_codons):
        codon = dna[ci*3:ci*3+3]
        aa = str(Seq(codon).translate())
        alts = raw_table.get(aa, {})
        if not alts:
            continue
        best_score = max(alts.values())
        cur_score = alts.get(codon, 0)
        if cur_score < best_score * 0.9:  # significantly downgraded
            downgraded.add(ci)

    if not downgraded:
        return []

    # For each downgraded codon, check if it's in a high-GC window
    for ci in sorted(downgraded):
        nt_pos = ci * 3
        # Check windows that contain this codon
        for ws in range(max(0, nt_pos - gc_window + 3), min(len(dna) - gc_window + 1, nt_pos + 1)):
            w = dna[ws:ws + gc_window]
            gc_frac = _gc_count(w) / gc_window
            # Window is near the GC max — this codon was likely downgraded for GC reasons
            if gc_frac >= gc_max - 0.10:
                regions.append((ws, ws + gc_window))
                break

    # Merge overlapping regions
    if not regions:
        return []
    regions.sort()
    merged = [regions[0]]
    for s, e in regions[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def refine_gc_windows(
    dna: str,
    raw_table: dict[str, dict[str, float]],
    gc_min: float,
    gc_max: float,
    gc_window: int,
    homo_gc: int,
    homo_at: int,
    avoid_patterns: list[str] | None = None,
    max_iterations: int = 10,
) -> str:
    """Refine codon choices to minimize cost of GC constraint satisfaction.

    Runs after DNA Chisel. For each GC-constrained window where codons were
    downgraded, uses DP to find the minimum total score loss that keeps GC%
    within bounds.
    """
    if len(dna) < gc_window:
        return dna

    result = list(dna)
    n_codons = len(dna) // 3
    aa_seq = [str(Seq(dna[i*3:i*3+3]).translate()) for i in range(n_codons)]

    # Build alternatives: aa -> [(codon, per_thousand)] sorted by score desc
    alt_cache: dict[str, list[tuple[str, float]]] = {}
    for aa, codons in raw_table.items():
        alt_cache[aa] = sorted(codons.items(), key=lambda x: -x[1])

    # Find regions to refine (windows with downgraded codons near GC limit)
    seq_str = "".join(result)
    regions = _find_suboptimal_windows(seq_str, raw_table, gc_max, gc_window)

    # Also add any remaining GC violations
    for i in range(len(seq_str) - gc_window + 1):
        w = seq_str[i:i + gc_window]
        gc = _gc_count(w) / gc_window
        if gc < gc_min or gc > gc_max:
            regions.append((i, i + gc_window))

    if not regions:
        return dna

    # Merge all regions
    regions.sort()
    merged = [regions[0]]
    for s, e in regions[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))

    for region_start, region_end in merged:
        seq_str = "".join(result)

        # Codons overlapping this region
        c_start = max(0, region_start // 3)
        c_end = min(n_codons, (region_end + 2) // 3)
        nuc_start = c_start * 3
        nuc_end = c_end * 3
        num_codons = c_end - c_start

        if num_codons == 0:
            continue

        # Get alternatives for each codon in the region
        codon_alts: list[list[tuple[str, float]]] = []
        for ci in range(c_start, c_end):
            aa = aa_seq[ci]
            alts = alt_cache.get(aa, [])
            if not alts:
                current = "".join(result[ci*3:ci*3+3])
                alts = [(current, 0.0)]
            codon_alts.append(alts)

        # For each codon, compute the best possible score
        best_scores = []
        for ci_local in range(num_codons):
            best_scores.append(codon_alts[ci_local][0][1])

        # DP: minimize total loss (sum of best_score - chosen_score)
        # State: dp[gc_count] = (min_total_loss, choices)
        max_gc = num_codons * 3
        INF = float("inf")
        dp: list[tuple[float, list[str]]] = [(INF, []) for _ in range(max_gc + 1)]
        dp[0] = (0.0, [])

        for ci_local in range(num_codons):
            new_dp: list[tuple[float, list[str]]] = [(INF, []) for _ in range(max_gc + 1)]

            for prev_gc in range(max_gc + 1):
                if dp[prev_gc][0] >= INF:
                    continue
                prev_loss, prev_choices = dp[prev_gc]

                for alt_codon, alt_score in codon_alts[ci_local]:
                    alt_gc = _gc_count(alt_codon)
                    new_gc = prev_gc + alt_gc
                    if new_gc > max_gc:
                        continue

                    loss = best_scores[ci_local] - alt_score
                    total_loss = prev_loss + loss

                    if total_loss < new_dp[new_gc][0]:
                        new_dp[new_gc] = (total_loss, prev_choices + [alt_codon])

            dp = new_dp

        # Find the best solution that satisfies all GC windows + other constraints
        # Get current total loss for comparison
        current_loss = 0.0
        for ci_local in range(num_codons):
            ci = c_start + ci_local
            current_codon = "".join(result[ci*3:ci*3+3])
            current_score = raw_table.get(aa_seq[ci], {}).get(current_codon, 0)
            current_loss += best_scores[ci_local] - current_score

        best_solution_loss = current_loss
        best_solution_choices: list[str] | None = None

        for gc_total in range(max_gc + 1):
            if dp[gc_total][0] >= INF:
                continue
            loss, choices = dp[gc_total]

            # Only consider solutions that are strictly better
            if loss >= best_solution_loss:
                continue

            # Build candidate and validate
            candidate = list(result)
            for ci_local, codon in enumerate(choices):
                ci = c_start + ci_local
                for j, base in enumerate(codon):
                    candidate[ci*3 + j] = base

            cand_str = "".join(candidate)

            # Check ALL GC windows overlapping the modified region
            check_start = max(0, nuc_start - gc_window + 1)
            check_end = min(len(dna) - gc_window + 1, nuc_end)
            gc_ok = True
            for ws in range(check_start, check_end):
                w = cand_str[ws:ws + gc_window]
                gc_frac = _gc_count(w) / gc_window
                if gc_frac < gc_min or gc_frac > gc_max:
                    gc_ok = False
                    break
            if not gc_ok:
                continue

            # Check homopolymer constraints
            context_start = max(0, nuc_start - max(homo_gc, homo_at))
            context_end = min(len(dna), nuc_end + max(homo_gc, homo_at))
            context = cand_str[context_start:context_end]
            homo_ok = True
            for base in "GC":
                if base * homo_gc in context:
                    homo_ok = False
                    break
            if homo_ok:
                for base in "AT":
                    if base * homo_at in context:
                        homo_ok = False
                        break
            if not homo_ok:
                continue

            # Check avoid_patterns
            if avoid_patterns:
                import re
                pattern_ok = True
                for pat in avoid_patterns:
                    if re.search(pat, cand_str[nuc_start:nuc_end]):
                        pattern_ok = False
                        break
                if not pattern_ok:
                    continue

            best_solution_loss = loss
            best_solution_choices = choices

        # Apply if we found something better
        if best_solution_choices is not None:
            for ci_local, codon in enumerate(best_solution_choices):
                ci = c_start + ci_local
                for j, base in enumerate(codon):
                    result[ci*3 + j] = base

    return "".join(result)
