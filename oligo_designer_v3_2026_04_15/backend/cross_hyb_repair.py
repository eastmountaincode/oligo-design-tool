"""
Cross-hybridization repair pass.

Runs AFTER codon optimization and oligo tiling. When
`check_cross_hybridization` flags one or more overlap pairs whose
heterodimer Tm exceeds `avg_overlap_tm - 20°C`, this module searches
for single-codon (and, if needed, paired) synonymous substitutions
inside the flagged overlap regions that drop the heterodimer below
threshold without introducing any new constraint violations.

Why run this as a separate stage
--------------------------------
The heterodimer threshold is ``avg_overlap_tm - 20°C``, and
``avg_overlap_tm`` is a property of the tiled design — it can only be
computed after codon optimization AND tiling are complete. That
circular dependency rules out expressing cross-hybridization as a
codon-level constraint. Instead we let the threshold settle after the
first tile pass (where it's stable to within a fraction of a degree
under single-codon edits), then use it as the evaluation target for a
bounded, auditable repair loop.

Design notes
------------
* Cut positions are FIXED during repair. We only change bases, never
  the tile layout. This guarantees the repair is a minimum perturbation
  of the original design and keeps the threshold calculation stable.
* Every candidate swap is validated against the SAME structural
  constraints used during codon optimization (GC window, homopolymer
  limits, avoid-patterns) evaluated on the FULL construct
  (flank + insert + flank). Swaps that would introduce a new violation
  are rejected.
* Synonyms are drawn from the user's codon table and must stay at or
  above `min_codon_frequency` (per-thousand).
* Codons inside gBlock regions are never touched — gBlocks are fixed
  segments the user has committed to.
* Ranking: among swaps that fix at least one cross-hyb pair without
  adding a new one, prefer the one that minimizes |ΔGC content|, then
  the one with the smallest codon-frequency drop. This matches the
  "do the least damage" intuition: keep the sequence as close to the
  originally-optimized version as possible.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from oligo_designer import (
    Overlap,
    OverlapIssue,
    ReactionConditions,
    check_overlap,
    check_cross_hybridization,
    reverse_complement,
    gc_content,
)


# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------

@dataclass
class RepairLogEntry:
    """One step of the repair loop — either an applied swap or a failure note."""
    iteration: int
    kind: str                           # "applied" | "no_move_found" | "skipped"
    # For applied swaps:
    codon_index: Optional[int] = None   # 0-indexed codon position in the protein
    insert_position: Optional[int] = None  # base position in insert (codon_index * 3)
    amino_acid: Optional[str] = None
    old_codon: Optional[str] = None
    new_codon: Optional[str] = None
    old_freq_per_k: Optional[float] = None
    new_freq_per_k: Optional[float] = None
    gc_before: Optional[float] = None
    gc_after: Optional[float] = None
    fixed_pairs: list[tuple[int, int]] = field(default_factory=list)   # 1-indexed
    # For either:
    remaining_pairs: list[tuple[int, int]] = field(default_factory=list)  # 1-indexed
    remaining_max_tm: Optional[float] = None
    threshold_tm: Optional[float] = None
    notes: str = ""


@dataclass
class RepairResult:
    """Outcome of running `repair_cross_hybridization`."""
    new_insert_dna: str
    new_overlaps: list[Overlap]
    log: list[RepairLogEntry]
    threshold_tm: float
    issues_before: int
    issues_after: int


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def repair_cross_hybridization(
    *,
    insert_dna: str,
    protein: str,
    codon_freqs: dict[str, dict[str, float]],
    min_codon_frequency: float,
    overlaps: list[Overlap],
    plasmid_upstream: str,
    plasmid_downstream: str,
    conditions: ReactionConditions,
    gc_min: float,
    gc_max: float,
    gc_window: int,
    avoid_homopolymers_gc: int,
    avoid_homopolymers_at: int,
    avoid_patterns: list[str],
    gblock_ranges_insert: list[tuple[int, int]] | None = None,
    max_iterations: int = 20,
    try_pairs_if_single_fails: bool = True,
    pair_search_cap: int = 800,
) -> RepairResult:
    """Fix overlap cross-hybridization via minimum-cost synonymous codon swaps.

    Parameters
    ----------
    insert_dna
        The codon-optimized DNA for the insert (in-frame with `protein`).
    protein
        Amino acid sequence (single-letter), with no trailing '*'.
    codon_freqs
        Per-thousand codon usage table: ``{aa_1letter: {codon: per_k}}``.
    min_codon_frequency
        Per-thousand floor: synonyms below this are ineligible.
    overlaps
        Current tiled overlap list (from `tile_sequence` / `tile_sequence_with_gblocks`).
        Their `start` / `end` fields are treated as fixed during repair.
    plasmid_upstream / plasmid_downstream
        Flank sequences. The full construct used for constraint checks is
        ``plasmid_upstream + insert_dna + plasmid_downstream``.
    conditions
        Reaction conditions for thermodynamic calls.
    gc_min / gc_max / gc_window
        Structural GC-window constraint evaluated on the full construct.
    avoid_homopolymers_gc / avoid_homopolymers_at
        Homopolymer run limits evaluated on the full construct.
    avoid_patterns
        Extra regex patterns (e.g. restriction sites) evaluated on the full
        construct.
    gblock_ranges_insert
        Optional list of ``(start_bp, end_bp)`` ranges (in INSERT coordinates)
        describing gBlock segments. Codons whose bases overlap a gBlock range
        are never swapped.
    max_iterations
        Upper bound on repair loop iterations.
    try_pairs_if_single_fails
        If single-codon search finds no valid move for a flagged pair,
        try two-codon combinations before giving up.
    pair_search_cap
        Hard cap on the number of two-codon candidates evaluated per
        iteration. 800 ~ 0.5-1s on a typical machine.
    """
    thermo = conditions.make_thermo()

    # Work on copies — leave the caller's objects untouched.
    current_dna = insert_dna
    current_overlaps = [_clone_overlap(o) for o in overlaps]

    # Threshold is fixed based on the ORIGINAL tiling's avg Tm. Single-codon
    # swaps move this value by a fraction of a degree, so treating it as
    # constant across the repair loop is safe and makes comparison stable.
    overlap_tms = [o.tm for o in current_overlaps if o.tm is not None]
    if not overlap_tms:
        return RepairResult(
            new_insert_dna=current_dna,
            new_overlaps=current_overlaps,
            log=[RepairLogEntry(iteration=0, kind="skipped",
                                notes="No overlap Tms available; skipping repair.")],
            threshold_tm=float("nan"),
            issues_before=0,
            issues_after=0,
        )
    avg_tm = sum(overlap_tms) / len(overlap_tms)
    threshold_tm = avg_tm - 20.0

    issues_before = _count_cross_hyb_pairs(current_overlaps, thermo, threshold_tm)

    # Pre-compute: for each codon position, the list of valid synonyms
    # (codon -> per_k). Excludes codons inside gBlock regions and
    # codons whose AA has no eligible synonym.
    candidate_codons = _enumerate_candidate_codons(
        insert_dna=current_dna,
        protein=protein,
        codon_freqs=codon_freqs,
        min_codon_frequency=min_codon_frequency,
        gblock_ranges_insert=gblock_ranges_insert or [],
    )

    # Per-codon original frequency lookup (used for log + cost ranking).
    def _orig_freq(codon_index: int, current_codon: str) -> float:
        aa = protein[codon_index] if codon_index < len(protein) else ""
        return float(codon_freqs.get(aa, {}).get(current_codon.upper(), 0.0))

    log: list[RepairLogEntry] = []

    for it in range(1, max_iterations + 1):
        flagged = _find_flagged_pairs(current_overlaps, thermo, threshold_tm)
        if not flagged:
            break

        # Codon positions we're allowed to explore this iteration: any codon
        # whose bases touch ANY flagged overlap (in insert coords). Swapping
        # outside the flagged overlaps can't fix the pair (the bases we'd
        # change aren't in either overlap).
        target_positions = _codon_positions_touching_overlaps(
            overlaps=current_overlaps,
            flagged=flagged,
            insert_offset=len(plasmid_upstream),
            insert_len=len(current_dna),
            candidate_codons=candidate_codons,
        )

        # Single-codon search first.
        best = _best_single_codon_swap(
            current_dna=current_dna,
            current_overlaps=current_overlaps,
            target_positions=target_positions,
            candidate_codons=candidate_codons,
            orig_freq_fn=_orig_freq,
            flagged=flagged,
            threshold_tm=threshold_tm,
            thermo=thermo,
            plasmid_upstream=plasmid_upstream,
            plasmid_downstream=plasmid_downstream,
            gc_min=gc_min, gc_max=gc_max, gc_window=gc_window,
            homo_gc=avoid_homopolymers_gc,
            homo_at=avoid_homopolymers_at,
            avoid_patterns=avoid_patterns,
        )

        if best is None and try_pairs_if_single_fails:
            best = _best_two_codon_swap(
                current_dna=current_dna,
                current_overlaps=current_overlaps,
                target_positions=target_positions,
                candidate_codons=candidate_codons,
                orig_freq_fn=_orig_freq,
                flagged=flagged,
                threshold_tm=threshold_tm,
                thermo=thermo,
                plasmid_upstream=plasmid_upstream,
                plasmid_downstream=plasmid_downstream,
                gc_min=gc_min, gc_max=gc_max, gc_window=gc_window,
                homo_gc=avoid_homopolymers_gc,
                homo_at=avoid_homopolymers_at,
                avoid_patterns=avoid_patterns,
                search_cap=pair_search_cap,
            )

        if best is None:
            # Give up on the remaining pairs, record what's left.
            log.append(RepairLogEntry(
                iteration=it,
                kind="no_move_found",
                remaining_pairs=_pairs_as_1indexed(flagged, current_overlaps),
                remaining_max_tm=_max_hetero_tm(flagged, current_overlaps, thermo),
                threshold_tm=threshold_tm,
                notes=(f"Tried all eligible single-codon swaps"
                       + (" and a paired search" if try_pairs_if_single_fails else "")
                       + " within the flagged overlaps; none resolved the pair(s)"
                       " without breaking another constraint."),
            ))
            break

        # Apply the move.
        current_dna, current_overlaps = _apply_move(
            current_dna=current_dna,
            current_overlaps=current_overlaps,
            move=best,
            plasmid_upstream=plasmid_upstream,
            plasmid_downstream=plasmid_downstream,
            thermo=thermo,
            conditions=conditions,
            threshold_tm=threshold_tm,
        )

        log.append(_move_to_log_entry(
            iteration=it,
            move=best,
            current_overlaps=current_overlaps,
            thermo=thermo,
            threshold_tm=threshold_tm,
            protein=protein,
        ))

    issues_after = _count_cross_hyb_pairs(current_overlaps, thermo, threshold_tm)

    return RepairResult(
        new_insert_dna=current_dna,
        new_overlaps=current_overlaps,
        log=log,
        threshold_tm=threshold_tm,
        issues_before=issues_before,
        issues_after=issues_after,
    )


# ---------------------------------------------------------------------------
# Move representation
# ---------------------------------------------------------------------------

@dataclass
class _Move:
    """A candidate or applied swap."""
    codon_swaps: list[tuple[int, str, str, float, float]]
    # Each entry: (codon_index, old_codon, new_codon, old_freq, new_freq)
    new_dna: str
    new_overlaps: list[Overlap]
    fixed_pairs: list[tuple[int, int]]           # 0-indexed
    remaining_pairs: list[tuple[int, int]]       # 0-indexed
    gc_before: float
    gc_after: float
    freq_drop: float       # sum of (old_freq - new_freq), larger = worse
    remaining_max_tm: float


# ---------------------------------------------------------------------------
# Candidate enumeration
# ---------------------------------------------------------------------------

def _enumerate_candidate_codons(
    *,
    insert_dna: str,
    protein: str,
    codon_freqs: dict[str, dict[str, float]],
    min_codon_frequency: float,
    gblock_ranges_insert: list[tuple[int, int]],
) -> dict[int, list[tuple[str, float]]]:
    """For each codon index, list (codon, per_k) synonyms eligible for swap.

    Excludes the original codon itself. Excludes codons whose bases fall
    inside any gBlock range (gBlocks are fixed).
    """
    out: dict[int, list[tuple[str, float]]] = {}
    n_codons = min(len(protein), len(insert_dna) // 3)
    for ci in range(n_codons):
        bp_start = ci * 3
        bp_end = bp_start + 3
        if _range_overlaps_any(bp_start, bp_end, gblock_ranges_insert):
            continue
        aa = protein[ci]
        orig = insert_dna[bp_start:bp_end].upper()
        syns = codon_freqs.get(aa, {})
        candidates = [
            (codon, freq)
            for codon, freq in syns.items()
            if codon != orig and freq >= min_codon_frequency
        ]
        if candidates:
            out[ci] = candidates
    return out


def _codon_positions_touching_overlaps(
    *,
    overlaps: list[Overlap],
    flagged: list[tuple[int, int]],
    insert_offset: int,
    insert_len: int,
    candidate_codons: dict[int, list[tuple[str, float]]],
) -> list[int]:
    """Codon indices whose bases intersect ANY flagged overlap."""
    ranges_insert: set[tuple[int, int]] = set()
    for i, j in flagged:
        for idx in (i, j):
            o = overlaps[idx]
            s = max(0, o.start - insert_offset)
            e = min(insert_len, o.end - insert_offset)
            if e > s:
                ranges_insert.add((s, e))
    positions: set[int] = set()
    for s, e in ranges_insert:
        first = s // 3
        last = (e - 1) // 3
        for ci in range(first, last + 1):
            if ci in candidate_codons:
                positions.add(ci)
    return sorted(positions)


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

def _best_single_codon_swap(
    *,
    current_dna: str,
    current_overlaps: list[Overlap],
    target_positions: list[int],
    candidate_codons: dict[int, list[tuple[str, float]]],
    orig_freq_fn,
    flagged: list[tuple[int, int]],
    threshold_tm: float,
    thermo,
    plasmid_upstream: str,
    plasmid_downstream: str,
    gc_min: float, gc_max: float, gc_window: int,
    homo_gc: int, homo_at: int,
    avoid_patterns: list[str],
) -> Optional[_Move]:
    best: Optional[_Move] = None
    gc_before = gc_content(current_dna)
    for ci in target_positions:
        orig_codon = current_dna[ci * 3:ci * 3 + 3]
        orig_freq = orig_freq_fn(ci, orig_codon)
        for new_codon, new_freq in candidate_codons[ci]:
            move = _evaluate_swap(
                current_dna=current_dna,
                current_overlaps=current_overlaps,
                swaps=[(ci, orig_codon, new_codon, orig_freq, new_freq)],
                flagged=flagged,
                threshold_tm=threshold_tm,
                thermo=thermo,
                plasmid_upstream=plasmid_upstream,
                plasmid_downstream=plasmid_downstream,
                gc_min=gc_min, gc_max=gc_max, gc_window=gc_window,
                homo_gc=homo_gc, homo_at=homo_at,
                avoid_patterns=avoid_patterns,
                gc_before=gc_before,
            )
            if move is None:
                continue
            if best is None or _move_cost(move) < _move_cost(best):
                best = move
    return best


def _best_two_codon_swap(
    *,
    current_dna: str,
    current_overlaps: list[Overlap],
    target_positions: list[int],
    candidate_codons: dict[int, list[tuple[str, float]]],
    orig_freq_fn,
    flagged: list[tuple[int, int]],
    threshold_tm: float,
    thermo,
    plasmid_upstream: str,
    plasmid_downstream: str,
    gc_min: float, gc_max: float, gc_window: int,
    homo_gc: int, homo_at: int,
    avoid_patterns: list[str],
    search_cap: int,
) -> Optional[_Move]:
    """Bounded two-codon search. Budget-capped: caller sets the limit."""
    best: Optional[_Move] = None
    gc_before = gc_content(current_dna)
    per_pos_options: list[tuple[int, str, str, float, float]] = []
    for ci in target_positions:
        orig_codon = current_dna[ci * 3:ci * 3 + 3]
        orig_freq = orig_freq_fn(ci, orig_codon)
        for new_codon, new_freq in candidate_codons[ci]:
            per_pos_options.append((ci, orig_codon, new_codon, orig_freq, new_freq))
    n = len(per_pos_options)
    budget = search_cap
    for a in range(n):
        for b in range(a + 1, n):
            if per_pos_options[a][0] == per_pos_options[b][0]:
                # Same codon position — skip (already covered by single-codon search)
                continue
            budget -= 1
            if budget <= 0:
                return best
            move = _evaluate_swap(
                current_dna=current_dna,
                current_overlaps=current_overlaps,
                swaps=[per_pos_options[a], per_pos_options[b]],
                flagged=flagged,
                threshold_tm=threshold_tm,
                thermo=thermo,
                plasmid_upstream=plasmid_upstream,
                plasmid_downstream=plasmid_downstream,
                gc_min=gc_min, gc_max=gc_max, gc_window=gc_window,
                homo_gc=homo_gc, homo_at=homo_at,
                avoid_patterns=avoid_patterns,
                gc_before=gc_before,
            )
            if move is None:
                continue
            if best is None or _move_cost(move) < _move_cost(best):
                best = move
    return best


def _evaluate_swap(
    *,
    current_dna: str,
    current_overlaps: list[Overlap],
    swaps: list[tuple[int, str, str, float, float]],
    flagged: list[tuple[int, int]],
    threshold_tm: float,
    thermo,
    plasmid_upstream: str,
    plasmid_downstream: str,
    gc_min: float, gc_max: float, gc_window: int,
    homo_gc: int, homo_at: int,
    avoid_patterns: list[str],
    gc_before: float,
) -> Optional[_Move]:
    """Apply `swaps` to `current_dna` and score. Returns None if:
      * structural constraints fail after the swap, or
      * the swap doesn't strictly reduce the flagged-pair set.
    """
    new_dna = _apply_codon_swaps(current_dna, swaps)
    # Structural-constraint check on the full construct.
    full = plasmid_upstream.upper() + new_dna + plasmid_downstream.upper()
    ok, _reason = _passes_structural_constraints(
        full,
        gc_min=gc_min, gc_max=gc_max, gc_window=gc_window,
        homo_gc=homo_gc, homo_at=homo_at,
        avoid_patterns=avoid_patterns,
    )
    if not ok:
        return None

    # Recompute overlap sequences/Tms at SAME start/end positions.
    new_overlaps = _recompute_overlaps(
        overlaps=current_overlaps, full_seq=full, thermo=thermo,
    )

    # Cross-hyb pair set after swap.
    new_flagged = _find_flagged_pairs(new_overlaps, thermo, threshold_tm)

    # Must strictly reduce the flagged set (and not add anything new).
    old_set = set(flagged)
    new_set = set(new_flagged)
    fixed = sorted(old_set - new_set)
    appeared = new_set - old_set
    if not fixed:
        return None
    if appeared:
        return None

    gc_after = gc_content(new_dna)
    freq_drop = sum(max(0.0, old_f - new_f) for _ci, _oc, _nc, old_f, new_f in swaps)
    remaining_max = _max_hetero_tm(new_flagged, new_overlaps, thermo) if new_flagged else 0.0

    return _Move(
        codon_swaps=list(swaps),
        new_dna=new_dna,
        new_overlaps=new_overlaps,
        fixed_pairs=fixed,
        remaining_pairs=sorted(new_flagged),
        gc_before=gc_before,
        gc_after=gc_after,
        freq_drop=freq_drop,
        remaining_max_tm=remaining_max,
    )


def _move_cost(m: _Move) -> tuple:
    """Sort key: fewer remaining pairs, then smaller GC delta, then smaller freq drop."""
    return (
        len(m.remaining_pairs),
        round(abs(m.gc_after - m.gc_before), 4),
        round(m.freq_drop, 2),
        len(m.codon_swaps),
    )


# ---------------------------------------------------------------------------
# State transitions
# ---------------------------------------------------------------------------

def _apply_move(
    *,
    current_dna: str,
    current_overlaps: list[Overlap],
    move: _Move,
    plasmid_upstream: str,
    plasmid_downstream: str,
    thermo,
    conditions: ReactionConditions,
    threshold_tm: float,
) -> tuple[str, list[Overlap]]:
    """Apply an already-validated move and return fresh state."""
    return move.new_dna, move.new_overlaps


def _apply_codon_swaps(dna: str, swaps: list[tuple]) -> str:
    out = list(dna)
    for ci, _old, new_codon, *_rest in swaps:
        bp = ci * 3
        out[bp:bp + 3] = list(new_codon)
    return "".join(out)


def _recompute_overlaps(
    *,
    overlaps: list[Overlap],
    full_seq: str,
    thermo,
) -> list[Overlap]:
    """Rebuild overlaps in place at their existing start/end positions.

    The returned overlaps carry their intrinsic issues (homopolymer, GC, etc.)
    and Tm but do NOT carry cross-hybridization issues — those are inspected
    separately by `_find_flagged_pairs` during the search loop, and attached
    once at the end by the caller via `check_cross_hybridization`.
    """
    out = []
    for o in overlaps:
        new_seq = full_seq[o.start:o.end]
        new_ovl = Overlap(seq=new_seq, start=o.start, end=o.end)
        new_ovl.tm = thermo.calc_tm(new_seq)
        new_ovl.issues = check_overlap(new_ovl, full_seq)
        out.append(new_ovl)
    return out


# ---------------------------------------------------------------------------
# Cross-hyb helpers
# ---------------------------------------------------------------------------

def _find_flagged_pairs(
    overlaps: list[Overlap],
    thermo,
    threshold_tm: float,
) -> list[tuple[int, int]]:
    """Return 0-indexed (i,j) overlap pairs whose heterodimer Tm exceeds threshold."""
    flagged = []
    for i in range(len(overlaps)):
        for j in range(i + 1, len(overlaps)):
            si = overlaps[i].seq.upper()
            sj = overlaps[j].seq.upper()
            rcj = reverse_complement(sj)
            hit = False
            for target in (sj, rcj):
                if thermo.calc_heterodimer(si, target).tm > threshold_tm:
                    hit = True
                    break
            if hit:
                flagged.append((i, j))
    return flagged


def _count_cross_hyb_pairs(overlaps, thermo, threshold_tm: float) -> int:
    return len(_find_flagged_pairs(overlaps, thermo, threshold_tm))


def _max_hetero_tm(
    pairs: list[tuple[int, int]],
    overlaps: list[Overlap],
    thermo,
) -> float:
    best = 0.0
    for i, j in pairs:
        si = overlaps[i].seq.upper()
        sj = overlaps[j].seq.upper()
        rcj = reverse_complement(sj)
        for target in (sj, rcj):
            tm = thermo.calc_heterodimer(si, target).tm
            if tm > best:
                best = tm
    return best


def _pairs_as_1indexed(
    pairs: list[tuple[int, int]],
    overlaps: list[Overlap],
) -> list[tuple[int, int]]:
    # Overlap "index" in the API response is i+1 of its list position;
    # we keep repair log pairs in the same 1-indexed form for the UI.
    return [(i + 1, j + 1) for (i, j) in pairs]


# ---------------------------------------------------------------------------
# Structural constraints
# ---------------------------------------------------------------------------

def _passes_structural_constraints(
    full_seq: str,
    *,
    gc_min: float, gc_max: float, gc_window: int,
    homo_gc: int, homo_at: int,
    avoid_patterns: list[str],
) -> tuple[bool, str]:
    """Fast local re-check of the constraints the beam enforced.

    We don't call DNA Chisel here — we only need a boolean pass/fail for
    candidate scoring, and the checks below are cheap to run in the inner
    loop. Semantics match `_build_constraints` in api.py.
    """
    import re
    s = full_seq.upper()

    # Homopolymers.
    for base, thr in [('G', homo_gc), ('C', homo_gc), ('A', homo_at), ('T', homo_at)]:
        if re.search(f'{base}{{{thr},}}', s):
            return False, f"homopolymer {base}x{thr}+"

    # Avoid patterns (interpreted as DNA Chisel would: substrings or
    # ambiguous-code patterns. We compile each as a plain regex after
    # translating IUPAC ambiguity codes.)
    for pat in avoid_patterns or []:
        if not pat:
            continue
        rx = _iupac_to_regex(pat.upper())
        if re.search(rx, s):
            return False, f"avoid pattern {pat}"
        rc = reverse_complement(pat.upper())
        rx_rc = _iupac_to_regex(rc)
        if re.search(rx_rc, s):
            return False, f"avoid pattern {pat} (rc)"

    # GC sliding window.
    if len(s) >= gc_window and gc_window > 0:
        # Count walk to avoid O(len*window).
        win = gc_window
        gc_count = sum(1 for c in s[:win] if c in "GC")
        lo = gc_count / win
        hi = gc_count / win
        for i in range(win, len(s)):
            if s[i] in "GC":
                gc_count += 1
            if s[i - win] in "GC":
                gc_count -= 1
            frac = gc_count / win
            if frac < lo:
                lo = frac
            if frac > hi:
                hi = frac
        if lo < gc_min or hi > gc_max:
            return False, f"GC window {lo:.2f}-{hi:.2f} outside [{gc_min},{gc_max}]"

    return True, ""


_IUPAC = {
    "A": "A", "T": "T", "C": "C", "G": "G", "U": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def _iupac_to_regex(pat: str) -> str:
    import re
    out = []
    for ch in pat:
        if ch in _IUPAC:
            out.append(_IUPAC[ch])
        else:
            # If the user supplied regex metacharacters, keep them literal-safe.
            out.append(re.escape(ch))
    return "".join(out)


# ---------------------------------------------------------------------------
# Misc
# ---------------------------------------------------------------------------

def _clone_overlap(o: Overlap) -> Overlap:
    new = Overlap(seq=o.seq, start=o.start, end=o.end)
    new.tm = o.tm
    new.issues = [OverlapIssue(i.kind, i.message, i.severity, i.start, i.end) for i in o.issues]
    return new


def _range_overlaps_any(
    start: int, end: int, ranges: list[tuple[int, int]]
) -> bool:
    for rs, re_ in ranges:
        if start < re_ and end > rs:
            return True
    return False


# ---------------------------------------------------------------------------
# Log formatting
# ---------------------------------------------------------------------------

def _move_to_log_entry(
    *,
    iteration: int,
    move: _Move,
    current_overlaps: list[Overlap],
    thermo,
    threshold_tm: float,
    protein: str = "",
) -> RepairLogEntry:
    if len(move.codon_swaps) == 1:
        ci, old_codon, new_codon, old_freq, new_freq = move.codon_swaps[0]
        return RepairLogEntry(
            iteration=iteration,
            kind="applied",
            codon_index=ci,
            insert_position=ci * 3,
            amino_acid=(protein[ci] if 0 <= ci < len(protein) else None),
            old_codon=old_codon,
            new_codon=new_codon,
            old_freq_per_k=old_freq,
            new_freq_per_k=new_freq,
            gc_before=round(move.gc_before, 4),
            gc_after=round(move.gc_after, 4),
            fixed_pairs=[(i + 1, j + 1) for (i, j) in move.fixed_pairs],
            remaining_pairs=[(i + 1, j + 1) for (i, j) in move.remaining_pairs],
            remaining_max_tm=round(move.remaining_max_tm, 1) if move.remaining_pairs else None,
            threshold_tm=round(threshold_tm, 1),
        )
    # Pair swap — flatten into a single log row with concatenated notes.
    parts = []
    for ci, old_codon, new_codon, _old_f, _new_f in move.codon_swaps:
        parts.append(f"codon#{ci} {old_codon}→{new_codon}")
    return RepairLogEntry(
        iteration=iteration,
        kind="applied",
        codon_index=move.codon_swaps[0][0],
        insert_position=move.codon_swaps[0][0] * 3,
        old_codon="+".join(s[1] for s in move.codon_swaps),
        new_codon="+".join(s[2] for s in move.codon_swaps),
        gc_before=round(move.gc_before, 4),
        gc_after=round(move.gc_after, 4),
        fixed_pairs=[(i + 1, j + 1) for (i, j) in move.fixed_pairs],
        remaining_pairs=[(i + 1, j + 1) for (i, j) in move.remaining_pairs],
        remaining_max_tm=round(move.remaining_max_tm, 1) if move.remaining_pairs else None,
        threshold_tm=round(threshold_tm, 1),
        notes="paired swap: " + ", ".join(parts),
    )
