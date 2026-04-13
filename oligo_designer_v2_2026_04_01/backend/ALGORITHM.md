# Oligo Designer Algorithm

How the tool goes from a DNA sequence to a set of oligos you can order from IDT.

## Overview

The input is a DNA sequence (sense strand, 5'->3'). The output is a set of
single-stranded oligos that, when mixed together, will anneal at their overlaps
and be extended by polymerase to produce the full double-stranded DNA product.

The oligos alternate between sense and antisense strands:

```
Sense:       ====1====>          ====3====>          ====5====>
                   <====2====         <====4====
Antisense:
```

Adjacent oligos share an overlap region (~20bp) where they anneal to each other.
The non-overlapping portions (gaps) are filled in by polymerase during assembly.

## Step 1: Build the full construct

If plasmid flanking sequences are provided, they're prepended/appended:

```
[upstream flank] + [insert DNA] + [downstream flank]
```

The flanks are ~20bp of plasmid sequence on either side of the insertion site.
These let the assembled product be inserted into the vector in a subsequent
step (e.g., Gibson assembly or restriction cloning). The flanks are optional —
without them, you just get the linear DNA product.

## Step 2: Generate initial cut points

The algorithm walks along the construct and places "cut points" at regular
intervals. A cut point is where one oligo's unique region ends and the overlap
with the next oligo begins.

```
step = oligo_length - overlap_length
```

With defaults (oligo_length=45, overlap_length=20), step = 25bp. So a cut
point is placed every 25bp along the construct.

Each oligo then spans from one cut point to the next, PLUS the overlap region.
This means adjacent oligos share `overlap_length` bp of identical sequence.

## Step 3: Optimize cut points

This is the key step. Each initial cut point is shifted up to `max_shift`
positions (default: 5bp) in either direction, and the position that produces
the best-scoring overlap is chosen.

### Spacing constraint

Adjacent cut points must remain at least `overlap_length` apart after
optimization. Without this constraint, aggressive shifts could cause
same-strand oligos to overlap each other (e.g., two sense oligos sharing
the same region), which is biologically meaningless. The optimizer skips
any candidate position that would violate this minimum spacing.

### Design principle: only score what cut placement can influence

The scoring function only considers factors where shifting the cut point by
a few bp actually makes a difference. Problems that are equally bad regardless
of where the cut falls (like G-quadruplexes) are flagged as warnings but don't
influence cut placement.

### What the scoring function values (higher = better):

**GC content near 50%**
- Penalty: `abs(GC - 0.5) * 10`
- Why: Overlaps with balanced GC content anneal more reliably. Very AT-rich
  overlaps have low melting temperature (weak annealing). Very GC-rich overlaps
  can form secondary structures.

This is currently the **only** scoring factor. The optimizer picks the
boundary position (within the shift window) that gives the overlap GC%
closest to 50%.

### What the scoring function does NOT consider (and why):

- **Homopolymer placement:** Putting a homopolymer in the overlap avoids
  polymerase slippage in the gap, but forces IDT to synthesize through it
  in **two** oligos instead of one — potentially doubling synthesis cost and
  risk. The tradeoff is unclear, so we don't bias cut placement for this.
  Homopolymers are flagged as warnings instead.

- **G-quadruplexes:** Equally bad in overlaps and gaps — they cause problems
  at IDT synthesis, during annealing, AND during polymerase extension. Moving
  the cut point just moves the problem. Flagged as a warning instead; the real
  fix is codon optimization upstream.

- **Tm (melting temperature):** Slim said not to worry about Tm *matching*
  across overlaps. However, Tm is now computed for each overlap using the
  primer3 nearest-neighbor model (SantaLucia 1998 stacking thermodynamics
  with Owczarzy 2008 Mg²⁺ salt correction) and reported in the output.
  This also enables thermodynamic cross-hybridization checking (see Step 5).

## Step 4: Build oligos from optimized cut points

Once the cut points are finalized:

1. `boundaries = [0] + [cut points] + [end of construct]`
2. Each oligo spans from `boundaries[i]` to `boundaries[i+1] + overlap_length`
3. Odd-numbered oligos (0-indexed) are antisense → reverse-complemented
4. The first and last oligos may be shorter (they don't extend beyond the
   construct ends)

## Step 5: Quality checks on overlaps

### Per-overlap checks

Each overlap is checked for:

| Check | Threshold | Why |
|-------|-----------|-----|
| Homopolymer runs | 4+ consecutive same base | Annealing register slip; IDT synthesis purity |
| G-quadruplex | 4+ G-runs of 3+ in 30bp | Stable secondary structure prevents annealing |

Note: Per-overlap GC% is no longer flagged as an issue — the overlap Tm
(computed via primer3) captures the same concern more precisely. GC% is
still displayed as informational context.

### Thermodynamic analysis (primer3)

Each overlap's **melting temperature (Tm)** is computed using the primer3
nearest-neighbor model with **SantaLucia 1998** stacking thermodynamics and
the **Owczarzy 2008** salt correction for Mg²⁺. This uses the actual stacking
energies of each dinucleotide pair, not a simple %GC formula. The Owczarzy
correction accounts for Mg²⁺ in the assembly buffer, which the SantaLucia
Na⁺-only correction would miss.

### Cross-hybridization check (replaces string-matching uniqueness)

All overlap pairs are checked for off-target binding using primer3's
`calc_heterodimer()`. This does a full dynamic-programming alignment between
two sequences and computes the ΔG and Tm of the best possible duplex they
could form — accounting for mismatches, bulges, dangling ends, and stacking.

An overlap pair is flagged as an **error** if the heterodimer Tm exceeds the
average overlap Tm minus 20°C. This threshold follows Gene2Oligo (Rouillard
et al. 2004, Nucleic Acids Res. 32:W176-W180), the only published tool for
gene synthesis oligo design that specifies a cross-hybridization Tm threshold.
Both orientations are checked (overlap vs overlap, and overlap vs reverse
complement of the other overlap).

This replaces the previous exact-string-match uniqueness check. The
thermodynamic approach is more rigorous: it catches partial complementarity
that string matching would miss, and it ignores similarities that look
concerning on paper but wouldn't actually form stable duplexes at the
reaction temperature.

### Reaction conditions

The thermodynamic calculations require buffer conditions as input:

| Parameter | Default | Description |
|-----------|---------|-------------|
| mv_conc | 50.0 mM | Monovalent cations (Na⁺/K⁺) |
| dv_conc | 2.0 mM | Divalent cations (Mg²⁺) |
| dntp_conc | 0.8 mM | Total dNTPs (0.2 mM each × 4) |
| dna_conc | 250.0 nM | Oligo concentration |
| annealing_temp | 50.0°C | Annealing/reaction temperature |

These defaults are typical PCR/assembly conditions. The exact values depend
on the specific kit and protocol being used (e.g., NEBuilder HiFi has its
own buffer). All are configurable via the API.

Issues are reported as warnings or errors but don't change the tiling.

## Step 6: Whole-sequence complexity scan

Independent of the tiling, the full input sequence is scanned for:

| Feature | Threshold | Severity |
|---------|-----------|----------|
| ~~GC content~~ | ~~Removed~~ | Replaced by the visual GC% track in the assembly viewer. Overlap Tm captures annealing quality more precisely. |
| Homopolymer run | 5+ same base | Warning (5-7) / Error (8+) |
| G-quadruplex | 4+ G-runs of 3+ in 30bp | Error |
| Repeated 10-mers | Same 10bp appears 2+ times | Warning (2x) / Error (3x+) |

These are problems that exist in the input sequence itself, before any oligo
design decisions. Most can only be fixed by codon optimization upstream.

## What the tool CANNOT fix

- **Repetitive sequences** (e.g., GGGGS linkers): No oligo placement can make
  these safe. Must be ordered as pre-made dsDNA (gBlocks) or fixed by using
  different synonymous codons.
- **Homopolymer runs >8bp**: Even in an overlap, long runs cause IDT synthesis
  quality issues. Needs codon optimization.
- **Very low complexity regions**: Need codon optimization to introduce sequence
  diversity.

## Defaults

| Parameter | Default | Range | Notes |
|-----------|---------|-------|-------|
| oligo_length | 45 bp | 30-200 | Standard cheap IDT oligos are up to 45bp. 60bp and 200-mer cost more. |
| overlap_length | 20 bp | 15-40 | ~20bp is standard for overlap assembly. |
| max_shift | 5 bp | — | How far cut points can be shifted to optimize overlap quality. |

---

# Codon Optimization (DNA Chisel)

Upstream of oligo tiling, the tool can optimize a protein sequence into
codon-optimized DNA. This uses DNA Chisel (open-source Python library) as
a constraint solver.

## Pipeline

1. User provides a protein sequence and selects a codon table
2. DNA Chisel reverse-translates, then optimizes codon choices subject to
   constraints (homopolymers, GC content) and objectives (codon frequency)
3. The optimized DNA is displayed in an interactive codon track for manual
   review and editing
4. User can swap individual codons; each swap triggers a live constraint
   re-check
5. When satisfied, user sends the DNA to the oligo tiling tab

## Codon tables

The tool supports:

- **Built-in species tables** from `python_codon_tables` (e.g., H. sapiens,
  E. coli)
- **Custom tables** in Kazusa/GCG format, pasted or uploaded
- **Lab tables** bundled with the tool (e.g., Brian's HumColi hybrid table,
  which picks codons that work well in both human and E. coli expression
  systems)

## Constraints (enforced during optimization)

These are hard constraints — DNA Chisel will not produce a sequence that
violates them.

| Constraint | Default | Why |
|------------|---------|-----|
| EnforceTranslation | Always on | The DNA must encode the input protein |
| AvoidPattern (homopolymers) | 5+ consecutive same base | Homopolymer runs cause polymerase slippage during synthesis and assembly. Runs of 4 are borderline (only 4x G is genuinely concerning at length 4, covered by G-quad check). Default threshold is 5. |
| EnforceGCContent | 25-75% in 50bp windows | Extreme GC causes problems: too high → secondary structures stall polymerase; too low → weak base pairing, poor annealing |

## K-mer uniqueness: warning only, NOT a constraint

**Decision:** Repeated 10-mers (UniquifyAllKmers) are detected and shown as
warnings but are NOT enforced during optimization.

**Why:** Enforcing k-mer uniqueness forces DNA Chisel to swap away from the
best codon to break a repeat. This directly conflicts with the primary goal
of choosing optimal codons for expression. The tradeoff is not worth it
because:

1. Repeated k-mers are only a real problem if both copies land in overlap
   regions during oligo tiling, where they could cause mis-annealing
2. The oligo designer already guards against this — it checks overlap
   uniqueness and runs thermodynamic cross-hybridization analysis
3. Sacrificing codon quality for preemptive k-mer uniqueness is over-cautious

K-mer violations are shown as purple markers on the codon track so the user
can make an informed judgment call.

## Objectives (optimized, not enforced)

| Objective | Description |
|-----------|-------------|
| CodonOptimize | Maximize usage of high-frequency codons from the selected table |

## Two-stage optimization: DNA Chisel + beam search refinement

The pipeline runs two optimizers in sequence and keeps whichever produces lower
total loss.

**Stage 1 — DNA Chisel.** Reverse-translates the protein and resolves
constraints via local backtracking search. It treats codon frequency as a
soft objective. Its weakness is that it processes constraint breaches
left-to-right and commits to early codon choices, sometimes hammering a
single codon with a large frequency sacrifice when it could have spread a
smaller cost across several positions.

**Stage 2 — Beam search refinement** (`gc_refinement.py`). Re-optimizes
*from scratch*, considering every synonymous codon at every position. At each
codon position, candidate paths are extended one codon at a time and pruned
against the same hard constraints (GC window, homopolymers, avoid patterns,
min codon frequency). Paths are bucketed by `(window_gc_count,
last_codon_gc_count)` to maintain spatial diversity, and the top `BEAM_K`
paths per bucket survive.

**Why beam, not exact DP.** Exact DP would require a state of `(position,
last ~50 bases)` to check sliding-window GC. That's ~2^47 distinct states per
position — infeasible. Beam search is the practical approximation: the
buckets are a lossy form of state equivalence, and `BEAM_K` controls how many
representatives per equivalence class survive.

**Why both stages, not just beam.** DNA Chisel acts as a feasibility safety
net: if the beam prunes itself into a corner, we still have DNA Chisel's
output to fall back on. The final output picks whichever has lower loss; the
optimization report reports `source: "beam"` or `source:
"dna_chisel_fallback"` so this is observable.

**Empirical loss/runtime tradeoff** (test_beam_k_sweep.py, 230 aa protein,
GC 40-60%):

| K | time | beam loss | chisel loss |
|---|------|-----------|-------------|
| 5 | 0.6s | 1900 | 5020 |
| 15 | 0.4s | 1840 | 5420 |
| 30 | 0.7s | 1790 | 6020 |
| 60 | 1.5s | 1790 | 5380 |
| 120 | 3.0s | 1800 | 5760 |
| 250 | 6.3s | 1800 | 5660 |

Beam consistently beats DNA Chisel by ~3000 loss points. Quality saturates
around K=30; raising K further does not help on this test case. Default
`BEAM_K=15` is exposed as a user parameter (`beam_k` in the API).

## Progress streaming

`POST /api/codon-optimize-stream` returns Server-Sent Events with progress
ticks (`{type: "progress", stage: "dna_chisel"|"beam", current, total}`)
and a final `{type: "result", data: ...}` event. The frontend consumes this
to show a progress bar during long beam runs.

## Codon coloring in the track viewer

Each codon is colored relative to its synonymous alternatives for the same
amino acid:

| Color | Meaning |
|-------|---------|
| Green | Most frequent codon for this amino acid |
| Yellow | Intermediate frequency |
| Red | Least frequent codon for this amino acid |

**Why relative coloring?** There is no defensible universal threshold for
what makes a codon "rare" in absolute terms. The literature does not provide
a numeric cutoff — Kane 1995 identifies specific rare codons but gives no
threshold; GenScript uses relative adaptiveness < 0.3 (from the CAI framework,
Sharp & Li 1987) but that is a binary rare/not-rare classification.

Our three-tier coloring is fully defensible without arbitrary percentile
cutoffs: best = green, worst = red, everything else = yellow. For amino
acids with only 2 codons, there is no yellow — just green and red. For
amino acids with 1 codon (Met, Trp), it is always green.

**Reference for relative adaptiveness:** The ratio w = codon_freq / max_freq
for each amino acid comes from the Codon Adaptation Index (Sharp & Li 1987,
Nucleic Acids Research 15:1281-1295). GenScript's threshold of w < 0.3 for
"rare" codons is documented in their GenRCA tool:
https://www.genscript.com/gsfiles/tools/Index_Definition_of_GenRCA_Rare_Codon_Analysis_Tool.pdf

## Warning track markers

The codon track shows warnings as colored bars above the codons:

| Warning | Color | Source |
|---------|-------|--------|
| Homopolymer run (4+) | Red | `_scan_warnings()` — runs of 4 are informational, 5+ violate constraints |
| GC content out of range | Teal | DNA Chisel constraint evaluation — shows actual GC% and whether high or low |
| Repeated k-mer | Purple | DNA Chisel UniquifyAllKmers check — detection only, not enforced |

## Live constraint checking

When the user manually swaps a codon:

1. The edited DNA is sent to `/api/check-constraints`
2. All constraints are re-evaluated (including k-mer uniqueness for detection)
3. Failing constraint locations are extracted from DNA Chisel's evaluation
   objects and converted to warning markers on the track
4. The constraints summary and warnings update live

This ensures the user always sees the current constraint status of their
edited sequence.
