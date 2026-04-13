"""Empirical sweep over BEAM_K values to characterize loss/runtime tradeoff.

Runs the full pipeline (DNA Chisel + beam refinement) on a tight-constraint
test sequence at several K values and prints loss + wall time for each.
"""
import time
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from api import _run_codon_optimize, CodonOptRequest

# A protein with several tight regions: AT-rich stretches, GC-rich stretches,
# repeats. The example sequence from the UI ("Use example sequence" button).
PROTEIN = (
    "MKGLVTDDEGQPIPGATIKVNRNKKPVTSSARGEYWRLLLPGNYTLTASAQEYISLTISI"
    "HIPEKTDPDDWSKPAQRQDFRLVRKSTGVTKPLTNSNNPTIITNEKNVPTELSFPIEPS"
    "ASSINSPPFITTKSNLKPNFFPTELTAQTEEPRHFENKGKQNNKPFQNNSNKDIKNNEN"
    "PRIEATTRRSNQSRGRGGRRCKCANGKRRCKKCSKRSKTRNRNLRNAQKKQK"
)


TIGHT = dict(gc_min=0.40, gc_max=0.60, homo_gc=4, homo_at=5)
LOOSE = dict(gc_min=0.30, gc_max=0.70, homo_gc=4, homo_at=6)
PROFILE = TIGHT


def run_one(beam_k: int) -> dict:
    req = CodonOptRequest(
        protein_sequence=PROTEIN,
        species="h_sapiens_9606",
        avoid_homopolymers_gc=PROFILE["homo_gc"],
        avoid_homopolymers_at=PROFILE["homo_at"],
        gc_min=PROFILE["gc_min"],
        gc_max=PROFILE["gc_max"],
        gc_window=50,
        uniquify_kmers=10,
        avoid_patterns=[],
        min_codon_frequency=10.0,
        beam_k=beam_k,
    )
    t0 = time.perf_counter()
    response = _run_codon_optimize(req)
    elapsed = time.perf_counter() - t0
    rep = response.optimization_report
    return {
        "beam_k": beam_k,
        "elapsed_s": elapsed,
        "total_loss": rep.get("total_loss"),
        "chosen_score": rep.get("chosen_score"),
        "ideal_score": rep.get("ideal_score"),
        "codons_changed": rep.get("codons_changed"),
        "source": rep.get("source"),
        "beam_loss": rep.get("beam_loss"),
        "chisel_loss": rep.get("chisel_loss"),
    }


def main() -> None:
    print(f"Protein length: {len(PROTEIN)} aa")
    print(f"DNA length: {len(PROTEIN) * 3} bp\n")
    print(f"{'K':>5} {'time(s)':>10} {'loss':>10} {'beam_loss':>12} {'chisel_loss':>14} {'changed':>10} {'source':>22}")
    print("-" * 90)
    for k in [5, 15, 30, 60, 120, 250]:
        try:
            r = run_one(k)
            print(
                f"{r['beam_k']:>5} "
                f"{r['elapsed_s']:>10.2f} "
                f"{r['total_loss']:>10.1f} "
                f"{(r['beam_loss'] if r['beam_loss'] is not None else float('nan')):>12.1f} "
                f"{(r['chisel_loss'] if r['chisel_loss'] is not None else float('nan')):>14.1f} "
                f"{r['codons_changed']:>10} "
                f"{str(r['source']):>22}"
            )
        except Exception as e:
            print(f"{k:>5}  ERROR: {e}")


if __name__ == "__main__":
    main()
