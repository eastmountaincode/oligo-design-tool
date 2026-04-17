"""Benchmark: does extending the beam's bucket key to include homopolymer
state (run-base + run-length) meaningfully improve search quality?

For a battery of test proteins + constraint profiles, run the full codon-opt
pipeline twice — once with include_homo_in_bucket_key=True, once with False —
and report for each: beam source (did it succeed?), total loss, wall time,
whether gBlock fallback kicked in.

Run from the backend folder:
    python test_bucket_key_homo.py
"""
from __future__ import annotations
import time
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from api import _run_codon_optimize, CodonOptRequest


# -----------------------------------------------------------------------------
# Test proteins
# -----------------------------------------------------------------------------

# Small, well-behaved protein — expect both to succeed trivially.
PROTEIN_SMALL = (
    "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPT"
    "LVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTL"
    "VNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLA"
    "DHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
)

# Medium, irregular — example 1 from the UI. Has some challenging regions.
PROTEIN_MEDIUM = (
    "MKGLVTDDEGQPIPGATIKVNRNKKPVTSSARGEYWRLLLPGNYTLTASAQEYISLTISI"
    "HIPEKTDPDDWSKPAQRQDFRLVRKSTGVTKPLTNSNNPTIITNEKNVPTELSFPIEPS"
    "ASSINSPPFITTKSNLKPNFFPTELTAQTEEPRHFENKGKQNNKPFQNNSNKDIKNNEN"
    "PRIEATTRRSNQSRGRGGRRCKCANGKRRCKKCSKRSKTRNRNLRNAQKKQK"
)

# Hard, repeat-heavy — example 2 from the UI (Fc fusion with multiple GGGGS
# linkers). The linker regions force tight codon corners.
PROTEIN_LINKERED = (
    "MPMGSLQPLATLYLLGMLVASVLAGGGGSGGGGSRVRRHGEGTFTSDLSKQMEEEAVRLFI"
    "EWLKNGGPSSGAPPPSGGGGGSGGGGSGGGGSGGGGSTGEPKSCDKTHTCPPCPAPELLGG"
    "PSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNS"
    "TYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELT"
    "KNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQG"
    "NVFSCSVMHEALHNHYTQKSLSLSPGKGGGGSGGGGSENLYFQGGGGGSGGGGSGDPSQVQ"
    "LQQSGPELMKPGASVKISCKASGNVFSSSWMNWVKQRPGKGLEWIGRIYPGDGETNYNGKF"
    "KGKATLTADKSSSTAYMQLSSLTSEDSAVYFCARSDSNYDWYFDVWGTGTTVTVSSAGGGG"
    "SGGGGSGGGGSGGGGSQIVLTQSPAIMSASPGEKVTMTCSASSSVSYMYWYQQKPGSSPRL"
    "LIYDTSNLASGVPVRFSGSGSGTSYSLTISRMEAEDAATYYCQQWSSYPPTFGAGTKLELK"
    "RADAAPTVSIFPPSSEQLTSGGASVGS"
)

# Synthetic stress test — poly-glycine + poly-alanine regions force narrow
# codon choices with homopolymer-prone tails.
PROTEIN_STRESS = (
    "MGAAAAAAAAAPPPPPPPPKKKKKKKKLLLLLLLLMMMMMMMM"
    "GGGGGGGGGAAAAAAAAPPPPPPPPKKKKKKKKLLLLLLLLMMMMMMMM"
    "GGGGGGGGGAAAAAAAAPPPPPPPPKKKKKKKKLLLLLLLLMMMMMMMM"
)

TEST_PROTEINS: list[tuple[str, str]] = [
    ("small / GFP-like", PROTEIN_SMALL),
    ("medium / example1", PROTEIN_MEDIUM),
    ("hard / linker-heavy Fc fusion", PROTEIN_LINKERED),
    ("stress / synthetic homo-rich", PROTEIN_STRESS),
]


# -----------------------------------------------------------------------------
# Constraint profiles
# -----------------------------------------------------------------------------

PROFILES: list[tuple[str, dict]] = [
    ("loose",  dict(gc_min=0.30, gc_max=0.70, homo_gc=4, homo_at=6)),
    ("tight",  dict(gc_min=0.40, gc_max=0.60, homo_gc=4, homo_at=5)),
    ("tightest", dict(gc_min=0.45, gc_max=0.55, homo_gc=3, homo_at=4)),
]


# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------

def run_one(protein: str, profile: dict, *, include_homo: bool, beam_k: int, auto_detect: bool = False) -> dict:
    req = CodonOptRequest(
        protein_sequence=protein,
        species="h_sapiens_9606",
        avoid_homopolymers_gc=profile["homo_gc"],
        avoid_homopolymers_at=profile["homo_at"],
        gc_min=profile["gc_min"],
        gc_max=profile["gc_max"],
        gc_window=50,
        uniquify_kmers=10,
        avoid_patterns=[],
        min_codon_frequency=10.0,
        beam_k=beam_k,
        include_homo_in_bucket_key=include_homo,
        auto_detect_gblocks=auto_detect,
    )
    t0 = time.perf_counter()
    try:
        response = _run_codon_optimize(req)
        elapsed = time.perf_counter() - t0
        rep = response.optimization_report or {}
        return {
            "ok": True,
            "elapsed_s": elapsed,
            "source": rep.get("source", "?"),
            "total_loss": rep.get("total_loss"),
            "chosen_score": rep.get("chosen_score"),
            "ideal_score": rep.get("ideal_score"),
            "gblocks_auto": bool(response.gblocks_auto_applied),
            "n_gblocks": len(response.suggested_gblocks or []),
            "gblock_bp": sum((g.end_bp - g.start_bp) for g in (response.suggested_gblocks or [])),
        }
    except Exception as e:
        elapsed = time.perf_counter() - t0
        return {
            "ok": False,
            "elapsed_s": elapsed,
            "error": str(e)[:80],
        }


def fmt_row(label: str, r: dict) -> str:
    if not r["ok"]:
        return f"  {label:<8}  FAIL  {r['elapsed_s']:>6.2f}s  {r.get('error', '')}"
    source = r["source"]
    loss = r["total_loss"]
    auto = "auto-gBlock" if r["gblocks_auto"] else ""
    gb = f"{r['n_gblocks']}gb/{r['gblock_bp']}bp" if r["n_gblocks"] else ""
    return (
        f"  {label:<8}  {source:<18}  "
        f"loss={loss!s:>6}  "
        f"{r['elapsed_s']:>6.2f}s  "
        f"{auto:<12} {gb}"
    )


def main() -> None:
    BEAM_K = 60  # a realistic mid-range value
    print(f"Beam K = {BEAM_K}")
    print(f"Comparing include_homo_in_bucket_key False vs True (auto_detect_gblocks=False: isolate beam behavior)\n")

    for pname, protein in TEST_PROTEINS:
        print(f"=== {pname} ({len(protein)} aa) ===")
        for profname, profile in PROFILES:
            print(f" Profile: {profname}  (GC {profile['gc_min']:.0%}-{profile['gc_max']:.0%}, homo_gc={profile['homo_gc']}, homo_at={profile['homo_at']})")
            r_off = run_one(protein, profile, include_homo=False, beam_k=BEAM_K, auto_detect=False)
            r_on  = run_one(protein, profile, include_homo=True,  beam_k=BEAM_K, auto_detect=False)
            print(fmt_row("OFF", r_off))
            print(fmt_row("ON",  r_on))

            if r_off["ok"] and r_on["ok"]:
                lo = r_off["total_loss"] if r_off["total_loss"] is not None else 0
                ln = r_on["total_loss"] if r_on["total_loss"] is not None else 0
                delta_loss = lo - ln
                delta_time = r_on["elapsed_s"] - r_off["elapsed_s"]
                flip = ""
                if r_off["source"] != r_on["source"]:
                    flip = f"   source flip: {r_off['source']} -> {r_on['source']}"
                print(f"   delta_loss (OFF-ON) = {delta_loss:+.1f}  "
                      f"delta_time (ON-OFF) = {delta_time:+.2f}s{flip}")
        print()


if __name__ == "__main__":
    main()
