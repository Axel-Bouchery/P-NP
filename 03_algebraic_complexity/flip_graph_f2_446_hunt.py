r"""
Chasse sur <2, 4, 4> (cible 26) puis <4, 4, 6> (collage <4,4,4> R=47 + <4,4,2> R<=26 = 73)
avec marches longues Arai et al. (plus_thresh = 10 000, steps = 1,5M/worker/ronde).
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from flip_graph_f2_beam import glue_along_n, rotate_cycle_scheme
from flip_graph_f2_bitvec import RESULTS_DIR, WITNESS_DIR, build_schoolbook_f2_bits
from flip_graph_f2_fast_inc import run_fast_inc_hunt


def main():
    t0 = time.time()
    print("====================================================================================")
    print("  CHASSE <2,4,4> -> <4,4,6> & REGLAGE LONGUE MARCHE SUR <4,4,5> / <3,5,5>")
    print("====================================================================================")

    # 1. <2, 4, 4> (naif = 32, record = 26)
    r244 = run_fast_inc_hunt(
        2, 4, 4, [build_schoolbook_f2_bits(2, 4, 4)],
        target_R=26, known_record=26, n_rounds=10, beam_width=32,
        n_workers=48, steps_per_round=400000, plus_thresh=2500, max_gap=3
    )
    # Rotation cyclique <2, 4, 4> -> <4, 4, 2>
    u442, v442, w442, a, b, c = rotate_cycle_scheme(r244["u"], r244["v"], r244["w"], 2, 4, 4)
    assert (a, b, c) == (4, 4, 2)

    # 2. Charger <4, 4, 4> R = 47 et coller avec <4, 4, 2> -> <4, 4, 6>
    d444 = np.load(os.path.join(WITNESS_DIR, "f2_best_4x4x4_R47.npz"))
    glued_446 = glue_along_n(d444["u"], d444["v"], d444["w"], 4, 4, 4, u442, v442, w442, 2)

    r446 = run_fast_inc_hunt(
        4, 4, 6, [glued_446],
        target_R=66, known_record=73, n_rounds=24, beam_width=48,
        n_workers=48, steps_per_round=1200000, plus_thresh=8000, max_gap=4
    )
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_4x4x6_R{r446['best_R']}.npz"), u=r446["u"], v=r446["v"], w=r446["w"], m=4, k=4, n=6, R=r446["best_R"])

    # 3. Re-tenter <3, 5, 5> depuis notre meilleur R=60 avec marches longues (plus_thresh = 8 000)
    d355 = np.load(os.path.join(WITNESS_DIR, "f2_best_3x5x5_R60.npz"))
    r355 = run_fast_inc_hunt(
        3, 5, 5, [(d355["u"], d355["v"], d355["w"])],
        target_R=57, known_record=58, n_rounds=20, beam_width=48,
        n_workers=48, steps_per_round=1200000, plus_thresh=8000, max_gap=4
    )
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_3x5x5_R{r355['best_R']}.npz"), u=r355["u"], v=r355["v"], w=r355["w"], m=3, k=5, n=5, R=r355["best_R"])

    out = {
        "2_4_4": int(r244["best_R"]),
        "4_4_6": int(r446["best_R"]),
        "3_5_5": int(r355["best_R"]),
        "total_time_s": round(time.time() - t0, 2),
    }
    out_path = os.path.join(RESULTS_DIR, "flip_graph_f2_446_hunt.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
