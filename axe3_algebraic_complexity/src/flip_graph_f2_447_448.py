r"""
Chasse sur formats rectangulaires moyens F_2 (64-bit) avec collage optimal + grand bassin :
  1. <3, 4, 7> : collage <3, 4, 4> (38) + <3, 4, 3> (29) = 67 !
  2. <4, 4, 7> : collage <4, 4, 4> (47) + <4, 4, 3> (38) = 85 !
  3. <4, 4, 8> : collage <4, 4, 4> (47) + <4, 4, 4> (47) = 94 !
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
from flip_graph_f2_nondestructive import run_nondest


def main():
    t0 = time.time()
    print("====================================================================================")
    print("  CHASSE <3,4,7>, <4,4,7> ET <4,4,8> SUR F_2 (COLLAGE OPTIMAL + GAUSS COPIE)")
    print("====================================================================================")

    # Recuperer <3, 3, 4> R = 29 -> rotation -> <3, 4, 3> R = 29
    r334 = run_fast_inc_hunt(3, 3, 4, [build_schoolbook_f2_bits(3, 3, 4)], target_R=29, known_record=29, n_rounds=8, beam_width=24, n_workers=32, steps_per_round=300000, plus_thresh=900, max_gap=4)
    u343, v343, w343, a, b, c = rotate_cycle_scheme(r334["u"], r334["v"], r334["w"], 3, 3, 4)
    assert (a, b, c) == (3, 4, 3)

    # Recuperer <3, 4, 4> R = 38 et <4, 4, 4> R = 47
    r344 = run_fast_inc_hunt(3, 4, 4, [build_schoolbook_f2_bits(3, 4, 4)], target_R=38, known_record=38, n_rounds=8, beam_width=24, n_workers=32, steps_per_round=400000, plus_thresh=900, max_gap=4)
    u443, v443, w443, a, b, c = rotate_cycle_scheme(r344["u"], r344["v"], r344["w"], 3, 4, 4)
    assert (a, b, c) == (4, 4, 3)

    d444 = np.load(os.path.join(WITNESS_DIR, "f2_best_4x4x4_R47.npz"))
    u444, v444, w444 = d444["u"], d444["v"], d444["w"]

    # 1. <3, 4, 7> : collage <3,4,4>(38) + <3,4,3>(29) = 67 (et <3,4,6>(56) + <3,4,1>(12) = 68)
    d346 = np.load(os.path.join(WITNESS_DIR, "f2_best_3x4x6_R56.npz"))
    u341, v341, w341 = build_schoolbook_f2_bits(3, 4, 1)
    g347_a = glue_along_n(r344["u"], r344["v"], r344["w"], 3, 4, 4, u343, v343, w343, 3)
    g347_b = glue_along_n(d346["u"], d346["v"], d346["w"], 3, 4, 6, u341, v341, w341, 1)
    r347 = run_nondest(3, 4, 7, [g347_a, g347_b], target_R=62, known_record=66, n_rounds=16, beam_width=48, n_workers=48, steps_per_round=800000, plus_thresh=6000, max_gap=5)
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_3x4x7_R{r347['best_R']}.npz"), u=r347["u"], v=r347["v"], w=r347["w"], m=3, k=4, n=7, R=r347["best_R"])

    # 2. <4, 4, 7> : collage <4,4,4>(47) + <4,4,3>(38) = 85
    g447 = glue_along_n(u444, v444, w444, 4, 4, 4, u443, v443, w443, 3)
    r447 = run_nondest(4, 4, 7, [g447], target_R=78, known_record=85, n_rounds=16, beam_width=48, n_workers=48, steps_per_round=800000, plus_thresh=6000, max_gap=5)
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_4x4x7_R{r447['best_R']}.npz"), u=r447["u"], v=r447["v"], w=r447["w"], m=4, k=4, n=7, R=r447["best_R"])

    # 3. <4, 4, 8> : collage <4,4,4>(47) + <4,4,4>(47) = 94
    g448 = glue_along_n(u444, v444, w444, 4, 4, 4, u444, v444, w444, 4)
    r448 = run_nondest(4, 4, 8, [g448], target_R=88, known_record=94, n_rounds=16, beam_width=48, n_workers=48, steps_per_round=800000, plus_thresh=6000, max_gap=5)
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_4x4x8_R{r448['best_R']}.npz"), u=r448["u"], v=r448["v"], w=r448["w"], m=4, k=4, n=8, R=r448["best_R"])

    out = {"3_4_7": int(r347["best_R"]), "4_4_7": int(r447["best_R"]), "4_4_8": int(r448["best_R"]), "total_time_s": round(time.time() - t0, 2)}
    out_path = os.path.join(RESULTS_DIR, "flip_graph_f2_447_448.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
