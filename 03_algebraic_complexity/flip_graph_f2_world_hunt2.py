r"""
Chasse aux records mondiaux sur F_2 (Moteur incremental O(R) + collage de blocs optimaux) :
  1. <3, 5, 5> : collage <3, 4, 5> (R=47) + <3, 1, 5> (15) = 62 (record connu = 58 -> cible 57 !)
  2. <4, 4, 5> : collage <4, 4, 4> (R=47) + <4, 4, 1> (16) = 63 (record connu = 56 -> cible 55 !)
  3. <4, 5, 5> : collage <4, 4, 5> (best) + <4, 1, 5> (20) (record connu = 73 -> cible 72 !)
  4. <4, 4, 6> : collage <4, 4, 4> (R=47) + <4, 4, 2> (R<=29, via <2,4,4>) = <= 76 !
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
    print("  CHASSE INTENSIVE AUX RECORDS MONDIAUX F_2 : <3,5,5>, <4,4,5>, <4,5,5>")
    print("====================================================================================")

    # Charger nos temoins exacts <3,4,5> R=47 et <4,4,4> R=47
    d345 = np.load(os.path.join(WITNESS_DIR, "f2_best_3x4x5_R47.npz"))
    u345, v345, w345 = d345["u"], d345["v"], d345["w"]

    d444 = np.load(os.path.join(WITNESS_DIR, "f2_best_4x4x4_R47.npz"))
    u444, v444, w444 = d444["u"], d444["v"], d444["w"]

    summary = []

    # 1. <3, 5, 5> (record connu = 58, cible = 57)
    # Collage sur le 2e axe : rotation <3,4,5> -> <4,5,3> -> <5,3,4>, coller <5,3,1>(15) -> <5,3,5>, rotation -> <3,5,5> !
    u453, v453, w453, a, b, c = rotate_cycle_scheme(u345, v345, w345, 3, 4, 5)
    u534, v534, w534, a, b, c = rotate_cycle_scheme(u453, v453, w453, a, b, c)
    assert (a, b, c) == (5, 3, 4)
    u531, v531, w531 = build_schoolbook_f2_bits(5, 3, 1)
    u535, v535, w535 = glue_along_n(u534, v534, w534, 5, 3, 4, u531, v531, w531, 1)
    u355_init, v355_init, w355_init, a, b, c = rotate_cycle_scheme(u535, v535, w535, 5, 3, 5)
    assert (a, b, c) == (3, 5, 5)

    r355 = run_fast_inc_hunt(
        3, 5, 5, [(u355_init, v355_init, w355_init)], target_R=57, known_record=58,
        n_rounds=20, beam_width=40, n_workers=48, steps_per_round=700000, plus_thresh=1400, max_gap=4
    )
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_3x5x5_R{r355['best_R']}.npz"), u=r355["u"], v=r355["v"], w=r355["w"], m=3, k=5, n=5, R=r355["best_R"])
    summary.append({"format": [3, 5, 5], "best_R": int(r355["best_R"]), "known_record": 58, "beat_record": bool(r355["best_R"] < 58)})

    # 2. <4, 4, 5> (record connu = 56 sur F2, cible = 55)
    u441, v441, w441 = build_schoolbook_f2_bits(4, 4, 1)
    glued_445 = glue_along_n(u444, v444, w444, 4, 4, 4, u441, v441, w441, 1)
    r445 = run_fast_inc_hunt(
        4, 4, 5, [glued_445], target_R=55, known_record=56,
        n_rounds=20, beam_width=40, n_workers=48, steps_per_round=700000, plus_thresh=1400, max_gap=4
    )
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_4x4x5_R{r445['best_R']}.npz"), u=r445["u"], v=r445["v"], w=r445["w"], m=4, k=4, n=5, R=r445["best_R"])
    summary.append({"format": [4, 4, 5], "best_R": int(r445["best_R"]), "known_record": 56, "beat_record": bool(r445["best_R"] < 56)})

    # 3. <4, 5, 5> depuis collage <4,4,5>(best) + <4,1,5>(20) (record connu = 73, cible = 72)
    # Rotation <4,4,5> -> <4,5,4> -> <5,4,4>, coller <5,4,1>(20) -> <5,4,5>, rotation -> <4,5,5> !
    u454, v454, w454, a, b, c = rotate_cycle_scheme(r445["u"], r445["v"], r445["w"], 4, 4, 5)
    u544, v544, w544, a, b, c = rotate_cycle_scheme(u454, v454, w454, a, b, c)
    assert (a, b, c) == (5, 4, 4)
    u541, v541, w541 = build_schoolbook_f2_bits(5, 4, 1)
    u545, v545, w545 = glue_along_n(u544, v544, w544, 5, 4, 4, u541, v541, w541, 1)
    u455_init, v455_init, w455_init, a, b, c = rotate_cycle_scheme(u545, v545, w545, 5, 4, 5)
    assert (a, b, c) == (4, 5, 5)

    r455 = run_fast_inc_hunt(
        4, 5, 5, [(u455_init, v455_init, w455_init)], target_R=72, known_record=73,
        n_rounds=20, beam_width=40, n_workers=48, steps_per_round=700000, plus_thresh=1600, max_gap=4
    )
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_4x5x5_R{r455['best_R']}.npz"), u=r455["u"], v=r455["v"], w=r455["w"], m=4, k=5, n=5, R=r455["best_R"])
    summary.append({"format": [4, 5, 5], "best_R": int(r455["best_R"]), "known_record": 73, "beat_record": bool(r455["best_R"] < 73)})

    out = {"summary": summary, "total_time_s": round(time.time() - t0, 2)}
    out_path = os.path.join(RESULTS_DIR, "flip_graph_f2_world_hunt2.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
