r"""
Flip-Graph F_2 par Recuit / Marche libre a grand bassin (sans reset premature) :
  - Chaque worker marche librement dans le graphe de flips jusqu'a `max_gap = 6`
    au-dessus de son meilleur rang, avec transitions PLUS probabilistes quand
    `since_reduce` depasse un seuil aleatoire (2000..12000 flips).
  - Detection Gauss de tranche non-destructive sur copie tous les 64 flips.
  - Test sur <3, 4, 6> (record 56), <4, 4, 6> (record 73), <3, 5, 5> (record 58).
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
from flip_graph_f2_nondestructive import run_nondest


def main():
    t0 = time.time()
    print("====================================================================================")
    print("  CHASSE GRAND BASSIN (max_gap=6) SUR <3,4,6>, <4,4,6> ET <3,5,5>")
    print("====================================================================================")

    # 1. Construire <3, 4, 6> depuis <3, 4, 5> (47) + <3, 4, 1> (12) = 59 ET <3,4,4>(38)+<3,4,2>(20)=58
    d345 = np.load(os.path.join(WITNESS_DIR, "f2_best_3x4x5_R47.npz"))
    u341, v341, w341 = build_schoolbook_f2_bits(3, 4, 1)
    g346_a = glue_along_n(d345["u"], d345["v"], d345["w"], 3, 4, 5, u341, v341, w341, 1)

    r346 = run_nondest(
        3, 4, 6, [g346_a, build_schoolbook_f2_bits(3, 4, 6)],
        target_R=55, known_record=56, n_rounds=16, beam_width=48,
        n_workers=48, steps_per_round=900000, plus_thresh=5000, max_gap=6
    )
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_3x4x6_R{r346['best_R']}.npz"), u=r346["u"], v=r346["v"], w=r346["w"], m=3, k=4, n=6, R=r346["best_R"])

    # 2. <4, 4, 6> avec max_gap = 6 depuis naif (96) + collage (73)
    d446 = np.load(os.path.join(WITNESS_DIR, "f2_best_4x4x6_R73.npz"))
    r446 = run_nondest(
        4, 4, 6, [(d446["u"], d446["v"], d446["w"])],
        target_R=70, known_record=73, n_rounds=16, beam_width=48,
        n_workers=48, steps_per_round=1000000, plus_thresh=6000, max_gap=6
    )
    np.savez(os.path.join(WITNESS_DIR, f"f2_best_4x4x6_R{r446['best_R']}.npz"), u=r446["u"], v=r446["v"], w=r446["w"], m=4, k=4, n=6, R=r446["best_R"])

    out = {
        "3_4_6": int(r346["best_R"]),
        "4_4_6": int(r446["best_R"]),
        "total_time_s": round(time.time() - t0, 2),
    }
    out_path = os.path.join(RESULTS_DIR, "flip_graph_f2_sa_hunt.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
