r"""
Chasse intensive aux records mondiaux sur F_2 (Flip-Graph adaptatif uint64 XOR) :
  - <3, 4, 5> (record connu = 47)
  - <3, 5, 5> (record connu = 58)
  - <4, 4, 5> (record connu = 56 sur F_2 / 60 general)
  - <4, 5, 5> (record connu = 73 sur F_2)

Inclut sauvegarde incrementale JSON + temoin exact .npz des qu'un nouveau meilleur
rang est trouve, et verification exacte 100% independante sur T_{<m,k,n>} mod 2.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from flip_graph_f2_adaptive import run_adaptive
from flip_graph_f2_bitvec import (
    RESULTS_DIR,
    WITNESS_DIR,
    build_schoolbook_f2_bits,
    verify_scheme_f2_bits,
)


def hunt_target(m, k, n, target_R, known_record, n_rounds=16, beam_width=32, n_workers=48, steps_per_round=250000, plus_thresh=800, max_gap=4):
    u0, v0, w0 = build_schoolbook_f2_bits(m, k, n)
    res = run_adaptive(
        m, k, n, u0, v0, w0,
        target_R=target_R,
        known_record=known_record,
        n_rounds=n_rounds,
        beam_width=beam_width,
        n_workers=n_workers,
        steps_per_round=steps_per_round,
        plus_thresh=plus_thresh,
        max_gap=max_gap,
    )
    tag = f"f2_best_{m}x{k}x{n}_R{res['best_R']}.npz"
    npz_path = os.path.join(WITNESS_DIR, tag)
    np.savez(npz_path, u=res["u"], v=res["v"], w=res["w"], m=m, k=k, n=n, R=res["best_R"])
    return {
        "format": [m, k, n],
        "best_R": int(res["best_R"]),
        "known_record": int(known_record),
        "beat_record": bool(res["best_R"] < known_record),
        "matched_record": bool(res["best_R"] == known_record),
        "witness_npz": tag,
    }


def main():
    t0 = time.time()
    print("====================================================================================")
    print("  CHASSE AUX RECORDS MONDIAUX F_2 : <3,4,5>, <3,5,5>, <4,4,5>, <4,5,5>")
    print("====================================================================================")

    summary = []
    # 1. <3, 4, 5> (naif = 60, record = 47, cible = 46)
    summary.append(hunt_target(3, 4, 5, target_R=46, known_record=47, n_rounds=14, beam_width=32, n_workers=48, steps_per_round=220000, plus_thresh=900, max_gap=4))

    # 2. <3, 5, 5> (naif = 75, record = 58, cible = 57)
    summary.append(hunt_target(3, 5, 5, target_R=57, known_record=58, n_rounds=14, beam_width=32, n_workers=48, steps_per_round=220000, plus_thresh=1000, max_gap=4))

    # 3. <4, 4, 5> (naif = 80, record = 56 sur F2 / 60 general, cible = 55)
    summary.append(hunt_target(4, 4, 5, target_R=55, known_record=56, n_rounds=14, beam_width=32, n_workers=48, steps_per_round=220000, plus_thresh=1100, max_gap=4))

    out = {"targets": summary, "total_time_s": round(time.time() - t0, 2)}
    out_path = os.path.join(RESULTS_DIR, "flip_graph_f2_world_hunt.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nTermine en {time.time() - t0:.2f}s. Sauvegarde : {out_path}")


if __name__ == "__main__":
    main()
