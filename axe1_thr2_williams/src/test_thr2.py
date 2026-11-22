from pathlib import Path
import sys, time
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from thr2_lab import make_instances, solve_batch, NF

ok = True
for sat in (False, True):
    for n in (8, 10, 12):
        Ws, ts, vs, Ts = make_instances(n, n, 20, seed=n + 100 * sat, sat=sat)
        for th in (np.zeros(NF), np.eye(NF)[1], np.random.default_rng(0).normal(size=NF)):
            st, nodes = solve_batch(Ws, ts, vs, Ts, th, 1 << 40)
            expect = 1 if sat else 0
            good = bool((st == expect).all())
            ok &= good
            print(f"sat={sat} n={n} statuts ok={good} noeuds moyens={nodes.mean():.0f} (2^(n+1)={2**(n+1)})")
t0 = time.time()
Ws, ts, vs, Ts = make_instances(18, 18, 12, seed=1)
t1 = time.time()
st, nodes = solve_batch(Ws, ts, vs, Ts, np.zeros(NF), 1 << 40)
print(f"n=18 : generation {t1-t0:.1f}s, resolution {time.time()-t1:.1f}s, noeuds moyens {nodes.mean():.0f}")
print("TOUT OK" if ok else "ECHEC")
