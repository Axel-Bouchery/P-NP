"""
Extension exacte des polynomes quadratiques P_3 sur {0,1}^m sans zero de poids <= 3
pour m = 5, 6, 7, 8, 9, 10 (avec canonisation ou filtre direct) !
"""

import time


def extend_no_wt3():
    # Etape m = 5 : trouver tous les (c, adj) sur m = 5 sans zero de poids 1, 2, 3
    # Fixer c_0 = 1 par multiplication par +-1 mod 3
    m = 5
    pairs5 = [(i, j) for i in range(m) for j in range(i + 1, m)]
    triples5 = [(i, j, k) for i in range(m) for j in range(i + 1, m) for k in range(j + 1, m)]

    bad = []
    for c_code in range(1 << 4):  # c_0 = 1 fixe
        c = (1,) + tuple(1 + ((c_code >> i) & 1) for i in range(4))
        allowed_d = [[v for v in range(3) if v != (-c[i] - c[j]) % 3] for (i, j) in pairs5]
        for mask in range(1 << 10):
            adj = [[0] * 5 for _ in range(5)]
            for idx, (i, j) in enumerate(pairs5):
                v = allowed_d[idx][(mask >> idx) & 1]
                adj[i][j] = adj[j][i] = v
            ok = True
            for (i, j, k) in triples5:
                if (c[i] + c[j] + c[k] + adj[i][j] + adj[i][k] + adj[j][k]) % 3 == 0:
                    ok = False
                    break
            if ok:
                bad.append((c, tuple(tuple(row) for row in adj)))

    print(f"[m=5] (c_0=1) total sans zero de poids <= 3 : {len(bad)}")

    # Etendre de m -> m+1 (on ne garde qu'un representant ou on compte)
    # Pour rester rapide, on dedoublonne par hash de degrès / canonisation simple ou on teste d'abord la taille
    cur = bad
    for m_next in range(6, 12):
        m_prev = m_next - 1
        nxt = []
        # Pre-calculer paires de m_prev
        prev_pairs = [(i, j) for i in range(m_prev) for j in range(i + 1, m_prev)]
        for (c, adj) in cur:
            for c_new in (1, 2):
                allowed = [[v for v in range(3) if v != (-c[i] - c_new) % 3] for i in range(m_prev)]
                for emask in range(1 << m_prev):
                    col = [allowed[i][(emask >> i) & 1] for i in range(m_prev)]
                    ok = True
                    for (i, j) in prev_pairs:
                        if (c[i] + c[j] + c_new + adj[i][j] + col[i] + col[j]) % 3 == 0:
                            ok = False
                            break
                    if ok:
                        new_adj = tuple(adj[i] + (col[i],) for i in range(m_prev)) + (tuple(col) + (0,),)
                        nxt.append((c + (c_new,), new_adj))
                        if len(nxt) >= 50000:
                            break
                if len(nxt) >= 50000:
                    break
            if len(nxt) >= 50000:
                break
        print(f"[m={m_next}] sans zero de poids <= 3 : {len(nxt)}" + (" (cap 50000)" if len(nxt) >= 50000 else " (EXACT)"))
        if len(nxt) == 0:
            break
        cur = nxt

    # Si cur non vide, afficher un exemple au dernier m atteint
    if cur:
        c, adj = cur[0]
        print(f"Exemple a m={len(c)} : c={c}")
        for row in adj:
            print("  ", row)


if __name__ == "__main__":
    t0 = time.time()
    extend_no_wt3()
    print(f"Termine en {time.time() - t0:.2f}s")
