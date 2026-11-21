r"""
Fermeture exhaustive EXACTE de N_gen^{(6)}((1, 2), OR) (cas P_2(x) = |x|_H mod 2,
i.e. polynomes quadratiques P_3 sur {0, 1}^n sans AUCUN zero non nul de poids PAIR) !
Propriete strictement hereditaire par suppression de sommet :
On enumere les classes canoniques mod S_k x {+-1} de k = 2 jusqu'a extinction !
"""

import time
from itertools import combinations


def canonical_key_p3(n, c, adj):
    """
    Canonisation 1-WL + quotient par vrais jumeaux (exacte, tres rapide)
    modulo S_n x {+-1} sur (c, adj) mod 3.
    """
    def canon_sign(cs, a):
        # 1-WL colors
        colors = [cs[i] for i in range(n)]
        for _ in range(n):
            sigs = []
            for i in range(n):
                nbrs = sorted((a[i][j], colors[j]) for j in range(n) if j != i)
                sigs.append((colors[i], tuple(nbrs)))
            uniq = {s: idx for idx, s in enumerate(sorted(set(sigs)))}
            new_colors = [uniq[s] for s in sigs]
            if new_colors == colors:
                break
            colors = new_colors
        # Grouper par couleur ; trier les classes et faire un backtracking exact sur les classes
        # de taille > 1 (en quotientant par les vrais jumeaux qui ont memes c_i et memes lignes adj)
        groups = {}
        for i, col in enumerate(colors):
            groups.setdefault(col, []).append(i)
        ordered_groups = [groups[col] for col in sorted(groups.keys())]

        # Verifier si chaque groupe est forme de vrais jumeaux (ou taille <= 1)
        # Si un groupe a des non-jumeaux, backtracking exact sur les permutations du groupe
        best = None

        def build_perm(g_idx, cur_perm):
            nonlocal best
            if g_idx == len(ordered_groups):
                c_p = tuple(cs[u] for u in cur_perm)
                e_p = tuple(a[cur_perm[i]][cur_perm[j]] for i in range(n) for j in range(i + 1, n))
                cand = (c_p, e_p)
                if best is None or cand < best:
                    best = cand
                return
            grp = ordered_groups[g_idx]
            if len(grp) == 1:
                build_perm(g_idx + 1, cur_perm + grp)
                return
            # Verifier si tous les sommets de grp sont des vrais jumeaux entre eux
            # (meme edge intra-groupe entre toutes paires, et memes edges vers l'exterieur)
            ref = grp[0]
            intra = a[grp[0]][grp[1]]
            is_twin_clique = True
            for x_i in range(len(grp)):
                u = grp[x_i]
                for x_j in range(x_i + 1, len(grp)):
                    if a[u][grp[x_j]] != intra:
                        is_twin_clique = False
                        break
                if not is_twin_clique:
                    break
                for w in range(n):
                    if w not in grp and a[u][w] != a[ref][w]:
                        is_twin_clique = False
                        break
                if not is_twin_clique:
                    break
            if is_twin_clique:
                build_perm(g_idx + 1, cur_perm + grp)
                return

            # Sinon, choisir chaque element du groupe et raffiner (simple backtracking)
            from itertools import permutations
            for perm in permutations(grp):
                build_perm(g_idx + 1, cur_perm + list(perm))

        build_perm(0, [])
        return best

    k1 = canon_sign(c, adj)
    c_neg = tuple((-v) % 3 for v in c)
    a_neg = tuple(tuple((-v) % 3 for v in row) for row in adj)
    k2 = canon_sign(c_neg, a_neg)
    return k1 if k1 < k2 else k2


def solve_exact_mixed_1_2():
    t0 = time.time()
    # Niveau k = 2 : c in {0,1,2}^2, d_{01} in {0,1,2} tel que c_0 + c_1 + d_{01} != 0 mod 3
    cur = {}
    for c0 in range(3):
        for c1 in range(3):
            for d01 in range(3):
                if (c0 + c1 + d01) % 3 == 0:
                    continue
                c = (c0, c1)
                adj = ((0, d01), (d01, 0))
                key = canonical_key_p3(2, c, adj)
                if key not in cur:
                    cur[key] = (c, adj)

    print(f"[k = 2] classes canoniques (sans zero de poids pair) : {len(cur)}")

    for k_next in range(3, 10):
        t1 = time.time()
        k_prev = k_next - 1
        # Pre-calculer les sous-ensembles de taille IMPAIRE de {0..k_prev-1}
        # (en ajoutant le nouveau sommet k_prev, ils donnent TOUS les nouveaux vecteurs de poids PAIR !)
        odd_masks = []
        for mask in range(1, 1 << k_prev):
            if (mask.bit_count() & 1) == 1:
                odd_masks.append(mask)

        nxt = {}
        tested = 0
        valid = 0
        for (c, adj) in cur.values():
            # Pre-calculer P_3(mask) pour tous les masques impairs de k_prev
            p3_odd = {}
            bits_of = {}
            for mask in odd_masks:
                bits = [i for i in range(k_prev) if (mask >> i) & 1]
                bits_of[mask] = bits
                val = sum(c[i] for i in bits)
                for a_i in range(len(bits)):
                    for b_i in range(a_i + 1, len(bits)):
                        val += adj[bits[a_i]][bits[b_i]]
                p3_odd[mask] = val % 3

            for c_new in range(3):
                # Pour chaque sommet i < k_prev, (e_i + e_{new}) est de poids 2 => c[i] + c_new + col[i] != 0 mod 3 !
                allowed = [[v for v in range(3) if (c[i] + c_new + v) % 3 != 0] for i in range(k_prev)]
                for emask in range(1 << k_prev):
                    tested += 1
                    col = [allowed[i][(emask >> i) & 1] for i in range(k_prev)]
                    ok = True
                    # Tester tous les masques impairs de taille >= 3 (taille 1 deja garantie par allowed !)
                    for mask in odd_masks:
                        bits = bits_of[mask]
                        if len(bits) == 1:
                            continue
                        s = p3_odd[mask] + c_new
                        for i in bits:
                            s += col[i]
                        if s % 3 == 0:
                            ok = False
                            break
                    if ok:
                        valid += 1
                        new_c = c + (c_new,)
                        new_adj = tuple(adj[i] + (col[i],) for i in range(k_prev)) + (tuple(col) + (0,),)
                        key = canonical_key_p3(k_next, new_c, new_adj)
                        if key not in nxt:
                            nxt[key] = (new_c, new_adj)

        dt = time.time() - t1
        print(f"[k = {k_next}] ({dt:.2f}s) tested={tested}, valid={valid}, classes={len(nxt)}")
        cur = nxt
        if len(cur) == 0:
            print(f"==> EXTINCTION EXACTE A k = {k_next} ! Donc N_gen^(6)((1, 2), OR) = {k_next - 1} EXACT !")
            break

    print(f"Temps total : {time.time() - t0:.2f}s")


if __name__ == "__main__":
    solve_exact_mixed_1_2()
