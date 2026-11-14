"""Independent check of the golden model.

  * full SW score  vs  Biopython PairwiseAligner (local, match 8, mismatch -5,
    open_gap_score -7, extend_gap_score -3)           pip install biopython
  * banded score    vs  a second, independently written banded implementation
  * every traceback is re-scored (check_ops) and must give the reported score

    python3 tools/validate_model.py [N]
"""
import os, sys, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from SW import align, check_ops
from Bio.Align import PairwiseAligner
al = PairwiseAligner(mode="local", match_score=8, mismatch_score=-5, open_gap_score=-7, extend_gap_score=-3)
B4 = "ACGT"
NEG = -10**9
def band_ref(S, T, B):          # independent banded implementation (dict based, column-major)
    H, E, F = {}, {}, {}
    g = lambda d, i, j, dv: d.get((i, j), dv)
    best = 0
    for j in range(len(T)):
        for i in range(len(S)):
            if not (-B <= j - i <= B - 1): continue
            E[i, j] = max(g(E, i, j-1, NEG) - 3, g(H, i, j-1, 0) - 7)
            F[i, j] = max(g(F, i-1, j, NEG) - 3, g(H, i-1, j, 0) - 7)
            H[i, j] = max(0, g(H, i-1, j-1, 0) + (8 if S[i] == T[j] else -5), E[i, j], F[i, j])
            best = max(best, H[i, j])
    return best
rng = random.Random(7)
n_full = n_band = 0
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
for it in range(N):
    L = rng.choice([1, 2, 3, 5, 8, 16, 33, 64, 128, 256]) if it < N - 200 else 256
    m = L if rng.random() < .8 else rng.randint(1, 64)
    S = [rng.randrange(4) for _ in range(L)]
    if rng.random() < .5:
        T = [rng.randrange(4) for _ in range(m)]
    else:   # mutated copy
        T = [b if rng.random() > .15 else rng.randrange(4) for b in S]
        for _ in range(rng.randint(0, 4)):
            p = rng.randrange(len(T) + 1); k = rng.randint(1, 8)
            T = T[:p] + [rng.randrange(4) for _ in range(k)] + T[p:] if rng.random() < .5 else T[:p] + T[p+k:]
        T = T or [0]
    r = align(S, T); check_ops(S, T, r)
    bio = al.score("".join(B4[x] for x in S), "".join(B4[x] for x in T))
    bio = max(0, int(bio))
    assert r["score"] == bio, (r["score"], bio, L)
    n_full += 1
    if len(S) == len(T) and L <= 128:
        for B in {1, 2, 3, 4, max(1, L // 3), L}:
            rb = align(S, T, B); check_ops(S, T, rb)
            assert rb["score"] == band_ref(S, T, B)
            n_band += 1
print("full SW vs Biopython:", n_full, "pairs OK;  banded vs independent impl:", n_band, "runs OK")
