"""
Golden model for (banded) Smith-Waterman with affine gap penalty and traceback.

This model is bit-exact with the RTL:
  * scoring : match +8, mismatch -5, gap open -7 (first gap base), gap extend -3
              E(i,j) = max(E(i,j-1) - 3, H(i,j-1) - 7)          (gap in S, "INS")
              F(i,j) = max(F(i-1,j) - 3, H(i-1,j) - 7)          (gap in T, "DEL")
              H(i,j) = max(0, H(i-1,j-1) + s(S[i],T[j]), E(i,j), F(i,j))
  * boundary / out-of-band cells: H = 0, E = F = -inf
  * band (Liao et al., "Adaptively Banded Smith-Waterman Algorithm for Long Reads
    and Its Hardware Accelerator", ASAP 2018): with B PEs on the anti-diagonal the
    PE front covers exactly the cells with  -B <= j - i <= B - 1.
  * tie-break rules (shared with the RTL):
      - best cell: highest score, then smallest i, then smallest j
      - H source : H == 0 -> STOP, else DIAG > UP(F) > LEFT(E)
      - E/F      : extend only if strictly better than open
  * traceback emits one op per step, from the alignment end back to its start.

All indices are 0-based sequence indices (cell (i,j) aligns S[i] with T[j]).
"""

MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 8, -5, -7, -3
NEG = -(1 << 30)

STOP, DIAG, UP, LEFT = 0, 1, 2, 3
OP_MATCH, OP_MISMATCH, OP_INS, OP_DEL = 0, 1, 2, 3
OP_NAME = {OP_MATCH: "M", OP_MISMATCH: "X", OP_INS: "I", OP_DEL: "D"}
BASE = "ACGT"


def in_band(i, j, band):
    return band is None or (-band <= j - i <= band - 1)


def align(S, T, band=None):
    """S, T: lists of ints 0..3. band: number of PEs B for banded SW, None = full SW.
    Returns dict(score, end, start, ops) where ops are in emission order (end -> start)."""
    n, m = len(S), len(T)
    H = [[0] * m for _ in range(n)]
    E = [[NEG] * m for _ in range(n)]
    F = [[NEG] * m for _ in range(n)]
    P = [[(STOP, 0, 0)] * m for _ in range(n)]

    def h(i, j):
        return H[i][j] if (i >= 0 and j >= 0 and in_band(i, j, band)) else 0

    def e(i, j):
        return E[i][j] if (i >= 0 and j >= 0 and in_band(i, j, band)) else NEG

    def f(i, j):
        return F[i][j] if (i >= 0 and j >= 0 and in_band(i, j, band)) else NEG

    best = (0, 0, 0)  # (score, i, j)
    for i in range(n):
        for j in range(m):
            if not in_band(i, j, band):
                continue
            e_open, e_extd = h(i, j - 1) + GAP_OPEN, e(i, j - 1) + GAP_EXTEND
            f_open, f_extd = h(i - 1, j) + GAP_OPEN, f(i - 1, j) + GAP_EXTEND
            ev, ex = (e_extd, 1) if e_extd > e_open else (e_open, 0)
            fv, fx = (f_extd, 1) if f_extd > f_open else (f_open, 0)
            d = h(i - 1, j - 1) + (MATCH if S[i] == T[j] else MISMATCH)
            hv = max(0, d, fv, ev)
            if hv == 0:
                src = STOP
            elif hv == d:
                src = DIAG
            elif hv == fv:
                src = UP
            else:
                src = LEFT
            H[i][j], E[i][j], F[i][j], P[i][j] = hv, ev, fv, (src, ex, fx)
            if hv > best[0] or (hv == best[0] and (i, j) < (best[1], best[2])):
                best = (hv, i, j)

    # ---------------- traceback ----------------
    score, i, j = best
    state, ops = "H", []
    while True:
        if state == "H" and (i < 0 or j < 0 or not in_band(i, j, band)):
            break
        src, ex, fx = P[i][j]
        act = state if state != "H" else {STOP: None, DIAG: "D", UP: "F", LEFT: "E"}[src]
        if act is None:
            break
        if act == "D":
            ops.append(OP_MATCH if S[i] == T[j] else OP_MISMATCH)
            i, j, state = i - 1, j - 1, "H"
        elif act == "F":
            ops.append(OP_DEL)
            i, state = i - 1, ("F" if fx else "H")
        else:
            ops.append(OP_INS)
            j, state = j - 1, ("E" if ex else "H")
    return dict(score=score, end=(best[1], best[2]), start=(i + 1, j + 1), ops=ops)


def check_ops(S, T, res):
    """Re-score the alignment from its ops; must equal the reported score."""
    i, j = res["start"]
    sc, prev = 0, None
    for op in reversed(res["ops"]):
        if op in (OP_MATCH, OP_MISMATCH):
            assert (S[i] == T[j]) == (op == OP_MATCH)
            sc += MATCH if op == OP_MATCH else MISMATCH
            i, j = i + 1, j + 1
        elif op == OP_DEL:
            sc += GAP_EXTEND if prev == OP_DEL else GAP_OPEN
            i += 1
        else:
            sc += GAP_EXTEND if prev == OP_INS else GAP_OPEN
            j += 1
        prev = op
    if res["ops"]:
        assert (i - 1, j - 1) == res["end"], ((i - 1, j - 1), res["end"])
    assert sc == res["score"], (sc, res["score"])
    return True


def pretty(S, T, res):
    i, j = res["start"]
    a, b, mid = "", "", ""
    for op in reversed(res["ops"]):
        if op in (OP_MATCH, OP_MISMATCH):
            a += BASE[S[i]]; b += BASE[T[j]]; mid += "|" if op == OP_MATCH else "."
            i += 1; j += 1
        elif op == OP_DEL:
            a += BASE[S[i]]; b += "-"; mid += " "; i += 1
        else:
            a += "-"; b += BASE[T[j]]; mid += " "; j += 1
    return a, mid, b


def read_dat(path):
    return [int(x, 16) for x in open(path).read().split()]


def write_dat(path, seq):
    with open(path, "w") as fp:
        fp.write("\n".join("%x" % v for v in seq) + "\n")


def write_golden(path, res):
    """hex words: score, end_i, end_j, start_i, start_j, n_ops, ops..."""
    words = [res["score"], res["end"][0], res["end"][1], res["start"][0], res["start"][1],
             len(res["ops"])] + res["ops"]
    with open(path, "w") as fp:
        fp.write("\n".join("%x" % w for w in words) + "\n")


def from_str(x):
    return [BASE.index(c) for c in x]


if __name__ == "__main__":
    import argparse, time
    ap = argparse.ArgumentParser(description="Smith-Waterman golden model (affine gap, traceback)")
    ap.add_argument("--s", help="S sequence .dat file (hex, one base per line)")
    ap.add_argument("--t", help="T sequence .dat file")
    ap.add_argument("--band", type=int, default=None, help="banded SW with B PEs (default: full SW)")
    args = ap.parse_args()
    if args.s and args.t:
        S, T = read_dat(args.s), read_dat(args.t)
    else:
        S = from_str("AAAGGCTGGGGACCACTGATCTAAATACACCAATAAAAAGAAAAAGATTGTAAGATTGGAGTTTAAAAGACCTGACTCTATACTGACCACAAATAAAAACC")
        T = from_str("AAAGGCTGGGGACCCTGATCTAAATACACCAATAAAAAGAAAAAGATTGTAAGATTGGATTTTAAAAGACCTGACTCTATACTGACCACAAAAAAACC")
    t0 = time.time()
    res = align(S, T, args.band)
    check_ops(S, T, res)
    print("execution time (ms): %.2f" % (1000 * (time.time() - t0)))
    print("score %d  S[%d..%d]  T[%d..%d]  %d ops" % (res["score"], res["start"][0], res["end"][0],
          res["start"][1], res["end"][1], len(res["ops"])))
    for line in pretty(S, T, res):
        print(line)
