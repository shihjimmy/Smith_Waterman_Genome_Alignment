"""Generate test patterns + golden files for the RTL testbench.

    python3 gen_pattern.py [--band B] [--length L] [--out DIR]          # regression set
    python3 gen_pattern.py --random N [--length L] [--seed S] --out DIR  # random stress set

Each case NNN produces caseNNN_s.dat / caseNNN_t.dat / caseNNN_golden.dat and a
line in DIR/cases.txt.  golden = hex words:
    score, end_i, end_j, start_i, start_j, n_ops, ops (end -> start)
ops: 0 MATCH, 1 MISMATCH, 2 INS (gap in S), 3 DEL (gap in T)

Directed cases use disjoint background alphabets (S from {A,C}, T from {G,T}) so
that the optimum is controlled by the planted motifs.
"""
import argparse, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from SW import align, check_ops, read_dat, write_dat, write_golden
    DEFAULT_BAND = None                       # SW repository: full Smith-Waterman
except ImportError:
    from BSW import align, check_ops, read_dat, write_dat, write_golden
    DEFAULT_BAND = 4                          # BSW repository: B = 4 PEs

L = 256


def rand_seq(n, rng, alphabet=(0, 1, 2, 3)):
    return [rng.choice(alphabet) for _ in range(n)]


def fit(seq, rng):
    seq = seq[:L]
    return seq + rand_seq(L - len(seq), rng)


def mutate(seq, rng, sub=0.05, indel=0.02, max_indel=3):
    out, k = [], 0
    while k < len(seq):
        r = rng.random()
        if r < indel / 2:                                   # deletion
            k += rng.randint(1, max_indel)
            continue
        if r < indel:                                       # insertion
            out += rand_seq(rng.randint(1, max_indel), rng)
        b = seq[k]
        if rng.random() < sub:
            b = (b + rng.randint(1, 3)) % 4
        out.append(b)
        k += 1
    return fit(out, rng)


# ---------------------------------------------------------------- random-like cases
def typical_cases(rng):
    c = []
    a = rand_seq(L, rng); c.append(("identical", a, list(a)))
    c.append(("all mismatch", [0] * L, [1] * L))
    c.append(("unrelated random", rand_seq(L, rng), rand_seq(L, rng)))
    a = rand_seq(L, rng); c.append(("5% subst + small indels", a, mutate(a, rng)))
    a = rand_seq(L, rng)
    c.append(("long gaps", a, fit(a[:L * 90 // 256] + rand_seq(L * 12 // 256, rng)
                                  + a[L * 90 // 256:L * 170 // 256] + a[L * 180 // 256:], rng)))
    core = rand_seq(L * 100 // 256, rng)
    c.append(("local region", fit(rand_seq(L * 60 // 256, rng) + core, rng),
              fit(rand_seq(L * 30 // 256, rng) + mutate(core + [0] * (L - len(core)), rng, 0.03, 0.02)[:len(core)], rng)))
    a = rand_seq(L, rng); c.append(("20% subst + indels", a, mutate(a, rng, 0.2, 0.06)))
    a = rand_seq(L, rng); c.append(("shifted copy", a, fit(rand_seq(3, rng) + a, rng)))
    a = [rng.choice([0, 0, 0, 1]) for _ in range(L)]
    c.append(("low complexity", a, mutate(a, rng, 0.1, 0.05)))
    return c


# ---------------------------------------------------------------- directed cases
def bg(rng):
    return rand_seq(L, rng, (0, 1)), rand_seq(L, rng, (2, 3))


def plant(seq, pos, motif):
    seq = list(seq)
    for k, b in enumerate(motif):
        if 0 <= pos + k < L:
            seq[pos + k] = b
    return seq


def directed_cases(rng, band):
    c = []
    m = max(4, L // 8)                        # motif length (32 for L = 256)
    q = L // 4

    def planted(name, ps, pt, motif=None, tmotif=None):
        motif = motif or rand_seq(m, rng)
        S, T = bg(rng)
        c.append((name, plant(S, ps, motif), plant(T, pt, tmotif or motif)))

    planted("motif in the middle", q, q + 3)
    planted("motif at the start (0,0)", 0, 0)
    planted("motif at the end (L-1,L-1)", L - m, L - m)
    planted("motif S start / T end", 0, L - m)
    planted("motif S end / T start", L - m, 0)
    S, T = bg(rng)
    c.append(("single matching base", plant(S, L // 3, [2]), plant(T, 2 * L // 3, [2])))

    mo = rand_seq(m, rng)
    S, T = bg(rng)
    c.append(("tie: two equal hits in S", plant(plant(S, q // 2, mo), q // 2 + L // 2, mo), plant(T, q // 2, mo)))
    S, T = bg(rng)
    c.append(("tie: two equal hits in T", plant(S, q // 2, mo), plant(plant(T, q // 2, mo), q // 2 + L // 2, mo)))
    S, T = bg(rng)
    c.append(("tie: 2x2 equal hits", plant(plant(S, 10, mo), L // 2, mo), plant(plant(T, 10, mo), L // 2, mo)))

    g = max(2, L // 12)                       # long gap (21 for L = 256)
    m1, m2 = rand_seq(m + 8, rng), rand_seq(m + 8, rng)
    S, T = bg(rng)
    c.append(("long INS (gap in S)", plant(S, q, m1 + m2), plant(T, q, m1 + rand_seq(g, rng, (2, 3)) + m2)))
    S, T = bg(rng)
    c.append(("long DEL (gap in T)", plant(S, q, m1 + rand_seq(g, rng, (0, 1)) + m2), plant(T, q, m1 + m2)))
    base = rand_seq(L // 2, rng)
    t = list(base)
    for pos, kind in [(10, "I"), (25, "D"), (40, "I"), (41, "I"), (60, "D"), (61, "D"), (80, "I"), (95, "D")]:
        pos = pos * L // 256
        if pos < len(t):
            t = t[:pos] + [rng.randrange(4)] + t[pos:] if kind == "I" else t[:pos] + t[pos + 1:]
    S, T = bg(rng)
    c.append(("mixed short indels", plant(S, q // 2, base), plant(T, q // 2, t)))
    mm = rand_seq(m, rng)
    mm_t = [b if k % 3 else (b + 1) % 4 for k, b in enumerate(mm)]
    S, T = bg(rng)
    c.append(("every 3rd base mismatch", plant(S, q, mm), plant(T, q, mm_t)))
    S, T = bg(rng)
    c.append(("bridge over 3 mismatches", plant(S, q, m1 + [0, 0, 0] + m2), plant(T, q, m1 + [3, 3, 3] + m2)))
    c.append(("all same base", [0] * L, [0] * L))
    c.append(("periodic AC vs CA", [k % 2 for k in range(L)], [(k + 1) % 2 for k in range(L)]))
    a = rand_seq(L, rng)
    c.append(("identical core, mutated ends", a, [(x + 1) % 4 if (k < L // 32 or k >= L - L // 32) else x
                                                  for k, x in enumerate(a)]))

    if band:                                  # band-edge cases (j - i offsets)
        for off, tag in [(band - 1, "inside, j-i = B-1"), (-band, "inside, j-i = -B"),
                         (band, "outside, j-i = B"), (-band - 1, "outside, j-i = -B-1")]:
            ps, pt = (q, q + off) if off >= 0 else (q - off, q)
            planted("band edge " + tag, ps, pt)
        a = rand_seq(L, rng)
        t = list(a)
        for pos in range(L // 8, L, L // 8):    # one inserted base every L/8 -> drift of 7
            t.insert(pos, rng.randrange(4))
        c.append(("indel drift beyond band", a, fit(t, rng)))
        S, T = bg(rng)
        gb = band + 2
        c.append(("gap longer than band", plant(S, q, m1 + m2), plant(T, q, m1 + rand_seq(gb, rng, (2, 3)) + m2)))
    return c


def emit(out, todo, band):
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "cases.txt"), "w") as idx:
        for n, (name, S, T) in enumerate(todo):
            assert len(S) == L and len(T) == L, name
            res = align(S, T, band)
            check_ops(S, T, res)
            write_dat(os.path.join(out, "case%03d_s.dat" % n), S)
            write_dat(os.path.join(out, "case%03d_t.dat" % n), T)
            write_golden(os.path.join(out, "case%03d_golden.dat" % n), res)
            line = "case%03d %-32s score %4d end %-10s start %-10s ops %d" % (
                n, name, res["score"], "(%d,%d)" % res["end"], "(%d,%d)" % res["start"], len(res["ops"]))
            idx.write(line + "\n")
            print(line)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", type=int, default=DEFAULT_BAND, help="banded SW with B PEs (0 = full SW)")
    ap.add_argument("--length", type=int, default=256, help="sequence length L")
    ap.add_argument("--random", type=int, default=0, help="generate N random cases only (stress test)")
    ap.add_argument("--seed", type=int, default=2023)
    ap.add_argument("--pairs", nargs="*", default=[], help="pairs of s/t .dat files to include first")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    L = args.length
    band = args.band or None
    out = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "testbench", "pattern")
    rng = random.Random(args.seed)
    if args.random:
        todo = []
        for _ in range(args.random):
            a = rand_seq(L, rng)
            kind = rng.randrange(4)
            if kind == 0:
                todo.append(("random", a, rand_seq(L, rng)))
            elif kind == 1 and L >= 8:
                todo.append(("random directed", *directed_cases(rng, band)[rng.randrange(17)][1:]))
            else:
                todo.append(("mutated", a, mutate(a, rng, rng.random() * 0.3, rng.random() * 0.2, rng.randint(1, 6))))
    else:
        todo = [("sample " + os.path.basename(args.pairs[k]), read_dat(args.pairs[k]), read_dat(args.pairs[k + 1]))
                for k in range(0, len(args.pairs), 2)]
        todo += typical_cases(rng) + directed_cases(rng, band)
    emit(out, todo, band)
