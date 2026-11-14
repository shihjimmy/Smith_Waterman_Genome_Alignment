"""Generate test patterns + golden files for the RTL testbench.

    python3 gen_pattern.py [--band B] [--out testbench/pattern]

Each case NN produces caseNN_s.dat / caseNN_t.dat / caseNN_golden.dat.
golden = hex words: score, end_i, end_j, start_i, start_j, n_ops, ops (end -> start)
"""
import argparse, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from SW import align, check_ops, read_dat, write_dat, write_golden
except ImportError:
    from BSW import align, check_ops, read_dat, write_dat, write_golden

L = 256


def rand_seq(n, rng):
    return [rng.randrange(4) for _ in range(n)]


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


def fit(seq, rng):
    seq = seq[:L]
    return seq + rand_seq(L - len(seq), rng)


def cases(rng, legacy):
    c = list(legacy)
    a = rand_seq(L, rng); c.append(("identical", a, list(a)))
    c.append(("all mismatch", [0] * L, [1] * L))
    c.append(("unrelated random", rand_seq(L, rng), rand_seq(L, rng)))
    a = rand_seq(L, rng); c.append(("5% subst + small indels", a, mutate(a, rng)))
    a = rand_seq(L, rng)
    b = a[:90] + rand_seq(12, rng) + a[90:170] + a[180:]      # long insertion + long deletion
    c.append(("long gaps", a, fit(b, rng)))
    core = rand_seq(100, rng)
    c.append(("local region", fit(rand_seq(60, rng) + core + rand_seq(96, rng), rng),
              fit(rand_seq(30, rng) + mutate(core + [0] * (L - 100), rng, 0.03, 0.02)[:100] + rand_seq(126, rng), rng)))
    a = rand_seq(L, rng); c.append(("20% subst + indels", a, mutate(a, rng, 0.2, 0.06)))
    a = rand_seq(L, rng); c.append(("shifted copy", a, fit(rand_seq(3, rng) + a, rng)))
    a = [rng.choice([0, 0, 0, 1]) for _ in range(L)]
    c.append(("low complexity", a, mutate(a, rng, 0.1, 0.05)))
    return c


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", type=int, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--legacy", nargs="*", default=[], help="pairs of legacy s/t .dat files")
    ap.add_argument("--seed", type=int, default=2023)
    ap.add_argument("--length", type=int, default=256, help="sequence length L")
    ap.add_argument("--random", type=int, default=0, help="generate N random cases only (stress test)")
    args = ap.parse_args()
    L = args.length
    out = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "testbench", "pattern")
    os.makedirs(out, exist_ok=True)
    legacy = []
    for k in range(0, len(args.legacy), 2):
        legacy.append(("legacy " + os.path.basename(args.legacy[k]), read_dat(args.legacy[k]), read_dat(args.legacy[k + 1])))
    rng = random.Random(args.seed)
    if args.random:
        todo = []
        for _ in range(args.random):
            a = rand_seq(L, rng)
            kind = rng.randrange(3)
            if kind == 0:
                todo.append(("random", a, rand_seq(L, rng)))
            else:
                todo.append(("mutated", a, mutate(a, rng, rng.random() * 0.3, rng.random() * 0.2, rng.randint(1, 6))))
    else:
        todo = cases(rng, legacy)
    for n, (name, S, T) in enumerate(todo):
        assert len(S) == L and len(T) == L
        res = align(S, T, args.band)
        check_ops(S, T, res)
        write_dat(os.path.join(out, "case%02d_s.dat" % n), S)
        write_dat(os.path.join(out, "case%02d_t.dat" % n), T)
        write_golden(os.path.join(out, "case%02d_golden.dat" % n), res)
        print("case%02d %-26s score %4d end %s start %s ops %d" % (n, name, res["score"], res["end"], res["start"], len(res["ops"])))
