#!/usr/bin/env python3
"""Full regression with Icarus Verilog (run inside testbench/).

    python3 run_regression.py [--quick] [--jobs N]

For every configuration: generate patterns (gen_pattern.py), compile tb with
iverilog, run, parse the per-case result file and write
regression_summary.md / regression_summary.csv.
"""
import argparse, concurrent.futures as cf, csv, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BANDED = os.path.exists(os.path.join(ROOT, "BSW.py"))
TB = "tb_BSW.v" if BANDED else "tb_SW.v"
WORK = os.path.join("/tmp", "regr_bsw" if BANDED else "regr_sw")


def configs(quick):
    c = []
    if not BANDED:
        # (label, L, NUM_PE, pattern-args, extra defines)
        for p in (256, 128, 64, 32):
            c.append(("regression set", 256, p, None, []))
        n = 40 if quick else 300
        c.append(("random", 256, 256, ["--random", str(n), "--seed", "11"], []))
        c.append(("random", 256, 128, ["--random", str(n), "--seed", "12"], []))
        for L, pes in [(4, (4, 2, 1)), (8, (8, 4, 2, 1)), (16, (16, 8, 4, 2, 1)), (32, (32, 8, 2)),
                       (64, (64, 16, 4)), (100, (100, 50, 25)), (128, (128, 32))]:
            for p in pes:
                c.append(("random", L, p, ["--random", "30" if quick else "150", "--seed", str(L * 1000 + p)], []))
    else:
        # (label, L, B, pattern-args, extra defines)
        for b in (4, 1, 2, 3, 8, 16, 32, 64, 256):
            c.append(("regression set", 256, b, None, []))
        c.append(("regression set, random valid gaps", 256, 4, None, ["RAND_GAP"]))
        n = 40 if quick else 300
        c.append(("random", 256, 4, ["--random", str(n), "--seed", "21"], []))
        c.append(("random", 256, 16, ["--random", str(n // 3), "--seed", "22"], []))
        for L, bs in [(4, (1, 2, 4)), (8, (1, 2, 3, 8)), (16, (1, 2, 4, 5, 16)), (20, (3, 7, 20)),
                      (33, (1, 4, 11)), (64, (2, 8, 13, 64)), (100, (4, 9, 50)), (128, (4, 16))]:
            for b in bs:
                c.append(("random", L, b, ["--random", "30" if quick else "150", "--seed", str(L * 1000 + b)], []))
    return c


def run(cfg):
    label, L, P, pat_args, defs = cfg
    tag = "L%d_%s%d%s" % (L, "B" if BANDED else "PE", P, "_gap" if defs else "")
    tag += "_rand" if pat_args else ""
    d = os.path.join(WORK, tag)
    os.makedirs(d, exist_ok=True)
    gen = [sys.executable, os.path.join(ROOT, "gen_pattern.py"), "--length", str(L), "--out", d]
    if BANDED:
        gen += ["--band", str(P)]
    if pat_args:
        gen += pat_args
    elif L == 256:
        gen += ["--legacy", os.path.join(HERE, "pattern", "case000_s.dat"), os.path.join(HERE, "pattern", "case000_t.dat")]
        if not BANDED:
            gen += [os.path.join(HERE, "pattern", "case001_s.dat"), os.path.join(HERE, "pattern", "case001_t.dat")]
    subprocess.run(gen, check=True, stdout=subprocess.DEVNULL)
    res = os.path.join(d, "tb_result.txt")
    vvp = os.path.join(d, "tb.vvp")
    cmd = ["iverilog", "-g2005", "-DRTL", "-DSEQ_L=%d" % L, "-DPAT_DIR=\"%s/\"" % d,
           "-DRESULT_FILE=\"%s\"" % res, "-DOPS_FILE=\"%s\"" % os.path.join(d, "tb_ops.txt"), "-o", vvp, TB]
    cmd += ["-DBAND=%d" % P] if BANDED else ["-DNUM_PE=%d" % P]
    cmd += ["-D" + x for x in defs]
    subprocess.run(cmd, check=True, cwd=HERE, stderr=subprocess.DEVNULL)
    t0 = time.time()
    out = subprocess.run(["vvp", "-n", vvp], cwd=HERE, capture_output=True, text=True).stdout
    dt = time.time() - t0
    rows = [l.split() for l in open(res) if not l.startswith("#")]
    rows = [list(map(int, r)) for r in rows]
    cases = rows[:-1]                     # last row = re-run after reset recovery
    ok = "ALL" in out and "PASS" in out and "Time out" not in out and "TIMEOUT" not in out
    expected = len([f for f in os.listdir(d) if f.endswith("_golden.dat")])
    assert len(cases) == expected, "%s: only %d of %d cases ran" % (tag, len(cases), expected)
    lat_no_ops = [r[10] - r[8] for r in cases]
    return dict(label=label, L=L, P=P, gap=bool(defs), n=len(cases), passed=sum(r[1] for r in cases),
                reset_ok=ok, ops=sum(r[8] for r in cases), max_score=max(r[2] for r in cases),
                lat=(min(lat_no_ops), max(lat_no_ops)), sim_s=dt, tag=tag)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args()
    cfgs = configs(a.quick)
    results = []
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        for r in ex.map(run, cfgs):
            results.append(r)
            print("%-36s L=%-4d %s=%-4d %4d/%-4d %s" % (r["label"], r["L"], "B" if BANDED else "NUM_PE", r["P"],
                  r["passed"], r["n"], "PASS" if r["reset_ok"] and r["passed"] == r["n"] else "FAIL"), flush=True)
    pname = "B" if BANDED else "NUM_PE"
    with open(os.path.join(HERE, "regression_summary.csv"), "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["config", "L", pname, "cases", "passed", "reset_recovery", "traceback_ops", "max_score",
                    "latency_min", "latency_max"])
        for r in results:
            w.writerow([r["label"], r["L"], r["P"], r["n"], r["passed"], int(r["reset_ok"]), r["ops"], r["max_score"],
                        r["lat"][0], r["lat"][1]])
    tot_c = sum(r["n"] for r in results)
    tot_p = sum(r["passed"] for r in results)
    tot_o = sum(r["ops"] for r in results)
    with open(os.path.join(HERE, "regression_summary.md"), "w") as fp:
        fp.write("| configuration | L | %s | cases | passed | reset recovery | traceback ops checked | latency w/o traceback (cycles) |\n" % pname)
        fp.write("|---|---|---|---|---|---|---|---|\n")
        for r in results:
            lat = "%d" % r["lat"][0] if r["lat"][0] == r["lat"][1] else "%d-%d" % r["lat"]
            fp.write("| %s | %d | %d | %d | %d | %s | %d | %s |\n" % (r["label"], r["L"], r["P"], r["n"], r["passed"],
                     "pass" if r["reset_ok"] else "FAIL", r["ops"], lat))
        fp.write("| **total** | | | **%d** | **%d** | | **%d** | |\n" % (tot_c, tot_p, tot_o))
    print("TOTAL %d / %d cases pass, %d traceback ops compared" % (tot_p, tot_c, tot_o))
    return 0 if tot_p == tot_c and all(r["reset_ok"] for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
