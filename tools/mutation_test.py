"""Mutation test of the testbench: inject small RTL bugs and check that the
regression set (testbench/pattern) detects every one of them.

    python3 tools/mutation_test.py            (SW: also run with SWV2=1 for NUM_PE = 128)
"""
import subprocess, sys, os
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
repo, tb, src = ROOT, "tb_SW.v", os.path.join(ROOT, "verilog", "SW.v")
SW = [("tie-break: E/F extend on equal (> -> >=)", "wire                 e_ext  = (e_extv > e_open);", "wire                 e_ext  = (e_extv >= e_open);"),
      ("H source priority: E checked before F", "if (f_new > h_new) begin h_new = f_new; hsrc = UP;   end\n        if (e_new > h_new) begin h_new = e_new; hsrc = LEFT; end", "if (e_new > h_new) begin h_new = e_new; hsrc = LEFT; end\n        if (f_new > h_new) begin h_new = f_new; hsrc = UP;   end"),
      ("best tie-break: larger i/j wins", "better = (a[2*IW-1:0] < b[2*IW-1:0]);     // smaller i, then smaller j", "better = (a[2*IW-1:0] > b[2*IW-1:0]);"),
      ("no boundary reset at column 0 (diag)", "wire signed [DW-1:0] h_diag = first ? {DW{1'b0}} : hd_r;", "wire signed [DW-1:0] h_diag = hd_r;"),
      ("gap-extend penalty -2 instead of -3", "GAP_EXT = -3;", "GAP_EXT = -2;"),
      ("traceback: ignore f_ext", "st    <= f_ext ? ST_F : ST_H;", "st    <= ST_H;"),
      ("recirculation buffer one stage short", "localparam D = L - NUM_PE;", "localparam D = L - NUM_PE - 1;")]
BSW = [("diag neighbour always PE k", "assign dg_v = (!mv_down && !mv_down_p) ? v2[p-1] : v2[p];", "assign dg_v = v2[p];"),
       ("left neighbour H after DOWN = PE k", "assign lf_h = mv_down ? h1[p+1] : h1[p];", "assign lf_h = h1[p];"),
       ("schedule: right move first", "r = (m + 1) / 2;", "r = m / 2;"),
       ("tie-break: E/F extend on equal", "wire                 e_ext  = (e_extv > e_open);", "wire                 e_ext  = (e_extv >= e_open);"),
       ("traceback ignores band check", "if (oob || !inband || hsrc == STOP) stop = 1'b1;", "if (oob || hsrc == STOP) stop = 1'b1;"),
       ("best tie-break: larger i/j wins", "better = (a[2*IW-1:0] < b[2*IW-1:0]);", "better = (a[2*IW-1:0] > b[2*IW-1:0]);"),
       ("traceback: ignore e_ext", "st    <= e_ext ? ST_E : ST_H;", "st    <= ST_H;")]
muts = SW if os.path.basename(src) == "SW.v" else BSW
if os.environ.get("ONLY"):
    muts = [m for m in muts if m[0] in os.environ["ONLY"].split("|")]
orig = open(src).read()
for k, (name, a, b) in enumerate(muts):
    assert orig.count(a) >= 1, name
    d = "/tmp/mutation_test/%s_%d" % (os.path.basename(src)[:-2], k); os.makedirs(d, exist_ok=True)
    open(d + "/dut.v", "w").write(orig.replace(a, b))
    t = open(os.path.join(repo, "testbench", tb)).read().replace('"../verilog/%s"' % os.path.basename(src), '"%s/dut.v"' % d)
    open(d + "/tb.v", "w").write(t)
    pat = os.path.join(repo, "testbench", "pattern") + "/"
    extra = ["-DSW_v2"] if os.environ.get("SWV2") else []
    subprocess.run(["iverilog", "-g2005", "-DRTL"] + extra + [ "-DPAT_DIR=\"%s\"" % pat, "-DRESULT_FILE=\"%s/r.txt\"" % d,
                    "-DOPS_FILE=\"%s/o.txt\"" % d, "-o", d + "/a.vvp", d + "/tb.v"], check=True, stderr=subprocess.DEVNULL)
    out = subprocess.run(["vvp", "-n", d + "/a.vvp"], capture_output=True, text=True).stdout
    ncase = sum(1 for l in out.splitlines() if l.startswith("case"))
    nfail = sum(1 for l in out.splitlines() if l.startswith("case") and "FAIL" in l)
    res = [l for l in out.splitlines() if "ALL" in l or "FAIL :" in l or "Time out" in l]
    print("%-44s -> %d of %d runs fail | %s" % (name, nfail, ncase, res[-1].strip() if res else "no result"), flush=True)
