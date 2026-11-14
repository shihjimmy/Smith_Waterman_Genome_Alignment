#!/bin/sh
# Open-source flow (Icarus Verilog >= 11). Run inside testbench/.
#   sh run_iverilog.sh            -> NUM_PE = 256 and NUM_PE = 128 on testbench/pattern
#   sh run_iverilog.sh stress     -> additionally random stress tests on small L / many NUM_PE
set -e
run() {   # $1 = label, rest = iverilog defines
    label=$1; shift
    iverilog -g2005 -DRTL "$@" -o /tmp/tb_SW.vvp tb_SW.v
    printf '%-28s ' "$label"
    vvp -n /tmp/tb_SW.vvp | grep -E 'ALL .* PASS|FAIL :|Time out' || echo "no result"
}
run "L=256 NUM_PE=256"
run "L=256 NUM_PE=128 (SW_v2)" -DSW_v2
if [ "$1" = "stress" ]; then
    for cfg in "16 16" "16 8" "16 4" "16 1" "32 2" "64 64" "64 16"; do
        set -- $cfg
        d=/tmp/sw_stress_$1/
        python3 ../gen_pattern.py --length $1 --random 60 --seed $1 --out $d > /dev/null
        run "L=$1 NUM_PE=$2 (60 random)" -DSEQ_L=$1 -DNUM_PE=$2 -DPAT_DIR=\"$d\" -DNCASE=60
    done
fi
