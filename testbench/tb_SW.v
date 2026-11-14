`timescale 1ns/10ps
//============================================================================
// Testbench for SW (Smith-Waterman, affine gap, traceback)
//
//   +define+RTL            : simulate ../verilog/SW.v
//   +define+SDF            : simulate gate-level netlist ../verilog/SW_syn.v
//   +define+SW_v2          : NUM_PE = 128 (2 passes + recirculation buffer)
//                            (default NUM_PE = 256, one pass)
//   +define+FSDB           : dump waveform (Verdi)
//   +define+VCD=\"x.vcd\"   : dump VCD (open-source viewers, e.g. GTKWave)
//   (input must be contiguous: S[p] is broadcast to PE p exactly when PE p starts)
//   +define+SEQ_L=<n> +define+NUM_PE=<p> +define+PAT_DIR=\"dir/\" +define+NCASE=<n>
//                          : other configurations (see gen_pattern.py --length)
//
// Runs NCASE patterns back-to-back (without reset in between) and checks
// score, alignment end/start, alignment length and every traceback op
// against the golden model (SW.py / gen_pattern.py).
//============================================================================
`define clk_period 12
`ifndef END_CYCLE
    `define END_CYCLE 50000000   // global watchdog (each case also has its own)
`endif
`ifndef NCASE
    `define NCASE 1000         // upper limit; the case list ends at the first missing caseNNN_golden.dat
`endif
`ifndef SEQ_L
    `define SEQ_L 256
`endif
`ifndef RESULT_FILE
    `define RESULT_FILE "tb_result.txt"
`endif
`ifndef OPS_FILE
    `define OPS_FILE "tb_ops.txt"
`endif
`ifndef PAT_DIR
    `define PAT_DIR "pattern/"
`endif
`ifdef SW_v2
    `define NUM_PE 128
    `define SDFFILE "../verilog/SW_v2_syn.sdf"
`else
    `ifndef NUM_PE
        `define NUM_PE 256
    `endif
    `define SDFFILE "../verilog/SW_syn.sdf"
`endif

`ifdef SDF
    `include "../verilog/tsmc13.v"
    `ifdef SW_v2
        `include "../verilog/SW_v2_syn.v"
    `else
        `include "../verilog/SW_syn.v"
    `endif
`else
    `include "../verilog/SW.v"
`endif

module tb_SW;
    localparam L  = `SEQ_L;
    localparam IW = (L > 1) ? $clog2(L) : 1;
    localparam DW = $clog2(8 * L + 1) + 1;

    reg           clk, reset, valid;
    reg  [1:0]    data_s, data_t;
    wire          finish, aln_valid;
    wire [1:0]    aln_op;
    wire [DW-1:0] max;
    wire [IW-1:0] end_i, end_j, start_i, start_j;
    wire [IW+1:0] aln_len;

    reg  [1:0]    seq_s [0:L-1];
    reg  [1:0]    seq_t [0:L-1];
    reg  [15:0]   gold  [0:2*L+8];
    reg  [1:0]    ops   [0:2*L+8];
    integer       n_ops, ncase, i, err, case_err, total_err, fd, fd_ops, fchk, gap_seed, n_pass, n_case;
    reg  [8*256-1:0] fname;

`ifdef SDF
    SW SW0 (
`else
    SW #(.L(L), .NUM_PE(`NUM_PE)) SW0 (
`endif
        .clk(clk), .reset(reset), .valid(valid), .data_s(data_s), .data_t(data_t),
        .finish(finish), .max(max), .end_i(end_i), .end_j(end_j),
        .start_i(start_i), .start_j(start_j), .aln_len(aln_len),
        .aln_valid(aln_valid), .aln_op(aln_op));

    initial begin
`ifdef SDF
        $sdf_annotate(`SDFFILE, SW0);
`endif
`ifdef VCD
        $dumpfile(`VCD);
        $dumpvars(0, tb_SW);
`endif
`ifdef FSDB
        $fsdbDumpfile("SW.fsdb");
        $fsdbDumpvars(0, tb_SW, "+mda");
`endif
    end

    always #(`clk_period/2) clk = ~clk;

    // collect traceback ops
    always @(negedge clk) begin
        if (aln_valid) begin
            ops[n_ops] = aln_op;
            n_ops = n_ops + 1;
        end
    end

    task run_case(input integer c);
        integer cycles;
        begin
            $sformat(fname, "%scase%03d_s.dat", `PAT_DIR, c);      $readmemh(fname, seq_s);
            $sformat(fname, "%scase%03d_t.dat", `PAT_DIR, c);      $readmemh(fname, seq_t);
            for (i = 0; i < 2*L+9; i = i + 1) gold[i] = 16'hxxxx;
            $sformat(fname, "%scase%03d_golden.dat", `PAT_DIR, c); $readmemh(fname, gold);
            n_ops = 0;

            // stream S and T (optionally with random idle cycles between inputs)
            for (i = 0; i < L; i = i + 1) begin
`ifdef RAND_GAP
                while (($random(gap_seed) & 3) == 0) begin
                    @(negedge clk);
                    valid  = 1'b0;
                    data_s = $random(gap_seed);
                    data_t = $random(gap_seed);
                end
`endif
                @(negedge clk);
                valid  = 1'b1;
                data_s = seq_s[i];
                data_t = seq_t[i];
            end
            @(negedge clk);
            valid  = 1'b0;
            data_s = 2'd0;
            data_t = 2'd0;

            // wait for finish
            cycles = 0;
            while (finish !== 1'b1 && cycles < 16 * L + 64) begin
                @(negedge clk);
                cycles = cycles + 1;
            end
            if (finish !== 1'b1) begin
                $display("case%03d TIMEOUT: no finish within %0d cycles", c, cycles);
                total_err = total_err + 1;
            end

            case_err = 0;
            if (max     !== gold[0][DW-1:0]) case_err = case_err + 1;
            if (end_i   !== gold[1][IW-1:0]) case_err = case_err + 1;
            if (end_j   !== gold[2][IW-1:0]) case_err = case_err + 1;
            if (start_i !== gold[3][IW-1:0]) case_err = case_err + 1;
            if (start_j !== gold[4][IW-1:0]) case_err = case_err + 1;
            if (aln_len !== gold[5][IW+1:0] || n_ops != gold[5]) case_err = case_err + 1;
            err = 0;
            for (i = 0; i < gold[5] && i < n_ops; i = i + 1)
                if (ops[i] !== gold[6+i][1:0]) err = err + 1;
            case_err  = case_err + err;
            total_err = total_err + case_err;
            if (case_err == 0) n_pass = n_pass + 1;
            // hardware traceback ops (end -> start), for plotting the alignment path
            $fwrite(fd_ops, "%0d %0d %0d %0d %0d ", c, end_i, end_j, start_i, start_j);
            for (i = 0; i < n_ops; i = i + 1) $fwrite(fd_ops, "%0d", ops[i]);
            $fwrite(fd_ops, "\n");
            $fdisplay(fd, "%0d %0d %0d %0d %0d %0d %0d %0d %0d %0d %0d", c, (case_err == 0),
                      max, gold[0], end_i, end_j, start_i, start_j, n_ops, err, cycles);

            $display("case%03d %s  score %4d/%4d  end (%3d,%3d)/(%3d,%3d)  start (%3d,%3d)/(%3d,%3d)  ops %3d/%3d (%0d op errors)  latency %0d cycles",
                     c, (case_err == 0) ? "PASS" : "FAIL",
                     max, gold[0], end_i, end_j, gold[1], gold[2],
                     start_i, start_j, gold[3], gold[4], n_ops, gold[5], err, cycles);
            repeat (3) @(negedge clk);
        end
    endtask

    initial begin
        clk = 1'b0; reset = 1'b0; valid = 1'b0; data_s = 2'd0; data_t = 2'd0;
        total_err = 0;
        n_pass    = 0;
        gap_seed  = 1;
        fd = $fopen(`RESULT_FILE);
        fd_ops = $fopen(`OPS_FILE);
        $fdisplay(fd, "# case pass score golden end_i end_j start_i start_j n_ops op_err latency");
        #(`clk_period)   reset = 1'b1;
        #(`clk_period*2) reset = 1'b0;
        #(`clk_period*2);

        $display("------------------------------------------------------------");
        $display(" SW  L = %0d  NUM_PE = %0d", L, `NUM_PE);
        $display("            yours/golden");
        $display("------------------------------------------------------------");
        n_case = 0;
        for (ncase = 0; ncase < `NCASE; ncase = ncase + 1) begin
            $sformat(fname, "%scase%03d_golden.dat", `PAT_DIR, ncase);
            fchk = $fopen(fname, "r");
            if (fchk == 0)
                ncase = `NCASE;                 // no more cases
            else begin
                $fclose(fchk);
                run_case(ncase);
                n_case = n_case + 1;
            end
        end

        // reset-recovery test: abort a run in the middle of the scoring phase,
        // check that nothing comes out, then run case 0 again
        $display("--- reset in the middle of scoring, then re-run case00 ---");
        for (i = 0; i < L; i = i + 1) begin
            @(negedge clk);
            valid = 1'b1; data_s = seq_s[i]; data_t = seq_t[i];
        end
        @(negedge clk) valid = 1'b0;
        repeat (L / 2) @(negedge clk);
        #(`clk_period/4) reset = 1'b1;
        #(`clk_period)   reset = 1'b0;
        repeat (4 * L) begin
            @(negedge clk);
            if (finish === 1'b1 || aln_valid === 1'b1) total_err = total_err + 1;
        end
        run_case(0);

        $display("------------------------------------------------------------");
        if (total_err == 0)
            $display(" ALL %0d CASES PASS (+ reset recovery)", n_case);
        else
            $display(" FAIL : %0d mismatches, %0d / %0d runs pass", total_err, n_pass, n_case + 1);
        $fclose(fd);
        $fclose(fd_ops);
        $display("------------------------------------------------------------");
        $finish;
    end

    initial begin
        #(`clk_period*`END_CYCLE);
        $display("Time out!!");
        $finish;
    end
endmodule
