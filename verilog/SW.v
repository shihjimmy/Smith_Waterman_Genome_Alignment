//============================================================================
// Smith-Waterman local alignment accelerator (affine gap) with traceback
//
//   * Linear systolic array of NUM_PE processing elements.
//       NUM_PE == L      : one PE per row of S, a single pass
//       NUM_PE == L/2 ...: the array is reused for L/NUM_PE passes;
//                          the bottom row of one pass is fed back to
//                          PE0 through a (L-NUM_PE)-deep buffer
//   * Scoring : match +8, mismatch -5, gap open -7, gap extend -3
//       E(i,j) = max(E(i,j-1)-3, H(i,j-1)-7)
//       F(i,j) = max(F(i-1,j)-3, H(i-1,j)-7)
//       H(i,j) = max(0, H(i-1,j-1)+s, E(i,j), F(i,j))
//   * Every cell writes a 4-bit traceback pointer {hsrc[1:0], e_ext, f_ext}
//     into the pointer bank of its PE. After the scoring phase a traceback
//     FSM walks back from the best cell, emitting one alignment op per cycle.
//
// Interface
//   valid/data_s/data_t : S[k] and T[k] for k = 0..L-1 on L consecutive cycles
//   aln_valid/aln_op    : alignment ops, from the alignment END back to its START
//                         0 = MATCH, 1 = MISMATCH, 2 = INS (gap in S), 3 = DEL (gap in T)
//   finish (1 cycle)    : max = best score, (end_i,end_j) = last aligned pair,
//                         (start_i,start_j) = first aligned pair, aln_len = #ops
//   Tie-break: best cell = highest score, then smallest i, then smallest j;
//              H source priority DIAG > UP(F) > LEFT(E); gaps extend only if
//              strictly better than open.
//============================================================================
module SW #(
    parameter L      = 256,               // sequence length
    parameter NUM_PE = 256,               // number of PEs, must divide L (power of 2 recommended)
    // derived widths -- do not override
    parameter IW     = (L > 1) ? $clog2(L) : 1,       // index width
    parameter DW     = $clog2(8 * L + 1) + 1          // signed score width (8L fits)
)(
    input                    clk,
    input                    reset,
    input                    valid,
    input      [1:0]         data_s,
    input      [1:0]         data_t,
    output reg               finish,
    output reg [DW-1:0]      max,
    output reg [IW-1:0]      end_i,
    output reg [IW-1:0]      end_j,
    output reg [IW-1:0]      start_i,
    output reg [IW-1:0]      start_j,
    output reg [IW+1:0]      aln_len,
    output                   aln_valid,
    output     [1:0]         aln_op
);
    localparam Q  = L / NUM_PE;                     // number of passes
    localparam QW = (Q > 1) ? $clog2(Q) : 1;
    localparam BW = DW + 2 * IW;                    // {score, i, j}
    localparam TW = 1 + 1 + 1 + 2 + IW + QW + DW + DW;   // token width
    localparam signed [DW-1:0] NEG = -16;

    localparam S_IDLE = 2'd0, S_RUN = 2'd1, S_TB = 2'd2, S_DONE = 2'd3;

    reg  [1:0]    state_r;
    reg  [IW:0]   in_cnt_r;
    reg  [1:0]    seq_s [0:L-1];
    reg  [1:0]    seq_t [0:L-1];

    // ------------------------------------------------------------------
    // input capture
    // ------------------------------------------------------------------
    wire accept = valid && (state_r == S_IDLE || state_r == S_RUN) && (in_cnt_r < L);
    always @(posedge clk or posedge reset) begin
        if (reset)
            in_cnt_r <= 0;
        else if (state_r == S_DONE)
            in_cnt_r <= 0;
        else if (accept)
            in_cnt_r <= in_cnt_r + 1'b1;
    end

    // sequence storage (plain memory, no reset needed)
    always @(posedge clk) begin
        if (accept) begin
            seq_s[in_cnt_r[IW-1:0]] <= data_s;
            seq_t[in_cnt_r[IW-1:0]] <= data_t;
        end
    end

    // ------------------------------------------------------------------
    // token stream:  {valid, first, last, t, col, q, h_up, f_up}
    // ------------------------------------------------------------------
    wire [TW-1:0] tok   [0:NUM_PE];     // tok[p] = input of PE p
    wire [BW-1:0] best  [0:NUM_PE];     // best[p] = best-so-far input of PE p
    wire [3:0]    rd_ptr[0:NUM_PE-1];
    wire [TW-1:0] recirc_tok;

    wire [TW-1:0] in_tok = {1'b1, (in_cnt_r == 0), (in_cnt_r == L - 1), data_t,
                            in_cnt_r[IW-1:0], {QW{1'b0}}, {DW{1'b0}}, NEG};

    assign tok[0]  = accept ? in_tok : recirc_tok;
    assign best[0] = {BW{1'b0}};

    // traceback read address (broadcast to all PE pointer banks)
    wire signed [IW:0] tb_i, tb_j;
    wire [IW-1:0] rd_row = tb_i[IW-1:0];
    wire [IW-1:0] rd_col = tb_j[IW-1:0];
    wire [IW-1:0] rd_pe  = rd_row % NUM_PE;
    wire [QW-1:0] rd_q   = rd_row / NUM_PE;
    wire [QW+IW-1:0] rd_addr = rd_q * L + rd_col;       // pointer-bank address = pass*L + col

    wire clr = (state_r == S_DONE);

    genvar p;
    generate
        for (p = 0; p < NUM_PE; p = p + 1) begin : PE_array
            // S for this PE: arrives on data_s during pass 0, from seq_s afterwards
            wire [QW-1:0] tq = tok[p][TW-6-IW -: QW];
            wire [1:0] s_load;
            if (Q == 1) begin : g_s1
                assign s_load = data_s;
            end else begin : g_sq
                assign s_load = (tq == 0) ? data_s : seq_s[tq * NUM_PE + p];
            end

            SW_PE #(.L(L), .NUM_PE(NUM_PE), .PE_ID(p), .DW(DW), .IW(IW), .QW(QW), .Q(Q)) u_pe (
                .clk(clk), .reset(reset), .clr(clr),
                .tok_in(tok[p]), .s_load(s_load), .best_in(best[p]),
                .tok_out(tok[p+1]), .best_out(best[p+1]),
                .rd_addr(rd_addr), .rd_ptr(rd_ptr[p])
            );
        end
    endgenerate

    // ------------------------------------------------------------------
    // recirculation buffer: bottom row of pass q -> PE0 of pass q+1
    // ------------------------------------------------------------------
    wire [TW-1:0] last_tok = tok[NUM_PE];
    wire          last_v   = last_tok[TW-1];
    wire          last_l   = last_tok[TW-3];
    wire [QW-1:0] last_q   = last_tok[TW-6-IW -: QW];

    generate
        if (Q == 1) begin : g_norecirc
            assign recirc_tok = {TW{1'b0}};
        end else begin : g_recirc
            localparam D = L - NUM_PE;
            reg [TW-1:0] dly [0:D-1];
            wire [TW-1:0] push = (last_v && last_q != Q - 1)
                ? {last_tok[TW-1:TW-5], last_tok[TW-6 -: IW], last_q + 1'b1, last_tok[2*DW-1:0]}
                : {TW{1'b0}};
            integer d;
            always @(posedge clk or posedge reset) begin
                if (reset) begin
                    for (d = 0; d < D; d = d + 1) dly[d] <= {TW{1'b0}};
                end else begin
                    dly[0] <= push;
                    for (d = 1; d < D; d = d + 1) dly[d] <= dly[d-1];
                end
            end
            assign recirc_tok = dly[D-1];
        end
    endgenerate

    wire score_done = last_v && last_l && (last_q == Q - 1);

    // ------------------------------------------------------------------
    // traceback
    // ------------------------------------------------------------------
    wire [BW-1:0] fin_best = best[NUM_PE];
    reg  [IW-1:0] bi_r, bj_r;
    reg           tb_start;
    wire          tb_done;
    wire [IW:0]   tb_si, tb_sj;
    wire [IW+1:0] tb_len;

    SW_traceback #(.IW(IW)) u_tb (
        .clk(clk), .reset(reset),
        .start(tb_start), .best_i(bi_r), .best_j(bj_r),
        .cur_i(tb_i), .cur_j(tb_j),
        .ptr(rd_ptr[rd_pe]), .inband(1'b1),
        .s_base(seq_s[rd_row]), .t_base(seq_t[rd_col]),
        .op_valid(aln_valid), .op(aln_op),
        .done(tb_done), .start_i(tb_si), .start_j(tb_sj), .len(tb_len)
    );

    always @(posedge clk or posedge reset) begin
        if (reset) begin
            state_r  <= S_IDLE;
            finish   <= 1'b0;
            max      <= 0;
            end_i    <= 0;
            end_j    <= 0;
            start_i  <= 0;
            start_j  <= 0;
            aln_len  <= 0;
            bi_r     <= 0;
            bj_r     <= 0;
            tb_start <= 1'b0;
        end else begin
            finish   <= 1'b0;
            tb_start <= 1'b0;
            case (state_r)
                S_IDLE: if (accept) state_r <= S_RUN;
                S_RUN: if (score_done) begin
                    max      <= fin_best[BW-1 -: DW];
                    bi_r     <= fin_best[2*IW-1 -: IW];
                    bj_r     <= fin_best[IW-1:0];
                    tb_start <= 1'b1;
                    state_r  <= S_TB;
                end
                S_TB: if (tb_done) begin
                    end_i    <= bi_r;
                    end_j    <= bj_r;
                    start_i  <= tb_si[IW-1:0];
                    start_j  <= tb_sj[IW-1:0];
                    aln_len  <= tb_len;
                    finish   <= 1'b1;
                    state_r  <= S_DONE;
                end
                S_DONE: state_r <= S_IDLE;
            endcase
        end
    end

endmodule

//============================================================================
// Processing element: one DP cell per cycle for row  q*NUM_PE + PE_ID
//============================================================================
module SW_PE #(
    parameter L = 256, NUM_PE = 256, PE_ID = 0,
    parameter DW = 13, IW = 8, QW = 1, Q = 1,
    // derived widths -- do not override
    parameter TW = 1 + 1 + 1 + 2 + IW + QW + DW + DW,
    parameter BW = DW + 2 * IW
)(
    input                 clk,
    input                 reset,
    input                 clr,
    input  [TW-1:0]       tok_in,
    input  [1:0]          s_load,
    input  [BW-1:0]       best_in,
    output [TW-1:0]       tok_out,
    output [BW-1:0]       best_out,
    input  [QW+IW-1:0]    rd_addr,
    output [3:0]          rd_ptr
);
    localparam signed [DW-1:0] MATCH = 8, MISMATCH = -5, GAP_OPEN = -7, GAP_EXT = -3;
    localparam signed [DW-1:0] NEG = -16;
    localparam [1:0] STOP = 2'd0, DIAG = 2'd1, UP = 2'd2, LEFT = 2'd3;

    // ---- unpack input token
    wire                 v_in  = tok_in[TW-1];
    wire                 first = tok_in[TW-2];
    wire [1:0]           t_in  = tok_in[TW-4 -: 2];
    wire [IW-1:0]        col   = tok_in[TW-6 -: IW];
    wire [QW-1:0]        q     = tok_in[TW-6-IW -: QW];
    wire signed [DW-1:0] h_up  = tok_in[2*DW-1 -: DW];
    wire signed [DW-1:0] f_up  = tok_in[DW-1:0];

    // ---- state
    reg  [1:0]           s_r;
    reg  signed [DW-1:0] h_r, e_r, hd_r;     // H(i,j-1), E(i,j-1), H(i-1,j-1)
    reg  [TW-1:0]        tok_r;
    reg  [BW-1:0]        best_r;
    localparam AW = (Q * L > 1) ? $clog2(Q * L) : 1;
    reg  [3:0]           tbm [0:Q*L-1];      // traceback pointer bank, address = q*L + col
    wire [AW-1:0]        wr_addr = q * L + col;

    // ---- cell computation
    wire [1:0]           s_cur  = first ? s_load : s_r;
    wire signed [DW-1:0] h_left = first ? {DW{1'b0}} : h_r;
    wire signed [DW-1:0] e_left = first ? NEG : e_r;
    wire signed [DW-1:0] h_diag = first ? {DW{1'b0}} : hd_r;

    wire signed [DW-1:0] e_open = h_left + GAP_OPEN;
    wire signed [DW-1:0] e_extv = e_left + GAP_EXT;
    wire                 e_ext  = (e_extv > e_open);
    wire signed [DW-1:0] e_new  = e_ext ? e_extv : e_open;

    wire signed [DW-1:0] f_open = h_up + GAP_OPEN;
    wire signed [DW-1:0] f_extv = f_up + GAP_EXT;
    wire                 f_ext  = (f_extv > f_open);
    wire signed [DW-1:0] f_new  = f_ext ? f_extv : f_open;

    wire signed [DW-1:0] d_val  = h_diag + ((s_cur == t_in) ? MATCH : MISMATCH);

    reg  signed [DW-1:0] h_new;
    reg  [1:0]           hsrc;
    always @(*) begin
        h_new = 0;
        hsrc  = STOP;
        if (d_val > h_new) begin h_new = d_val; hsrc = DIAG; end
        if (f_new > h_new) begin h_new = f_new; hsrc = UP;   end
        if (e_new > h_new) begin h_new = e_new; hsrc = LEFT; end
    end

    // ---- best-so-far: merge(best_r, best_in, this cell)
    wire [IW-1:0] row = q * NUM_PE + PE_ID;
    wire [BW-1:0] cand = v_in ? {h_new, row, col} : {BW{1'b0}};

    function better;   // a strictly better than b ?
        input [BW-1:0] a, b;
        begin
            if ($signed(a[BW-1 -: DW]) != $signed(b[BW-1 -: DW]))
                better = $signed(a[BW-1 -: DW]) > $signed(b[BW-1 -: DW]);
            else
                better = (a[2*IW-1:0] < b[2*IW-1:0]);     // smaller i, then smaller j
        end
    endfunction

    wire [BW-1:0] m1 = better(best_in, best_r) ? best_in : best_r;
    wire [BW-1:0] m2 = (v_in && better(cand, m1)) ? cand : m1;

    always @(posedge clk or posedge reset) begin
        if (reset) begin
            s_r    <= 2'd0;
            h_r    <= 0;
            e_r    <= 0;
            hd_r   <= 0;
            tok_r  <= {TW{1'b0}};
            best_r <= {BW{1'b0}};
        end else begin
            tok_r  <= v_in ? {tok_in[TW-1:2*DW], h_new, f_new} : {TW{1'b0}};
            best_r <= clr ? {BW{1'b0}} : m2;
            if (v_in) begin
                s_r  <= s_cur;
                h_r  <= h_new;
                e_r  <= e_new;
                hd_r <= h_up;
            end
        end
    end

    always @(posedge clk) begin
        if (v_in) tbm[wr_addr] <= {hsrc, e_ext, f_ext};
    end

    assign tok_out  = tok_r;
    assign best_out = best_r;
    assign rd_ptr   = tbm[rd_addr[AW-1:0]];
endmodule

//============================================================================
// Traceback FSM (one op per cycle)
//   ptr = {hsrc[1:0], e_ext, f_ext} of cell (cur_i, cur_j), read combinationally
//============================================================================
module SW_traceback #(
    parameter IW = 8
)(
    input                    clk,
    input                    reset,
    input                    start,
    input      [IW-1:0]      best_i,
    input      [IW-1:0]      best_j,
    output reg signed [IW:0] cur_i,
    output reg signed [IW:0] cur_j,
    input      [3:0]         ptr,
    input                    inband,
    input      [1:0]         s_base,
    input      [1:0]         t_base,
    output reg               op_valid,
    output reg [1:0]         op,
    output reg               done,
    output     [IW:0]        start_i,
    output     [IW:0]        start_j,
    output reg [IW+1:0]      len
);
    localparam [1:0] ST_H = 2'd0, ST_E = 2'd1, ST_F = 2'd2;
    localparam [1:0] STOP = 2'd0, DIAG = 2'd1, UP = 2'd2, LEFT = 2'd3;
    localparam [1:0] OP_M = 2'd0, OP_X = 2'd1, OP_I = 2'd2, OP_D = 2'd3;

    reg       busy;
    reg [1:0] st;

    wire [1:0] hsrc  = ptr[3:2];
    wire       e_ext = ptr[1];
    wire       f_ext = ptr[0];
    wire       oob   = cur_i[IW] || cur_j[IW];            // index == -1

    // effective move for this cycle
    reg [1:0] act;     // ST_H means DIAG here
    reg       stop;
    always @(*) begin
        stop = 1'b0;
        act  = ST_H;
        case (st)
            ST_E: act = ST_E;
            ST_F: act = ST_F;
            default: begin
                if (oob || !inband || hsrc == STOP) stop = 1'b1;
                else if (hsrc == DIAG)              act  = ST_H;
                else if (hsrc == UP)                act  = ST_F;
                else                                act  = ST_E;
            end
        endcase
    end

    assign start_i = cur_i + 1'b1;
    assign start_j = cur_j + 1'b1;

    always @(posedge clk or posedge reset) begin
        if (reset) begin
            busy     <= 1'b0;
            st       <= ST_H;
            cur_i    <= 0;
            cur_j    <= 0;
            op_valid <= 1'b0;
            op       <= 2'd0;
            done     <= 1'b0;
            len      <= 0;
        end else begin
            op_valid <= 1'b0;
            done     <= 1'b0;
            if (start) begin
                busy  <= 1'b1;
                st    <= ST_H;
                cur_i <= {1'b0, best_i};
                cur_j <= {1'b0, best_j};
                len   <= 0;
            end else if (busy) begin
                if (stop) begin
                    busy <= 1'b0;
                    done <= 1'b1;
                end else begin
                    op_valid <= 1'b1;
                    len      <= len + 1'b1;
                    case (act)
                        ST_H: begin
                            op    <= (s_base == t_base) ? OP_M : OP_X;
                            cur_i <= cur_i - 1'b1;
                            cur_j <= cur_j - 1'b1;
                            st    <= ST_H;
                        end
                        ST_F: begin
                            op    <= OP_D;
                            cur_i <= cur_i - 1'b1;
                            st    <= f_ext ? ST_F : ST_H;
                        end
                        default: begin
                            op    <= OP_I;
                            cur_j <= cur_j - 1'b1;
                            st    <= e_ext ? ST_E : ST_H;
                        end
                    endcase
                end
            end
        end
    end
endmodule
