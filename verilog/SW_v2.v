module SW(
    input clk,
    input reset,
    input [1:0] data_s,
    input [1:0] data_t,
    input valid,
    output reg finish,
    output reg [11:0] max
);  
    parameter IDLE = 0;
    parameter CALC_1 = 1;
    parameter BUFF_1 = 2;
    parameter BUFF_2 = 3;
    parameter CALC_2 = 4;
    parameter CALC_3 = 5;
    
    integer i;
    reg [2:0]  state_r,state_w;
    reg [9:0]  count_r, count_w;

    //-------------------//
    wire [1:0]   T_out[0:127];
    wire [11:0]  Max_out[0:127];
    wire [11:0]  V_out[0:127];
    wire [11:0]  F_out[0:127];
    reg [1:0]    S_in[0:127];
    reg [1:0]    T_in[0:127];
    reg [11:0]   Max_in;
    reg [11:0]   V_in;
    reg [11:0]   F_in;
    reg [127:0]  PE_valid;
    wire [1:0]   data_s_out,data_t_out;
    wire [11:0]  Max,V,F;

    genvar j;
    generate
        PE PE0(clk,reset,(PE_valid[0]),S_in[0],T_in[0],Max_in,V_in,F_in,
                             T_out[0],Max_out[0],V_out[0],F_out[0]);
                        
        for(j=1;j<128;j=j+1) begin: PE_array
            PE PE1(clk,reset,(PE_valid[j]),S_in[j],T_in[j],Max_out[j-1],V_out[j-1],F_out[j-1],
                             T_out[j],Max_out[j],V_out[j],F_out[j]);
        end
    endgenerate

    buffer buf0(clk, reset, (count_r>=10'd128), data_s, T_out[127], Max_out[127], V_out[127], F_out[127],
                data_s_out, data_t_out, Max, V, F);
    //-------------------//

    always@ (*) begin
        case (state_r)
            IDLE:   state_w = (valid) ? CALC_1 : IDLE;
            CALC_1: state_w = (&count_r[6:0]) ? BUFF_1 : CALC_1;    //count_r = 127 (0-127)
            BUFF_1: state_w = (&count_r[7:0]) ? BUFF_2 : BUFF_1;    //count_r = 255 (128-255)
            BUFF_2: state_w = (count_r==10'd383) ? CALC_2 : BUFF_2; //count_r = 383 (256-383)
            CALC_2: state_w = (&count_r[8:0]) ? CALC_3 : CALC_2;    //count_r = 511 (384-511)
            CALC_3: state_w = (count_r==10'd638) ? IDLE : CALC_3;   //count_r = 638 (512-638)
            default : state_w = IDLE;
        endcase
    end

    // counter
    always@ (*) begin
        if(state_r!=IDLE || valid)
            count_w = count_r + 1'b1;
        else
            count_w = count_r;
    end

    // S_in and T_in signals
    always@ (*) begin
        for(i=0;i<128;i=i+1) begin
            S_in[i] = 0;
            T_in[i] = 0;
        end

        case(state_r) 
            IDLE: begin
                S_in[0] = (valid) ? data_s : 0;
                T_in[0] = (valid) ? data_t : 0;
            end
            CALC_1: begin
                T_in[0] = data_t;

                for(i=1;i<128;i=i+1) begin
                    S_in[i] = (PE_valid[i]) ? data_s : 0;
                    T_in[i] = T_out[i-1];
                end
            end
            BUFF_1: begin
                T_in[0] = data_t;

                for(i=1;i<128;i=i+1) begin
                    T_in[i] = T_out[i-1];
                end
            end
            BUFF_2: begin
                T_in[0] = data_t_out;

                for(i=1;i<128;i=i+1) begin
                    S_in[i] = (PE_valid[i]) ? data_s_out : 0;
                    T_in[i] = T_out[i-1];
                end
            end
            CALC_2: begin
                T_in[0] = data_t_out;

                for(i=1;i<128;i=i+1) begin
                    T_in[i] = T_out[i-1];
                end
            end
            CALC_3: begin
                T_in[0] = data_t_out;

                for(i=1;i<128;i=i+1) begin
                    T_in[i] = T_out[i-1];
                end
            end
        endcase
    end

    // PE-valid signals
    always@ (*) begin
        for(i=0;i<128;i=i+1) begin
            PE_valid[i] = 0;
        end

        case(state_r)
            IDLE: begin
                PE_valid[0] = (valid) ? 1 : 0;
            end
            CALC_1: begin
                for(i=1;i<128;i=i+1) begin
                    PE_valid[i] = (count_r==i);
                end
            end
            BUFF_2: begin
                for(i=0;i<128;i=i+1) begin
                    PE_valid[i] = (count_r==(i+10'd256));
                end
            end
        endcase
    end

    //Max_in V_in F_in signals
    always@ (*) begin
        Max_in = 0;
        V_in   = 0;
        F_in   = 0;

        if(state_r == BUFF_2 || state_r == CALC_2 || state_r == CALC_3) begin
            Max_in = Max;
            V_in   = V;
            F_in   = F;
        end
    end

    always@(posedge clk or posedge reset) begin
        if(reset) begin
            state_r <= IDLE;
            count_r <= 0;
            finish  <= 0;
            max     <= 0;
        end
        else begin
            state_r <= state_w;
            count_r <= count_w;
            if(state_r==IDLE && count_r==10'd639) begin
                finish  <= 1;
                max     <= Max_out[127]; 
            end
            else begin
                finish  <= 0;
                max     <= 0;
            end
        end
    end

endmodule

module buffer(
    input clk,
    input reset,
    input valid,
    input [1:0]  data_s,
    input [1:0]  data_t,
    input [11:0] Max_in,
    input [11:0] V_in,
    input [11:0] F_in,
    output [1:0] data_s_out,
    output [1:0] data_t_out,
    output [11:0] Max_out,
    output [11:0] V_out,
    output [11:0] F_out
);
    integer i;
    reg [1:0]  data_s_out_r[0:127], data_s_out_w[0:127];
    reg [1:0]  data_t_out_r[0:127], data_t_out_w[0:127];
    reg [11:0] Max_out_r[0:127], Max_out_w[0:127];
    reg [11:0] V_out_r[0:127],   V_out_w[0:127];
    reg [11:0] F_out_r[0:127],   F_out_w[0:127];

    assign data_s_out = data_s_out_r[127];
    assign data_t_out = data_t_out_r[127];
    assign Max_out = Max_out_r[127];
    assign V_out = V_out_r[127];
    assign F_out = F_out_r[127];

    always@ (*) begin
        for(i=0;i<128;i=i+1) begin
            Max_out_w[i] = Max_out_r[i];
            V_out_w[i] = V_out_r[i];
            F_out_w[i] = F_out_r[i];    
            data_s_out_w[i] = data_s_out_r[i];
            data_t_out_w[i] = data_t_out_r[i];
        end

        if(valid) begin
            Max_out_w[0] = Max_in;
            V_out_w[0] = V_in;
            F_out_w[0] = F_in;
            data_s_out_w[0] = data_s;
            data_t_out_w[0] = data_t;

            for(i=1;i<128;i=i+1) begin
                Max_out_w[i] = Max_out_r[i-1];
                V_out_w[i] = V_out_r[i-1];
                F_out_w[i] = F_out_r[i-1];
                data_s_out_w[i] = data_s_out_r[i-1];
                data_t_out_w[i] = data_t_out_r[i-1];
            end
        end
    end

    always@(posedge clk or posedge reset) begin
        if(reset) begin
            for(i=0;i<128;i=i+1) begin
                Max_out_r[i] <= 0;
                V_out_r[i] <= 0;
                F_out_r[i] <= 0;
                data_s_out_r[i] <= 0;
                data_t_out_r[i] <= 0;
            end
        end
        else begin
            for(i=0;i<128;i=i+1) begin
                Max_out_r[i] <= Max_out_w[i];
                V_out_r[i] <= V_out_w[i];
                F_out_r[i] <= F_out_w[i];
                data_s_out_r[i] <= data_s_out_w[i];
                data_t_out_r[i] <= data_t_out_w[i];
            end
        end
    end

endmodule

module PE(
    input clk,
    input reset,
    input valid,
    input [1:0]  S_in,
    input [1:0]  T_in,
    input [11:0] Max_in,
    input [11:0] V_in,
    input [11:0] F_in,
    output [1:0] T_out,
    output [11:0] Max_out,
    output [11:0] V_out,
    output [11:0] F_out
);
    localparam GAP_OPEN = -12'd7; 
    localparam GAP_EXTEND =  -12'd3; 
    localparam MATCH = 12'd8;
    localparam MISMATCH = -12'd5;

    parameter IDLE = 0;
    parameter RUN = 1;

    // For DFF parameters
    reg state_r ,state_w;
    reg [1:0]  S_out_r,S_out_w;
    reg [1:0]  T_out_r;
    reg [11:0] V_diag_r;
    reg [11:0] E_out_r,E_out_w;
    reg [11:0] Max_out_r,Max_out_w;
    reg [11:0] V_out_r,V_out_w;
    reg [11:0] F_out_r,F_out_w;

    assign Max_out = Max_out_r;
    assign F_out = F_out_r;
    assign V_out = V_out_r;
    assign T_out = T_out_r;

    // wire and regs
    reg [11:0] score;
    reg [11:0] max_out_temp;
    reg [11:0] max_E_1,max_E_2;
    reg [11:0] max_F_1,max_F_2;
    reg [11:0] max_V_1,max_V_2;

    always@ (*) begin
        if(state_r==IDLE)
            state_w = (valid) ? RUN : IDLE;
        else 
            state_w = state_r;
    end

    always@ (*) begin
        S_out_w = (valid) ? S_in : S_out_r;
        score   = (S_out_w==T_in) ? MATCH : MISMATCH;
        max_V_1 = $signed(V_diag_r) + $signed(score);
        max_F_1 = $signed(V_in) + $signed(GAP_OPEN);
        max_F_2 = $signed(F_in) + $signed(GAP_EXTEND);
        max_E_1 = $signed(E_out_r) + $signed(GAP_EXTEND);
        max_E_2 = $signed(V_out_r) + $signed(GAP_OPEN);

        if($signed(F_out_w) >= $signed(E_out_w))
            max_V_2 = F_out_w;
        else 
            max_V_2 = E_out_w;

        if($signed(V_out_w) >= $signed(Max_in))
            max_out_temp = V_out_w;
        else 
            max_out_temp = Max_in;
    end

    always@ (*) begin
        if(state_r==IDLE && (!valid)) begin
            Max_out_w = 0;
            V_out_w = 0;
            E_out_w = 0;
            F_out_w = 0;
        end
        else begin
            if($signed(max_out_temp) >= $signed(Max_out_r))
                Max_out_w = max_out_temp;
            else 
                Max_out_w = Max_out_r;

            if($signed(max_V_1) >= $signed(max_V_2))
                V_out_w = max_V_1;
            else 
                V_out_w = max_V_2;

            if($signed(max_E_1) >= $signed(max_E_2))
                E_out_w = max_E_1;
            else 
                E_out_w = max_E_2;
            
            if($signed(max_F_1) >= $signed(max_F_2))
                F_out_w = max_F_1;
            else 
                F_out_w = max_F_2;
        end
    end

    always@ (posedge clk or posedge reset) begin
        if(reset) begin
            state_r   <= 0;
            S_out_r   <= 0;
            T_out_r   <= 0;
            V_diag_r  <= 0;
            E_out_r   <= 0;
            Max_out_r <= 0;
            V_out_r   <= 0;
            F_out_r   <= 0;
        end
        else begin
            state_r   <= state_w;
            S_out_r   <= S_out_w;
            T_out_r   <= T_in;
            V_diag_r  <= V_in;
            E_out_r   <= E_out_w;
            Max_out_r <= Max_out_w;
            V_out_r   <= V_out_w;
            F_out_r   <= F_out_w;
        end
    end

endmodule
