module SW(
    input clk,
    input reset,
    input [1:0] data_s,
    input [1:0] data_t,
    input valid,
    output reg finish,
    output reg [11:0] max
);  
    parameter RUN_1 = 0;
    parameter RUN_2 = 1;
    
    reg state_r,state_w;
    reg [8:0] counter_r,counter_w;

    //-------------------//
    wire [1:0]  T_out[0:255];
    wire [11:0] Max_out[0:255];
    wire [11:0] V_out[0:255];
    wire [11:0] F_out[0:255];
    
    genvar j;
    generate
        PE PE0(clk,reset,(counter_r==0 && valid),(counter_r>=9'd256),data_s,data_t,12'b0,12'b0,12'b0,
                             T_out[0],Max_out[0],V_out[0],F_out[0]);
        for(j=1;j<256;j=j+1) begin: PE_array
            PE PE1(clk,reset,(counter_r==j),(counter_r>=(9'd256 + j)),data_s,T_out[j-1],Max_out[j-1],V_out[j-1],F_out[j-1],
                             T_out[j],Max_out[j],V_out[j],F_out[j]);
        end
    endgenerate
    //-------------------//

    always@ (*) begin
        if(state_r==RUN_1) 
            state_w = (&counter_r[7:0]) ? RUN_2 : RUN_1;
        else
            state_w = (&counter_r) ? RUN_1 : RUN_2;
    end

    always@ (*) begin
        if(state_r==RUN_1) begin
            if(!valid)
                counter_w = 0;
            else
                counter_w = counter_r + 1'b1;
        end
        else begin
            counter_w = counter_r + 1'b1;
        end
    end

    always@ (posedge clk or posedge reset) begin
        if(reset) begin
            state_r <= RUN_1;
            counter_r <= 0;
            finish <= 0;
            max <= 0;
        end
        else begin
            state_r <= state_w;
            counter_r <= counter_w;

            if(&counter_w) begin
                finish <= 1;
                max <= Max_out[255];
            end
            else begin
                finish <= 0;
                max <= 0;
            end
        end
    end

endmodule

module PE(
    input clk,
    input reset,
    input valid,
    input stop,
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
        score = (S_out_w==T_in) ? MATCH : MISMATCH;
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
            if(!stop) begin
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
            else begin
                Max_out_w = Max_out_r;
                V_out_w = V_out_r;
                E_out_w = E_out_r;
                F_out_w = F_out_r;
            end
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
