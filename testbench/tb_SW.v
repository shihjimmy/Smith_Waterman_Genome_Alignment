`timescale 1ns/10ps
`define clk_period 12              
`define END_CYCLE 10000 

`ifdef SW
    `define PAT_S "s1.dat"
    `define PAT_T "t1.dat"
    `define ANS "golden1.dat"
    `define SDFFILE "../verilog/SW_syn.sdf"
    `ifdef SDF
        `include "../verilog/SW_syn.v"
        `include "../verilog/tsmc13.v"
    `endif
    `ifdef RTL
        `include "../verilog/SW.v"
    `endif
`endif

`ifdef SW_v2
    `define PAT_S "s2.dat"
    `define PAT_T "t2.dat"
    `define ANS "golden2.dat"
    `define SDFFILE "../verilog/SW_v2_syn.sdf"
    `ifdef SDF
        `include "../verilog/SW_v2_syn.v"
        `include "../verilog/tsmc13.v"
    `endif
    `ifdef RTL
        `include "../verilog/SW_v2.v"
    `endif
`endif

module tb_SW;

    reg clk;
    reg reset;
    reg [1:0] seq_s [0:255];
    reg [1:0] seq_t [0:255];
    reg [11:0] max_num[0:1];
    reg valid;
    reg [1:0] data_s;
    reg [1:0] data_t;
    wire finish;
    wire [11:0] max;
    
    initial begin
        $readmemh(`PAT_S, seq_s);
        $readmemh(`PAT_T, seq_t);
        $readmemh(`ANS, max_num);
    end

    SW SW0(.clk(clk), .reset(reset), .valid(valid), .data_s(data_s), .data_t(data_t),
        .finish(finish), .max(max) ); 

    initial begin
        `ifdef SDF
            $sdf_annotate(`SDFFILE, SW0);
        `endif

        $fsdbDumpfile("SW.fsdb");
        $fsdbDumpvars(0, tb_SW, "+mda");
    end

    always begin 
        #(`clk_period/2) clk = ~clk; 
    end 
    
    initial begin
        reset = 1'b0;
        clk = 1'b0;
        #(`clk_period)   reset = 1'b1;
        #(`clk_period*2) reset = 1'b0;
    end

    integer i;
    initial begin
        valid = 1'b0;
        i = 0;
        #(`clk_period*5);

        @(negedge clk) begin
            valid = 1'b1;
            data_s = seq_s[i];
            data_t = seq_t[i];
            for (i=1;i<256;i=i+1) begin
                @(negedge clk) begin 
                    data_s = seq_s[i];
                    data_t = seq_t[i];
                end
            end
        end

        @(negedge clk) begin
            valid = 1'b0;
            data_s = 2'b0;
            data_t = 2'b0;
        end
    end

    always@ (negedge clk) begin
        if(finish == 1'b1) begin
            if (max_num[0] == max) begin    
                $display("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                $display("---------- The test result is ..... PASS -------------------\n");
                $display("Golden: %d\n",max_num[0]);
                $display("Yours: %d\n",max);
                $display("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                $finish;
            end 
            else begin
                $display("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                $display("---------- The test result is ..... FAIL -------------------\n");
                $display("Golden: %d\n",max_num[0]);
                $display("Yours: %d\n",max);
                $display("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
                $finish;
            end
        end
    end

    initial begin 
        #(`clk_period*`END_CYCLE);
        $display("Time out!!"); 
        $finish;
    end

endmodule
