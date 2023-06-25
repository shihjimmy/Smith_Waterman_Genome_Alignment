read_file -format verilog  "SW.v"
current_design [get_designs SW]

#You may modified the clock constraints 
#or add more constraints for your design
####################################################
set cycle  10
set_max_area 0       
####################################################

#The following are design spec. for synthesis
#You can NOT modify this seciton 
#####################################################
create_clock -name clk -period $cycle [get_ports clk]
set_fix_hold                          [get_clocks clk]
set_dont_touch_network                [get_clocks clk]
set_ideal_network                     [get_ports clk]
set_clock_uncertainty            0.1  [get_clocks clk] 
set_clock_latency                0.5  [get_clocks clk] 

set_max_fanout 6 [all_inputs] 

set_operating_conditions -min_library fast -min fast -max_library slow -max slow
set_wire_load_model -name tsmc13_wl10 -library slow  
set_drive        1     [all_inputs]
set_load         1     [all_outputs]
set t_in   0.1
set t_out  0.1
set_input_delay  $t_in  -clock clk [remove_from_collection [all_inputs] [get_ports clk]]
set_output_delay $t_out -clock clk [all_outputs]
#####################################################


#Compile and save files
#You may modified setting of compile 
#####################################################
compile
write_sdf -version 2.1 SW_syn.sdf
write -format verilog -hier -output SW_syn.v
write -format ddc     -hier -output SW_syn.ddc  
#####################################################  