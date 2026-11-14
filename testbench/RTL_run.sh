# VCS (run inside testbench/ after `source license.sh`)
vcs tb_SW.v +define+RTL       -full64 -R -debug_access+all +v2k +notimingcheck    # NUM_PE = 256
vcs tb_SW.v +define+RTL+SW_v2 -full64 -R -debug_access+all +v2k +notimingcheck    # NUM_PE = 128
# add +define+FSDB to dump SW.fsdb for Verdi
