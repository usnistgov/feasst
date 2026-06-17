#!/bin/bash
write="trials_per_write=1e4 output_file=lj"
feasst << EOF
Restart checkpoint_file=checkpoint.fst
Log $write.csv
Movie $write.xyz
Energy ${write}_en.csv
Metropolis trials_per_cycle=1e4 cycles_to_complete=1e2
Run until=complete
EOF
