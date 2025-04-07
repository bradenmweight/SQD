#!/bin/bash

for n in {0..5}; do
    echo $n
    cd TRAJ/traj-${n}/
    rm -r MD* output.slurm
    rm -r EL_STR*
    cp ../../NAMD.in .
    cp ../../submit.SQD .
    sbatch submit.SQD
    cd ../../
    sleep 0.1
done