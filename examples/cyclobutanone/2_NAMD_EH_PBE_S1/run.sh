#!/bin/bash

for n in {0..99}; do
    echo $n
    cd TRAJ/traj-${n}/
    rm -r MD* output.slurm
    cp ../../NAMD.in .
    cp ../../submit.SQD .
    sbatch submit.SQD
    cd ../../
    sleep 0.1
done