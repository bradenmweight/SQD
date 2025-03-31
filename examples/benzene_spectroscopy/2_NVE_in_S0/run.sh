#!/bin/bash

#for n in {0..19}; do
for n in {20..99}; do
    echo $n
    cd TRAJ/traj-${n}/
    rm -r MD.out output.slurm MD_OUTPUT
    cp ../../NAMD.in .
    cp ../../submit.SQD .
    sbatch submit.SQD
    cd ../../
    sleep 0.05
done