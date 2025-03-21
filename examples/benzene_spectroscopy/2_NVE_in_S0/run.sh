#!/bin/bash

for n in {0..99}; do
    echo $n
    cd TRAJ/traj-${n}/
    cp ../../NAMD.in .
    cp ../../submit.SQD .
    # sbatch submit.SQD
    cd ../../
done