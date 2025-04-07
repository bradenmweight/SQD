#!/bin/bash

#for n in {0..19}; do
#for n in {20..99}; do
#for n in {100..103}; do
#for n in {104..149}; do
#for n in {150..199}; do
#for n in {200..205}; do
#for n in {206..220}; do
#for n in {221..230}; do
#for n in {231..250}; do
for n in {251..260}; do
    echo $n
    cd TRAJ/traj-${n}/
    rm -r MD.out output.slurm MD_OUTPUT
    cp ../../NAMD.in .
    cp ../../submit.SQD .
    sbatch submit.SQD
    cd ../../
    sleep 0.05
done