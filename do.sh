#!/bin/bash
git pull
cc stats.c -o stats.exe
sbatch stats.batch
rm *.err *.out
squeue -u kb385496