#!/bin/bash
git pull
make
sbatch graph.batch
squeue -u kb385496
