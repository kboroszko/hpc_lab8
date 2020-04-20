txt="#!/bin/bash -l
#SBATCH --job-name many-${1}-${2}-job          # this will be shown in the queueing system
#SBATCH --output \"${1}-${2}.out\"   # stdout redirection
#SBATCH --error \"${1}-${2}.txt\"    # stderr redirection
#SBATCH --account \"GC80-33\"           # the number of our grant
#SBATCH --nodes 1                     # how many nodes we want
#SBATCH --tasks-per-node ${1}           # each node is 2 socket, 12 core, so we want 24 tasks on each node
#SBATCH --time 00:03:00               # if the job runs longer than this, itll be killed

srun floyd-warshall-par.exe   ${2}          # what command to run"

echo "${txt}"
