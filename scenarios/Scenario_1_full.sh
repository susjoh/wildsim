#!/bin/sh
#$ -cwd
#$ -l h_rt=06:00:00
#$ -V
#$ -pe sharedmem 8
#$ -l h_vmem=16G

. /etc/profile.d/modules.sh

module load R
R CMD BATCH Scenario_1_full.R




