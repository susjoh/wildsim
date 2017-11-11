#!/bin/sh
#$ -cwd
#$ -l h_rt=06:00:00
#$ -V
#$ -l h_vmem=8G

. /etc/profile.d/modules.sh

module load R
R CMD BATCH Scenario_4_full.R




