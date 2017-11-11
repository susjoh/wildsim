runtime <- "06:00:00"

writeLines(paste0("
#!/bin/sh
#$ -cwd
#$ -l h_rt=", 06:00:00
#$ -V
#$ -l h_vmem=128G

. /etc/profile.d/modules.sh

module load R
R CMD BATCH Scenario_3_full.R




