#
#
# Parse Simulation Output
# Single Population
# Susan Johnston, Lewis Spurgin
# October 2017
#
#
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Load libraries, functions and data information   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(GenABEL)
library(plyr)
library(ggplot2)
library(wildsim)


model.name <- "Scenario_0_10k_1000"
markerfile <- "p1_mrk_001.txt"
phenofile <- "p1_data_001.txt"

marker.prefix <- gsub(".txt", "", markerfile)
qtl.prefix <- gsub("mrk", "qtl", marker.prefix)

flanking.window <- 10

report <- readLines("report.txt")
params <- readLines(paste0(model.name, ".prm"))

qtl.positions <- read.table(paste0(gsub("p1_mrk", "lm_qtl", marker.prefix), ".txt"), header = T)
qtl.effects   <- read.table(paste0(gsub("p1_mrk", "p1_freq_qtl", marker.prefix), ".txt"), header = T, fill = T)
names(qtl.effects) <- c("ID", "Gen", "Chr", "Var", "Allele.Freq.1", "Allele.Freq.2")
qtl.positions$BP <- qtl.positions$Position * 1e6
qtl.positions <- join(qtl.positions, qtl.effects)

qtl.positions

rm(qtl.effects)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Prepare and convert files                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Format phenotype file and save

phenotab <- read.table(phenofile, header = T)
head(phenotab)

phenotab$id <- phenotab$Progeny
phenotab$sex <- ifelse(phenotab$Sex == "M", 1, 0)
phenotab$PlinkSex <- ifelse(phenotab$Sex == "M", 1, 2)

write.table(phenotab, paste0(marker.prefix, ".pheno"), row.names = F, sep = "\t", quote = F)

#~~ Create files to add sex and phenotypic informtion to the PLINK Files

pheno2 <- cbind(Family = 1, subset(phenotab, select = c(Progeny, Phen)))
head(pheno2)

write.table(pheno2, paste0(marker.prefix, ".phenochange.txt"), row.names = F, col.names = F, quote = F)

rm(pheno2)

pheno2 <- cbind(Family = 1, subset(phenotab, select = c(Progeny, PlinkSex)))
head(pheno2)

write.table(pheno2, paste0(marker.prefix, ".sexchange.txt"), row.names = F, col.names = F, quote = F)

rm(pheno2)

#~~ Format the marker file. Need to scrub first line, add family, and replace multiple spaces with a single space.

system("cmd", input = paste0("sed \"1d\" ", markerfile,
                             " | sed -e \"s/  */ /g\"  | sed -e \"s/^/1 /g\" | sed \"s/ / 0 0 0 /2\" > ",
                             marker.prefix, ".ped"))

#~~ Format qtlfile. Need to scrub first line, add family, and replace multiple spaces with a single space.


system("cmd", input = paste0("sed \"1d\" ", gsub("mrk", "qtl", markerfile),
                             " | sed -e \"s/  */ /g\"  | sed -e \"s/^/1 /g\" | sed \"s/ / 0 0 0 /2\" > ",
                             qtl.prefix, ".ped"))


#~~ Create map files - convert linkage map positions to approx bp positions (cM * 1e6)

# markers
maptab <- read.table(paste0(gsub("p1_m", "lm_m", marker.prefix), ".txt"), header = T)
head(maptab)
maptab <- maptab[,c("Chr", "ID", "Position")]
maptab$bp <- maptab$Position*1e6

write.table(maptab, paste0(marker.prefix, ".map"), row.names = F, col.names = F, quote = F)

# qtl
maptab <- read.table(paste0(gsub("p1_mrk", "lm_qtl", marker.prefix), ".txt"), header = T)
head(maptab)
maptab <- maptab[,c("Chr", "ID", "Position")]
maptab$bp <- maptab$Position*1e6

write.table(maptab, paste0(qtl.prefix, ".map"), row.names = F, col.names = F, quote = F)

#~~ Add sex to PLINK

RunPLINK(paste0("--file ", marker.prefix, " --no-pheno --update-sex ", marker.prefix, ".sexchange.txt --recode --out ", marker.prefix))

RunPLINK(paste0("--file ", qtl.prefix   , " --no-pheno --update-sex ", marker.prefix, ".sexchange.txt --recode --out ", qtl.prefix))


#~~ Add phenotypes to PLINK

RunPLINK(paste0("--file ", marker.prefix, " --pheno ", marker.prefix, ".phenochange.txt --recode --out ", marker.prefix))

RunPLINK(paste0("--file ", qtl.prefix, " --pheno ", marker.prefix, ".phenochange.txt --recode --out ", qtl.prefix))

#~~ Remove monomorphic SNPs from PLINK

RunPLINK(paste0("--file ", marker.prefix, " --maf 0.01 --recode --out ", marker.prefix))

#~~ Create a merged .ped file

RunPLINK(paste0("--file ", marker.prefix, " --merge ", qtl.prefix, ".ped ", qtl.prefix, ".map --recode --out ", marker.prefix, "_merged"))

#~~ Create GCTA objects for unmerged and merged files

RunPLINK(paste0("--file ", marker.prefix, " --make-bed --out ", marker.prefix))
RunGCTA(paste0("--bfile ", marker.prefix, " --make-grm-gz --out ", marker.prefix, "_GRM"))

RunPLINK(paste0("--file ", marker.prefix, "_merged --make-bed --out ", marker.prefix, "_merged"))
RunGCTA(paste0("--bfile ", marker.prefix, "_merged --make-grm-gz --out ", marker.prefix, "_merged_GRM"))

#~~ Convert to GenABEL format for unmerged and merged files

maptab <- read.table(paste0(marker.prefix, ".map"))
maptab <- maptab[,-3]
write.table(maptab, paste0(marker.prefix, ".map"), row.names = F, col.names = F, quote = F)

mergemaptab <- read.table(paste0(marker.prefix, "_merged.map"))
mergemaptab <- mergemaptab[,-3]
write.table(mergemaptab, paste0(marker.prefix, "_merged.map"), row.names = F, col.names = F, quote = F)

convert.snp.ped(pedfile = paste0(marker.prefix, ".ped"),
                mapfile = paste0(marker.prefix, ".map"),
                outfile = paste0(marker.prefix, ".gen"),
                traits = 1,
                mapHasHeaderLine = F)

convert.snp.ped(pedfile = paste0(marker.prefix, "_merged.ped"),
                mapfile = paste0(marker.prefix, "_merged.map"),
                outfile = paste0(marker.prefix, "_merged.gen"),
                traits = 1,
                mapHasHeaderLine = F)

abeldata <- load.gwaa.data(phenofile = paste0(marker.prefix, ".pheno"),
                           genofile = paste0(marker.prefix, ".gen"))

abeldata.merged <- load.gwaa.data(phenofile = paste0(marker.prefix, ".pheno"),
                                  genofile = paste0(marker.prefix, "_merged.gen"))

#~~ Clean up files

if(Sys.info()["sysname"] == "Windows") {

  system("del *.ped")
  system("del *.map")
  system("del *.sexchange.txt")
  system("del *.nosex")

} else {

  system("rm *.ped")
  system("rm *.map")
  system("rm *.sexchange.txt")
  system("rm *.nosex")
}



#~~ Create a directory for results

system("mkdir results")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Characterise phenotype and population structure. #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(phenotab)

# Phenotype distribution

ggplot(phenotab, aes(Phen)) + geom_histogram(bins = 20)

# QTL distribution

ggplot(qtl.positions, aes(Var)) + geom_histogram()


# Population structure

data1.gkin <- ibs(abeldata[, autosomal(abeldata)], weight="freq")
data1.dist <- as.dist(0.5 - data1.gkin)
data1.mds <- cmdscale(data1.dist)

png("results/Population_Structure.png", width = 4, height = 4, units = "in", res = 300)
plot(data1.mds)
dev.off()

rm(data1.mds, data1.dist)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Run the GWAS                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Model with marker data

model1 <- egscore(Phen, abeldata, kin = data1.gkin)
lambda(model1)
model1.results <- results(model1)
model1.results$SNP.Name <- row.names(model1.results)

png("results/GWAS_Markers_Only.png", width = 8, height = 4, units = "in", res = 300)
FullSummary(model1)
dev.off()

#~~ Model with merged data

model2 <- egscore(Phen, abeldata.merged, kin = data1.gkin)
lambda(model2)
model2.results <- results(model2)
model2.results$SNP.Name <- row.names(model2.results)

png("results/GWAS_Markers_with_QTL.png", width = 8, height = 4, units = "in", res = 300)
FullSummary(model2, qtl.ids = qtl.positions$ID)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Estimate the genomic heritability                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

RunGCTA(paste0("--grm-gz ", marker.prefix, "_GRM --pheno ",
               marker.prefix, ".phenochange.txt --reml --out ",
               marker.prefix, "_GRM_res"))

gcta.res <- read.table(paste0(marker.prefix, "_GRM_res.hsq"), header = T, sep = "\t", fill = T)
gcta.res

RunGCTA(paste0("--grm-gz ", marker.prefix, "_merged_GRM --pheno ",
               marker.prefix, ".phenochange.txt --reml --out ",
               marker.prefix, "_merged_GRM_res"))

gcta.res.merged <- read.table(paste0(marker.prefix, "_merged_GRM_res.hsq"), header = T, sep = "\t", fill = T)
gcta.res.merged


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Estimate variance attributed to QTL              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

qtl.positions
qtl.positions.merged <- qtl.positions

#~~ First, run the regions in the dataset that does not include the QTL.

qtl.positions$QTL.Position <- NA
qtl.positions$Vqtl <- NA
qtl.positions$Vqtl.SE <- NA
qtl.positions$Va <- NA
qtl.positions$Va.SE <- NA
qtl.positions$Vp <- NA
qtl.positions$Vp.SE <- NA
qtl.positions$qtl2 <- NA
qtl.positions$qtl2.SE <- NA
qtl.positions$h2 <- NA
qtl.positions$h2.SE <- NA
qtl.positions$LogL <- NA


for(i in 1:nrow(qtl.positions)){

  focal.snp <- which(model1.results$Chromosome == qtl.positions$Chr[i] &
                       model1.results$Position > qtl.positions$BP[i])[1]

  qtl.positions$QTL.Position[i] <- model1.results$SNP.Name[focal.snp]

  snp.list <- model1.results$SNP.Name[(focal.snp - flanking.window + 1):(focal.snp + flanking.window)]
  write.table(snp.list, paste0(marker.prefix, ".snplist.", i, ".txt"), row.names = F, col.names = F, quote = F)

  RunGCTA(paste0("--bfile ", marker.prefix, " --make-grm-gz --extract ", marker.prefix, ".snplist.", i, ".txt --out ", marker.prefix, ".snplist.", i , "_GRM"))
  RunGCTA(paste0("--bfile ", marker.prefix, " --make-grm-gz --exclude ", marker.prefix, ".snplist.", i, ".txt --out ", marker.prefix, ".snplist.", i , "_wo_GRM"))

  writeLines(paste0(marker.prefix, ".snplist.", i , "_GRM\n",
                    marker.prefix, ".snplist.", i , "_wo_GRM"),
             paste0(marker.prefix, ".snplist.", i , "_GRMs.txt"))

  RunGCTA(paste0("--mgrm-gz ", marker.prefix, ".snplist.", i , "_GRMs.txt --pheno ", marker.prefix, ".phenochange.txt --reml --out ", marker.prefix, "_GRMs_res"))

  gcta.res.temp <- read.table(paste0(marker.prefix, "_GRMs_res.hsq"), header = T, sep = "\t", fill = T)

  qtl.positions$Vqtl[i]    <- gcta.res.temp[1,2]
  qtl.positions$Vqtl.SE[i] <- gcta.res.temp[1,3]
  qtl.positions$Va[i]      <- gcta.res.temp[2,2]
  qtl.positions$Va.SE[i]   <- gcta.res.temp[2,3]
  qtl.positions$Vp[i]      <- gcta.res.temp[4,2]
  qtl.positions$Vp.SE[i]   <- gcta.res.temp[4,3]
  qtl.positions$qtl2[i]    <- gcta.res.temp[5,2]
  qtl.positions$qtl2.SE[i] <- gcta.res.temp[5,3]
  qtl.positions$h2[i]      <- gcta.res.temp[6,2]
  qtl.positions$h2.SE[i]   <- gcta.res.temp[6,3]
  qtl.positions$LogL[i]    <- gcta.res.temp[8,2]

  system(paste0("rm ", marker.prefix, ".snplist.", i , "*"))

}


#~~ Second, run the regions in the dataset including the QTL.

qtl.positions.merged$QTL.Position <- NA
qtl.positions.merged$Vqtl <- NA
qtl.positions.merged$Vqtl.SE <- NA
qtl.positions.merged$Va <- NA
qtl.positions.merged$Va.SE <- NA
qtl.positions.merged$Vp <- NA
qtl.positions.merged$Vp.SE <- NA
qtl.positions.merged$qtl2 <- NA
qtl.positions.merged$qtl2.SE <- NA
qtl.positions.merged$h2 <- NA
qtl.positions.merged$h2.SE <- NA
qtl.positions.merged$LogL <- NA


for(i in 1:nrow(qtl.positions.merged)){

  focal.snp <- which(model2.results$SNP.Name == qtl.positions.merged$ID[i])

  qtl.positions.merged$QTL.Position[i] <- model2.results$SNP.Name[focal.snp]

  snp.list <- model1.results$SNP.Name[(focal.snp + 1 - flanking.window):(focal.snp + flanking.window)]
  write.table(snp.list, paste0(marker.prefix, ".snplist.", i, ".txt"), row.names = F, col.names = F, quote = F)

  RunGCTA(paste0("--bfile ", marker.prefix, "_merged --make-grm-gz --extract ", marker.prefix, ".snplist.", i, ".txt --out ", marker.prefix, ".snplist.", i , "_merged_GRM"))
  RunGCTA(paste0("--bfile ", marker.prefix, "_merged --make-grm-gz --exclude ", marker.prefix, ".snplist.", i, ".txt --out ", marker.prefix, ".snplist.", i , "_merged_wo_GRM"))

  writeLines(paste0(marker.prefix, ".snplist.", i , "_merged_GRM\n",
                    marker.prefix, ".snplist.", i , "_merged_wo_GRM"),
             paste0(marker.prefix, ".snplist.", i , "_merged_GRMs.txt"))

  RunGCTA(paste0("--mgrm-gz ", marker.prefix, ".snplist.", i , "_merged_GRMs.txt --pheno ", marker.prefix, ".phenochange.txt --reml --out ", marker.prefix, "_merged_GRMs_res"))

  gcta.res.temp <- read.table(paste0(marker.prefix, "_merged_GRMs_res.hsq"), header = T, sep = "\t", fill = T)

  qtl.positions.merged$Vqtl[i]    <- gcta.res.temp[1,2]
  qtl.positions.merged$Vqtl.SE[i] <- gcta.res.temp[1,3]
  qtl.positions.merged$Va[i]      <- gcta.res.temp[2,2]
  qtl.positions.merged$Va.SE[i]   <- gcta.res.temp[2,3]
  qtl.positions.merged$Vp[i]      <- gcta.res.temp[4,2]
  qtl.positions.merged$Vp.SE[i]   <- gcta.res.temp[4,3]
  qtl.positions.merged$qtl2[i]    <- gcta.res.temp[5,2]
  qtl.positions.merged$qtl2.SE[i] <- gcta.res.temp[5,3]
  qtl.positions.merged$h2[i]      <- gcta.res.temp[6,2]
  qtl.positions.merged$h2.SE[i]   <- gcta.res.temp[6,3]
  qtl.positions.merged$LogL[i]    <- gcta.res.temp[8,2]

  system(paste0("rm ", marker.prefix, ".snplist.", i , "*"))

}

#~~ Plot the observed and expected effects.

qtl.positions

max.var <- max(qtl.positions$Var, qtl.positions$qtl2) + 0.01

p1 <- ggplot(qtl.positions, aes(Var, qtl2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Simulated qtl2", y = "Observed qtl2") +
  coord_cartesian(ylim = c(0, max.var), xlim = c(0, max.var)) +
  ggtitle("QTL Effects on Marker Data only")

qtl.positions.merged

max.var <- max(qtl.positions.merged$Var, qtl.positions.merged$qtl2) + 0.01

p2 <- ggplot(qtl.positions.merged, aes(Var, qtl2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Simulated qtl2", y = "Observed qtl2") +
  coord_cartesian(ylim = c(0, max.var), xlim = c(0, max.var))+
  ggtitle("QTL Effects on with QTL Included")

png("results/QTL_Effect_Sizes.png", width = 8, height = 4, units = "in", res = 300)
multiplot(p1, p2, cols = 2)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Write to file                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

marker.only.model <- model1
marker.qtl.model <- model2

save(marker.only.model, marker.qtl.model, data1.gkin, gcta.res, gcta.res.merged, abeldata, abeldata,
     phenotab, qtl.positions, qtl.positions.merged, flanking.window, params, report,
     file = paste0("results/", model.name, ".RData"))

#~~ Remove files that were generated above.

remove.vec <- dir()
remove.vec <- remove.vec[grep(".bed|.fam|.bim|.phenochange|.hsq|.gen|.gz|.id|.log", remove.vec)]

if(Sys.info()["sysname"] == "Windows") {
  system(paste0("del ", paste0(remove.vec, collapse = " ")))
} else {
  system(paste0("rm ", paste0(remove.vec, collapse = " ")))
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7. Generate report                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



