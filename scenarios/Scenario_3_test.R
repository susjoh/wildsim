paramfile <- "Scenario_3_full.prm"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up working environment    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Load libraries

library(wildsim)
library(GenABEL)

#~~ Path to the parameter file

savefile <- gsub(".prm", ".RData", paramfile)

#~~ Parse the parameter file

paramtab <- parseParamFile(paramfile = paramfile)
paramtab$Value <- as.character(paramtab$Value)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run and Parse QMSim                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Run QMSim on the parameter file

RunQMSim(paramfile)

#~~ In this scenario, there are multiple populations. 

populations <- paramtab$Value[grep("begin_pop", paramtab$Parameter)]

#~~ Navigate to the marker and phenotype files

resdir <- paramtab$Value[which(paramtab$Parameter == "output_folder")]
temp <- strsplit(resdir, split = "")[[1]]
if(temp[length(temp)] != "/") resdir <- paste0(resdir, "/")
rm(temp)

markerfile.vec <- paste0(resdir, populations, "_mrk_001.txt")
phenofile.vec  <- paste0(resdir, populations, "_data_001.txt")

#~~ Parse the marker files into PLINK and GenABEL format.

for(i in 1:length(populations)){
  parseQMSim(markerfile = markerfile.vec[i], phenofile  = phenofile.vec[i], createGCTA = F)
}

#~~ Load the GenABEL datasets and merge them all together (keep seperate ones for GWAS)

for(i in 1:length(populations)){
  
  abeldata <- load.gwaa.data(phenofile = gsub(".txt", ".pheno", markerfile.vec[i]),
                             genofile = gsub(".txt", ".gen", markerfile.vec[i]))
  
  abeldata.merged <- load.gwaa.data(phenofile = gsub(".txt", ".pheno", markerfile.vec[i]),
                                    genofile = gsub(".txt", "_merged.gen", markerfile.vec[i]))
  
  abeldata <- add.phdata(abeldata, data.frame(id = phdata(abeldata)$id, population = populations[i], stringsAsFactors = F))
  abeldata.merged <- add.phdata(abeldata.merged, data.frame(id = phdata(abeldata.merged)$id, population = populations[i], stringsAsFactors = F))
  
  eval(parse(text = paste0("abeldata.", populations[i], " <- abeldata")))
  eval(parse(text = paste0("abeldata.merged.", populations[i], " <- abeldata.merged")))
  
  rm(abeldata, abeldata.merged)
  
}

eval(parse(text = paste0("abeldata <- abeldata.", populations[1])))
eval(parse(text = paste0("abeldata.merged <- abeldata.merged.", populations[1])))

for(i in 2:length(populations)){
  eval(parse(text = paste0("abeldata <- merge.gwaa.data(abeldata, abeldata.", populations[i], ")")))
  eval(parse(text = paste0("abeldata.merged <- merge.gwaa.data(abeldata.merged, abeldata.merged.", populations[i], ")")))
}

#~~ Merge PLINK binary files into one file:

# --bmerge [binary fileset prefix]

markerprefix.vec <- gsub(".txt", "", markerfile.vec)

markerfile <- gsub("p1", "p0", markerprefix.vec[1])
system(paste0("cp ", markerprefix.vec[1], ".bed ", markerfile, ".bed"))
system(paste0("cp ", markerprefix.vec[1], ".bim ", markerfile, ".bim"))
system(paste0("cp ", markerprefix.vec[1], ".fam ", markerfile, ".fam"))


for(i in 2:length(markerprefix.vec)){
  RunPLINK(paste0("--bfile ", markerfile, " --bmerge ", markerprefix.vec[i], " --make-bed --out ", markerprefix.vec[1]))
  system(paste0("rm ", resdir, "*~"))
}


#~~ how many polymorphic SNPs and ids?

snp.info <- data.frame(nsnps = nsnps(abeldata),
                       nids = nids(abeldata),
                       population = "all")


for(i in 1:length(populations)){
  
  eval(parse(text = paste0("snp.info <- rbind(snp.info, data.frame(nsnps = nsnps(abeldata.", populations[i], "),
             nids = nids(abeldata.", populations[i], "),
             population = \"", populations[i], "\"))")))
  
  
}

#~~ remove the abeldata files


# eval(parse(text = paste0("rm(", 
#                          paste0("abeldata.", populations, collapse = ", "), ", ",
#                          paste0("abeldata.merged.", populations, collapse = ", "),
#                          ")")))


#~~ Calculate the kinship matrix

system.time(kin.mat <- createKinshipMatrix(abeldata[,sample(snpnames(abeldata), 10000, replace = F)]))

#~~ What are the QTL IDs?

qtl.ids <- snpnames(abeldata.merged)[which(!snpnames(abeldata.merged) %in% snpnames(abeldata))]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Full Dataset                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ GWAS with marker data

model1 <- egscore(Phen, abeldata, kin = kin.mat)
model1.results <- results(model1)
model1.results$SNP.Name <- row.names(model1.results)
FullGwasSummary(model1, corrected = T)

#~~ GWAS with merged data

model2 <- egscore(Phen, abeldata.merged, kin = kin.mat)
model2.results <- results(model2)
model2.results$SNP.Name <- row.names(model2.results)
FullGwasSummary(model2, qtl.ids = qtl.ids, corrected = T)


#~~ Add inflation information to snp.info

snp.info <- cbind(snp.info, data.frame(mod1.estimate = NA,
                                       mod1.se = NA, 
                                       mod2.estimate = NA,
                                       mod2.se = NA))
snp.info[1,4:5] <- lambda(model1)
snp.info[1,6:7] <- lambda(model2)

#~~ Save significant hits and QTL hits

sig.results <- unique(rbind(subset(model2.results, P1df < 0.05/nsnps(abeldata)),
                            subset(model2.results, SNP.Name %in% qtl.ids)))

sig.results$Sig           <- ifelse(sig.results$P1df  < 0.05/nsnps(abeldata), "yes", "no")
sig.results$Sig.Corrected <- ifelse(sig.results$Pc1df < 0.05/nsnps(abeldata), "yes", "no")


#~~ Estimate genomic heritability

gcta.res <- estimateGenomicHeritability(markerfile = markerfile)
gcta.res.merged <- estimateGenomicHeritability(markerfile = markerfile, merged = T)

#~~ Estimate QTL effects

qtl.effects        <- estimateQTLEffects(markerfile = markerfile, flanking.window = 20, merged = F)
qtl.effects.merged <- estimateQTLEffects(markerfile = markerfile, flanking.window = 20, merged = T)

qtl.effects <- rbind(cbind(qtl.effects, Merged = F),
                     cbind(qtl.effects.merged, Merged = T))

rm(qtl.effects.merged)

ggplot(qtl.effects, aes(Var, qtl2)) + 
  geom_point() + 
  labs(x = "Simulated qtl2", y = "Observed qtl2") +
  geom_abline(intercept = 0, slope1) +
  coord_cartesian(ylim = c(0, max(c(qtl.effects$Var, qtl.effects$qtl2))),
                  xlim = c(0, max(c(qtl.effects$Var, qtl.effects$qtl2)))) +
  facet_wrap(~Merged)

#~~ Examine LD

system.time(ld.profile <- estimateLD(markerfile))

ld.full <- estimateLD(markerfile, full.results = T)
ggplot(ld.full, aes(Diff, R2)) + geom_point(alpha = 0.1) + stat_smooth()

#~~ Save the seed

seed <- readLines(paste0(resdir, "seed"))

#~~ Save output

save(paramtab, gcta.res, gcta.res.merged, ld.profile, qtl.effects, snp.info, sig.results, seed,
     file = savefile)

#~~ Clean up workspace

if(Sys.info()["sysname"] == "Windows") {
  
  system("cmd", input = paste0("del /F /Q ", gsub("/", "\\\\", resdir), "*"))
  system("cmd", input = paste0("rmdir ", gsub("/", "\\\\", resdir)))
  
} else {
  
  system(paste0("cp Rplots.pdf ", gsub(".prm", ".pdf", paramfile)))
  system("rm Rplots*")
  system(paste0("rm ", resdir, "/*"))
  system(paste0("rm -r ", resdir))
  
}

