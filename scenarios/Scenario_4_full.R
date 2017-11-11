paramfile <- "Scenario_4_full.prm"

#~~ Load libraries

library(wildsim)
library(GenABEL)

#~~ Path to the parameter file

savefile <- gsub(".prm", ".RData", paramfile)

#~~ Parse the parameter file

paramtab <- parseParamFile(paramfile = paramfile)
paramtab$Value <- as.character(paramtab$Value)

#~~ Run QMSim on the parameter file

RunQMSim(paramfile)

#~~ Navigate to the marker and phenotype files

resdir <- paramtab$Value[which(paramtab$Parameter == "output_folder")]
temp <- strsplit(resdir, split = "")[[1]]
if(temp[length(temp)] != "/") resdir <- paste0(resdir, "/")
rm(temp)

markerfile <- paste0(resdir, "p1_mrk_001.txt")
phenofile  <- paste0(resdir, "p1_data_001.txt")

#~~ Parse the marker file into PLINK and GenABEL format.

parseQMSim(markerfile = markerfile, phenofile  = phenofile)     #288.15

#~~ Load the GenABEL datasets

abeldata <- load.gwaa.data(phenofile = gsub(".txt", ".pheno", markerfile),
                           genofile = gsub(".txt", ".gen", markerfile))

abeldata.merged <- load.gwaa.data(phenofile = gsub(".txt", ".pheno", markerfile),
                                  genofile = gsub(".txt", "_merged.gen", markerfile))

#~~ how many polymorphic SNPs and ids?

snp.info <- data.frame(nsnps = nsnps(abeldata),
                       nids = nids(abeldata))


#~~ Calculate the kinship matrix

system.time(kin.mat <- createKinshipMatrix(abeldata))

#~~ What are the QTL IDs?

qtl.ids <- snpnames(abeldata.merged)[which(!snpnames(abeldata.merged) %in% snpnames(abeldata))]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Full Dataset                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ GWAS with marker data

model1 <- egscore(Phen, abeldata, kin = kin.mat)
model1.results <- results(model1)
model1.results$SNP.Name <- row.names(model1.results)
FullGwasSummary(model1)

#~~ GWAS with merged data

model2 <- egscore(Phen, abeldata.merged, kin = kin.mat)
model2.results <- results(model2)
model2.results$SNP.Name <- row.names(model2.results)
FullGwasSummary(model2, qtl.ids = qtl.ids)


#~~ Add inflation information to snp.info

snp.info <- cbind(snp.info, lambda(model1))
names(snp.info)[3:4] <- c("mod1.estimate", "mod1.se")

snp.info <- cbind(snp.info, lambda(model2))
names(snp.info)[5:6] <- c("mod2.estimate", "mod2.se")

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

