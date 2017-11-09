#' parseQMSim: Parse files in a QMSIM directory. This function converst the
#' simulated data to GenABEL and PLINK format files and creates relatedness
#' matrices for the full dataset using GCTA (optional)
#' @param markerfile string. Default is "p1_mrk_001.txt".
#' @param phenofile string. Default is "p1_data_001.txt"
#' @param createGCTA boolean. Default TRUE. Creates a relatedness matrix from
#'   the full data set.
#' @param createQTLMerged boolean. Default TRUE. Creates files with ("_merged")
#'   and without the simulated quantitative trait loci.
#' @export
#'

parseQMSim <- function(markerfile = "p1_mrk_001.txt",
                        phenofile = "p1_data_001.txt",
                        createGCTA = T,
                        createQTLMerged = T){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 0. Load libraries, functions and data information   #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


  require(GenABEL)
  require(plyr)
  require(ggplot2)

  current.wd <- getwd()

  marker.dir <- gsub("\\\\", "/", markerfile)
  marker.dir <- strsplit(marker.dir, split = "/")[[1]]
  if(length(marker.dir > 1)) marker.dir <- paste(marker.dir[1:(length(marker.dir)-1)], collapse = "/")

  pheno.dir <- gsub("\\\\", "/", phenofile)
  pheno.dir <- strsplit(pheno.dir, split = "/")[[1]]
  if(length(pheno.dir > 1)) pheno.dir <- paste(pheno.dir[1:(length(pheno.dir)-1)], collapse = "/")

  if(marker.dir != pheno.dir) stop("markerfile and phenofile should be in the same directory.")

  setwd(marker.dir)

  markerfile <- strsplit(gsub("\\\\", "/", markerfile), split = "/")[[1]]
  markerfile <- markerfile[length(markerfile)]

  phenofile <- strsplit(gsub("\\\\", "/", phenofile), split = "/")[[1]]
  phenofile <- phenofile[length(phenofile)]

  marker.prefix <- gsub(".txt", "", markerfile)
  qtl.prefix <- gsub("mrk", "qtl", marker.prefix)

  report <- readLines("report.txt")
  params <- readLines(grep(".prm", dir(), value = T))



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

  if(Sys.info()["sysname"] == "Windows") {

    system("cmd", input = paste0("sed \"1d\" ", markerfile,
                                 " | sed -e \"s/  */ /g\"  | sed -e \"s/^/1 /g\" | sed \"s/ / 0 0 0 /2\" > ",
                                 marker.prefix, ".ped"), show.output.on.console = F)
  } else {

    system(paste0("sed \"1d\" ", markerfile,
                          " | sed -e \"s/  */ /g\"  | sed -e \"s/^/1 /g\" | sed \"s/ / 0 0 0 /2\" > ",
                          marker.prefix, ".ped"))

  }

  #~~ Format qtlfile. Need to scrub first line, add family, and replace multiple spaces with a single space.


  if(Sys.info()["sysname"] == "Windows") {

    system("cmd", input = paste0("sed \"1d\" ", gsub("mrk", "qtl", markerfile),
                                 " | sed -e \"s/  */ /g\"  | sed -e \"s/^/1 /g\" | sed \"s/ / 0 0 0 /2\" > ",
                                 qtl.prefix, ".ped"), show.output.on.console = F)

  } else{

    system(paste0("sed \"1d\" ", gsub("mrk", "qtl", markerfile),
                          " | sed -e \"s/  */ /g\"  | sed -e \"s/^/1 /g\" | sed \"s/ / 0 0 0 /2\" > ",
                          qtl.prefix, ".ped"))
  }


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
  if(createGCTA) RunGCTA(paste0("--bfile ", marker.prefix, " --make-grm-gz --out ", marker.prefix, "_GRM"))

  RunPLINK(paste0("--file ", marker.prefix, "_merged --make-bed --out ", marker.prefix, "_merged"))
  if(createGCTA) RunGCTA(paste0("--bfile ", marker.prefix, "_merged --make-grm-gz --out ", marker.prefix, "_merged_GRM"))

  #~~ Convert to GenABEL format for unmerged and merged files

  maptab <- read.table(paste0(marker.prefix, ".map"))
  maptab <- maptab[,-3]
  write.table(maptab, paste0(marker.prefix, ".map"), row.names = F, col.names = F, quote = F)

  convert.snp.ped(pedfile = paste0(marker.prefix, ".ped"),
                  mapfile = paste0(marker.prefix, ".map"),
                  outfile = paste0(marker.prefix, ".gen"),
                  traits = 1,
                  mapHasHeaderLine = F)

  if(createQTLMerged){
    mergemaptab <- read.table(paste0(marker.prefix, "_merged.map"))
    mergemaptab <- mergemaptab[,-3]
    write.table(mergemaptab, paste0(marker.prefix, "_merged.map"), row.names = F, col.names = F, quote = F)



    convert.snp.ped(pedfile = paste0(marker.prefix, "_merged.ped"),
                    mapfile = paste0(marker.prefix, "_merged.map"),
                    outfile = paste0(marker.prefix, "_merged.gen"),
                    traits = 1,
                    mapHasHeaderLine = F)

  }


  if(Sys.info()["sysname"] == "Windows") {

    system("cmd", input = "del *.ped")
    system("cmd", input = "del *.map")
    system("cmd", input = "del *.sexchange.txt")
    system("cmd", input = "del *.nosex")

  } else {

    system("rm *.ped")
    system("rm *.map")
    system("rm *.sexchange.txt")
    system("rm *.nosex")
  }



  setwd(current.wd)

}


