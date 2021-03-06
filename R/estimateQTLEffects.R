#' estimateQTLEffects: estimate variance attributed to simulated QTLs
#' @param markerfile string
#' @param phenofile string
#' @param gwaa.data GenABEL gwaa.data set, optional. If specified, map positions
#'   are taken from here.
#' @param idvec Vector of IDs to include (all others excluded). Optional
#' @param flanking.window integer. Number of SNPs flanking the QTL for
#'   estimation of effect sizes
#' @param merged Boolean, FALSE. If TRUE, include the quantitative trait
#'   nucleotide from the simulation.
#' @export
#'

estimateQTLEffects <- function(markerfile, phenofile, gwaa.data = NULL, idvec = NULL, flanking.window, outfile = NULL, merged = F, variance.threshold = 0){

  require(plyr)
  require(evaluate)

  marker.prefix <- gsub(".txt", "", markerfile)
  temp <- gsub("\\\\", "/", marker.prefix)
  temp <- strsplit(temp, "/")[[1]]
  temp2 <- strsplit(temp[length(temp)], split = "_")[[1]]
  temp2 <- paste0("lm_qtl_", temp2[length(temp2)])
  temp <- paste(temp[-length(temp)], collapse = "/")
  tempname <- paste(c(temp, temp2), collapse = "/")

  if(is.null(outfile)) outfile <- marker.prefix
  if(!is.null(idvec )) write.table(data.frame(Family = 1, ID = idvec),
                                   paste0(outfile, ".idvec"),
                                   row.names = F, col.names = F, quote = F)

  extra.piece <- paste(c(ifelse(!is.null(snpvec), paste0(" --extract ", outfile, ".snpvec"), ""),
                         ifelse(!is.null(idvec), paste0(" --keep ", outfile, ".idvec"), "")), collapse = " ")



  qtl.positions <- read.table(paste0(tempname, ".txt"), header = T, stringsAsFactors = F)
  qtl.effects   <- read.table(paste0(gsub("_mrk", "_freq_qtl", marker.prefix), ".txt"), header = T, fill = T)
  names(qtl.effects) <- c("ID", "Gen", "Chr", "Var", "Allele.Freq.1", "Allele.Freq.2")
  qtl.positions$BP <- qtl.positions$Position * 1e6
  qtl.positions <- join(qtl.positions, qtl.effects)

  qtl.positions
  rm(temp, temp2)


  #~~ Are there multiple generations in the QTL effect file? At this stage, just take the weighted mean...

  if(any(table(qtl.positions$ID) > 1)){

    phenofile <- read.table(phenofile, header = T)
    head(phenofile)
    temp <- data.frame(tapply(phenofile$QTL, phenofile$G, length))
    names(temp) <- "N"
    temp$Gen <- row.names(temp)

    qtl.positions <- join(qtl.positions, temp)
    rm(temp)

    new.qtl <- NULL

    qtl.positions$Allele.Freq.1 <- as.character(qtl.positions$Allele.Freq.1)
    qtl.positions$Allele.Freq.2 <- as.character(qtl.positions$Allele.Freq.2)

    two.vec <- grep("2:", qtl.positions$Allele.Freq.1)
    if(length(two.vec > 0)){
      qtl.positions$Allele.Freq.2[two.vec] <- qtl.positions$Allele.Freq.1[two.vec]
      qtl.positions$Allele.Freq.1[two.vec] <- ""
    }

    for(i in unique(qtl.positions$ID)){

      temp <- subset(qtl.positions, ID == i)
      temp$Allele.Freq.1 <- ifelse(temp$Allele.Freq.1 == "", 0, temp$Allele.Freq.1)
      temp$Allele.Freq.2 <- ifelse(temp$Allele.Freq.2 == "", 0, temp$Allele.Freq.2)
      temp2 <- data.frame(Var = mean(temp$Var, weights = temp$N),
                          Allele.Freq.1 = mean(as.numeric(gsub("1:", "", temp$Allele.Freq.1)), weights = temp$N),
                          Allele.Freq.2 = mean(as.numeric(gsub("2:", "", temp$Allele.Freq.2)), weights = temp$N),
                          TotalN = sum(temp$N))
      temp2 <- cbind(temp[1, c("ID", "Chr", "Position", "BP", "Gen")], temp2)
      new.qtl <- rbind(new.qtl, temp2)
      rm(temp, temp2)

    }

    qtl.positions <- new.qtl
    rm(new.qtl)
  }

  #~~ First, run the regions in the dataset that does not include the QTL.

  qtl.positions$Flanking.Window <- flanking.window
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
  qtl.positions$qtl.error <- NA


  if(is.null(gwaa.data)){
    RunPLINK(paste0("--bfile ", marker.prefix, " --recode --out ", marker.prefix))
    map.file <- read.table(paste0(marker.prefix, ".map"))
    system(paste0("rm ", marker.prefix, ".ped"))
    system(paste0("rm ", marker.prefix, ".map"))

    head(map.file)
    map.file <- map.file[,c(2, 1, 3, 4)]
    names(map.file) <- c("ID", "Chr", "Position", "BP")
  } else {
    map.file <- data.frame(ID = snpnames(gwaa.data),
                           Chr = chromosome(gwaa.data),
                           Position = map(gwaa.data)/1e6,
                           BP = map(gwaa.data))

  }

  if(merged){
    map.file <- rbind(map.file, qtl.positions[,1:4])
    map.file <- arrange(map.file, Chr, Position)
  }

  map.file <- unique(map.file)

  for(i in 1:nrow(qtl.positions)){


    x <- subset(map.file, Chr == qtl.positions$Chr[i])

    x <- arrange(x, Position)

    focal.snp <- which(x$BP > qtl.positions$BP[i])[1]

    qtl.positions$QTL.Position[i] <- as.character(x$ID[focal.snp])

    if(qtl.positions$Var[i] > variance.threshold){


      if(merged){
        start.pos <- (focal.snp - flanking.window)
      } else {
        start.pos <- (focal.snp - flanking.window + 1)
      }
      stop.pos  <- (focal.snp + flanking.window)


      if(start.pos < 1){
        start.pos <- 1
        stop.pos <- flanking.window*2
      }

      if(stop.pos > nrow(x)){
        start.pos <- nrow(x) - (flanking.window*2)
        stop.pos <- nrow(x)
      }

      snp.list <- as.character(x$ID[start.pos:stop.pos])
      if(merged == F & length(grep("Q", snp.list) > 0)) snp.list <- snp.list[-grep("Q", snp.list)]


      rest.list <- unique(as.character(subset(map.file, !ID %in% snp.list)$ID))
      if(merged == F & length(grep("Q", rest.list) > 0)) rest.list <- rest.list[-grep("Q", rest.list)]


      snp.list.rng <- paste0(".", paste0(sample(letters, 10, replace = T), collapse = ""), collapse = "")

      write.table(snp.list, paste0(marker.prefix, snp.list.rng, ".snplist.", i, ".txt"), row.names = F, col.names = F, quote = F)
      write.table(rest.list, paste0(marker.prefix, snp.list.rng, ".restlist.", i, ".txt"), row.names = F, col.names = F, quote = F)


      RunGCTA(paste0("--bfile ", marker.prefix, "_merged --make-grm-gz --keep ", outfile, ".idvec --extract ", marker.prefix, snp.list.rng, ".snplist.", i, ".txt --out ", marker.prefix, snp.list.rng, ".snplist.", i , "_GRM"))
      RunGCTA(paste0("--bfile ", marker.prefix, "_merged --make-grm-gz --keep ", outfile, ".idvec --extract ", marker.prefix, snp.list.rng, ".restlist.", i, ".txt --out ", marker.prefix, snp.list.rng, ".snplist.", i , "_wo_GRM"))

      writeLines(paste0(marker.prefix, snp.list.rng, ".snplist.", i , "_GRM\n",
                        marker.prefix, snp.list.rng, ".snplist.", i , "_wo_GRM"),
                 paste0(marker.prefix, snp.list.rng, ".snplist.", i , "_GRMs.txt"))


      y <- evaluate("RunGCTA(paste0(\"--mgrm-gz \", marker.prefix, snp.list.rng, \".snplist.\", i , \"_GRMs.txt --reml-maxit 1000  --pheno \", marker.prefix, \".phenochange.txt --reml --out \", marker.prefix, snp.list.rng,  \"_GRMs_res\"))")

      if(length(y) == 1) qtl.positions$qtl.error[i] <- "none"

      if(length(y) == 2){
        y <- evaluate("RunGCTA(paste0(\"--mgrm-gz \", marker.prefix, snp.list.rng,  \".snplist.\", i , \"_GRMs.txt --reml-no-constrain --reml-maxit 1000  --pheno \", marker.prefix, \".phenochange.txt --reml --out \", marker.prefix, snp.list.rng,  \"_GRMs_res\"))")
      }

      if(length(y) == 1) qtl.positions$qtl.error[i] <- "none"
      if(length(y) == 2) qtl.positions$qtl.error[i] <- "boundary"



      if(file.exists(paste0(marker.prefix, snp.list.rng, "_GRMs_res.hsq"))){
        gcta.res.temp <- read.table(paste0(marker.prefix, snp.list.rng, "_GRMs_res.hsq"), header = T, sep = "\t", fill = T)

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

      }

      system("rm Rplots*")
      system(paste0("rm ", marker.prefix, snp.list.rng, ".snplist.", i , "*"))
      system(paste0("rm ", marker.prefix, snp.list.rng, ".restlist.", i , "*"))
      system(paste0("rm ", marker.prefix, snp.list.rng, "_GRMs_res.hsq"))
    }
  }
  qtl.positions

}


