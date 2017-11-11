#' estimateGenomicHeritability: create a kinship matrix using ibs in R
#' @param markerfile Path to markerfile
#' @param merged Default FALSE. If merged, will include the QTL markers.
#' @import reshape
#' @import plyr
#' @export
#'

estimateGenomicHeritability <- function(markerfile, merged = F){

  require(reshape)
  require(plyr)

  marker.prefix <- gsub(".txt", "", markerfile)

  if(!merged){

    if(!file.exists(paste0(marker.prefix, "_GRM.grm.gz"))){

      RunGCTA(paste0("--bfile ", marker.prefix, " --make-grm-gz --out ", marker.prefix, "_GRM"))

    }

    RunGCTA(paste0("--grm-gz ", marker.prefix, "_GRM --pheno ",
                   marker.prefix, ".phenochange.txt --reml --out ",
                   marker.prefix, "_GRM_res"))

    gcta.res <- read.table(paste0(marker.prefix, "_GRM_res.hsq"), header = T, sep = "\t", fill = T)
    gcta.res
  } else {

    if(!file.exists(paste0(marker.prefix, "_merged_GRM.grm.gz"))){

      RunGCTA(paste0("--bfile ", marker.prefix, "_merged --make-grm-gz --out ", marker.prefix, "_merged_GRM"))

    }

    RunGCTA(paste0("--grm-gz ", marker.prefix, "_merged_GRM --pheno ",
                   marker.prefix, ".phenochange.txt --reml --out ",
                   marker.prefix, "_merged_GRM_res"))

    gcta.res <- read.table(paste0(marker.prefix, "_merged_GRM_res.hsq"), header = T, sep = "\t", fill = T)
    gcta.res

  }

  gcta.res <- melt(gcta.res, id.vars = "Source")
  gcta.res$Source <- as.character(gcta.res$Source)
  gcta.res <- subset(gcta.res, !is.na(value))
  gcta.res$Source[which(gcta.res$variable == "SE")] <- paste0(gcta.res$Source[which(gcta.res$variable == "SE")], ".SE")
  gcta.res <- arrange(gcta.res, Source)
  gcta.res <- gcta.res[, -2]
  name.vec <- gcta.res[,1]
  gcta.res <- data.frame(t(gcta.res[,2]))
  names(gcta.res) <- name.vec
  gcta.res

}
