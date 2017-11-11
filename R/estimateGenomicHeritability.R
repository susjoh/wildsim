#' estimateGenomicHeritability: create a kinship matrix using ibs in R
#' @param markerfile Path to markerfile
#' @param snpvec Vector of SNPs to include (all others excluded). Optional
#' @param idvec Vector of IDs to include (all others excluded). Optional
#' @param merged Default FALSE. If merged, will include the QTL markers.
#' @import reshape
#' @import plyr
#' @export
#'

estimateGenomicHeritability <- function(markerfile, snpvec = NULL, idvec = NULL, outfile = NULL, merged = F){

  require(reshape)
  require(plyr)

  marker.prefix <- gsub(".txt", "", markerfile)

  if(is.null(outfile)) outfile <- marker.prefix
  if(!is.null(snpvec)) writeLines(snpvec, paste0(marker.prefix, ".snpvec"))
  if(!is.null(idvec )) write.table(data.frame(Family = 1, ID = idvec),
                                   paste0(marker.prefix, ".idvec"),
                                   row.names = F, col.names = F, quote = F)

  extra.piece <- paste(c(ifelse(!is.null(snpvec), paste0(" --extract ", marker.prefix, ".snpvec"), ""),
                         ifelse(!is.null(idvec), paste0(" --keep ", marker.prefix, ".idvec"), "")), collapse = " ")

  if(!merged){

    if(!file.exists(paste0(outfile, "_GRM.grm.gz"))){

      RunGCTA(paste0("--bfile ", marker.prefix, extra.piece, " --make-grm-gz --out ", outfile, "_GRM"))

    }

    RunGCTA(paste0("--grm-gz ", outfile, "_GRM --pheno ",
                   marker.prefix, ".phenochange.txt --reml --out ",
                   outfile, "_GRM_res"))

    gcta.res <- read.table(paste0(outfile, "_GRM_res.hsq"), header = T, sep = "\t", fill = T)
    gcta.res
  } else {

    if(!file.exists(paste0(outfile, "_merged_GRM.grm.gz"))){

      RunGCTA(paste0("--bfile ", marker.prefix, "_merged ", extra.piece, " --make-grm-gz --out ", outfile, "_merged_GRM"))

    }

    RunGCTA(paste0("--grm-gz ", outfile, "_merged_GRM --pheno ",
                   marker.prefix, ".phenochange.txt --reml --out ",
                   outfile, "_merged_GRM_res"))

    gcta.res <- read.table(paste0(outfile, "_merged_GRM_res.hsq"), header = T, sep = "\t", fill = T)
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
