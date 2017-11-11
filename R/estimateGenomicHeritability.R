#' estimateGenomicHeritability: create a kinship matrix using ibs in R
#' @param markerfile Path to markerfile
#' @param merged Default FALSE. If merged, will include the QTL markers.
#' @export
#'



estimateGenomicHeritability <- function(markerfile, merged = F){

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

    gcta.res.merged <- read.table(paste0(marker.prefix, "_merged_GRM_res.hsq"), header = T, sep = "\t", fill = T)
    gcta.res.merged

  }
}
