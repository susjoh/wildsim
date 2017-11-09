#' estimateAdjacentLD: Estimated LD (r2) between adjacent markers.
#' @param markerfile
#' @param full.results Default FALSE. If TRUE, prints the table of LD measures
#'   between adjacent loci.
#' @param ld.window.kb Default 1000. Maximum kb distance between which to
#'   measure LD. Higher values will ensure all loci are captured, but will take
#'   longer.
#' @export
#'


estimateLD <- function(markerfile, full.results = F, ld.window.kb = 1000){

  marker.prefix <- gsub(".txt", "", markerfile)

  RunPLINK(paste0("--bfile ", marker.prefix, " --ld-window-r2 0 --ld-window-kb ", ld.window.kb, " --r2 yes-really --out ", marker.prefix))

  ld.tab <- read.table(paste0(marker.prefix, ".ld"), header = T, stringsAsFactors = F)
  ld.tab$SNP_A <- as.numeric(gsub("M", "", ld.tab$SNP_A))
  ld.tab$SNP_B <- as.numeric(gsub("M", "", ld.tab$SNP_B))

  ld.tab$Diff <- ld.tab$BP_B - ld.tab$BP_A
  ld.tab$Diff <- ifelse(ld.tab$Diff > 0, ld.tab$Diff, ld.tab$Diff * -1)
  ld.tab$LogR2 <- -log10(ld.tab$R2)
  ld.tab$LogR2[which(ld.tab$LogR2 == Inf)] <- NA

  temp <- lm(ld.tab$LogR2 ~ ld.tab$Diff)

  if(full.results) {
    return(ld.tab)
  } else {
    return(data.frame(mean.LD = mean(ld.tab$R2),
                      median.LD = median(ld.tab$R2),
                      mean.BP = mean(ld.tab$Diff),
                      median.BP = median(ld.tab$Diff),
                      N.ld.measures = nrow(ld.tab),
                      log10.intercept = temp$coefficients[1],
                      log10.slope = temp$coefficients[2]))
  }


  system(paste0("rm ", marker.prefix, ".ld"))

}
