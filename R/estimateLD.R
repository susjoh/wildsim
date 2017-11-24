#' estimateLD: Estimated LD (r2) between adjacent markers.
#' @param markerfile Path to markerfile
#' @param snpvec Vector of SNPs to include (all others excluded). Optional
#' @param idvec Vector of IDs to include (all others excluded). Optional
#' @param full.results Default FALSE. If TRUE, prints the table of LD measures
#'   between adjacent loci.
#' @param ld.window.kb Default 1000. Maximum kb distance between which to
#'   measure LD. Higher values will ensure all loci are captured, but will take
#'   longer.
#' @export
#'


estimateLD <- function(markerfile, snpvec = NULL, idvec = NULL, full.results = F, ld.window.kb = 1000, outfile = NULL){

  marker.prefix <- gsub(".txt", "", markerfile)

  if(is.null(outfile)) outfile <- marker.prefix
  if(!is.null(snpvec)) writeLines(snpvec, paste0(outfile, ".snpvec"))
  if(!is.null(idvec )) write.table(data.frame(Family = 1, ID = idvec),
                                   paste0(outfile, ".idvec"),
                                   row.names = F, col.names = F, quote = F)

  extra.piece <- paste(c(ifelse(!is.null(snpvec), paste0(" --extract ", outfile, ".snpvec"), ""),
                         ifelse(!is.null(idvec), paste0(" --keep ", outfile, ".idvec"), "")), collapse = " ")

  snp.list.rng <- paste0(sample(letters, 10, replace = T), collapse = "")


  RunPLINK(paste0("--bfile ", marker.prefix, " --ld-window-r2 0 --ld-window-kb ", ld.window.kb, extra.piece, " --r2 yes-really --out ", snp.list.rng))

  ld.tab <- read.table(paste0(snp.list.rng, ".ld"), header = T, stringsAsFactors = F)
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


  system(paste0("rm ", snp.list.rng, ".ld"))

}
