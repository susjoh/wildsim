#' CumuPos: Determine the cumulative positions of scan.gwaa.object for plotting.
#' @param scan.gwaa.object GenABEL scan.gwaa.object object containing GWAS results.
#' @import GenABEL
#' @export
#'


CumuPos<-function(scan.gwaa.object) {

  #~~ Extract results and calculate cumulative positions
  gwa_res <- results(scan.gwaa.object)

  gwa_res$Chromosome <- as.character(gwa_res$Chromosome)

  if("X" %in% gwa_res$Chromosome){
    chrvec <- unique(gwa_res$Chromosome)

    if(length(chrvec > 1)){
      chrvec <- as.numeric(chrvec[-which(chrvec == "X")])
      gwa_res$Chromosome[which(gwa_res$Chromosome == "X")] <- max(chrvec) + 1
    }
  }

  gwa_res$Chromosome <- as.numeric(gwa_res$Chromosome)

  gwa_res <- gwa_res[with(gwa_res, order(Chromosome,Position)), ]

  gwa_res$Diff <- c(0,diff(gwa_res$Position))
  gwa_res$Diff[gwa_res$Diff < 0] <- 1

  gwa_res$Cumu <- cumsum(gwa_res$Diff)

  gwa_res$Cumu2 <- gwa_res$Cumu + (25000000 * gwa_res$Chromosome)

  #~~ Plot the uncorrected p-values
  return(gwa_res)
}
