#' FullSummary: Create a Manhatten and P-P plots and return top results from a
#' GenABEL scan.gwaa.object.
#' @param scan.gwaa.object GenABEL scan.gwaa.object object containing GWAS
#'   results.
#' @param corrected logical. If true, GWAS results will be reported after
#'   correction for genomic control.
#' @param bonf numeric. Threshold P-value for genome-wide significance. If not
#'   defined, will be calculated using a Bonferroni correction.
#' @param qtl.ids vector. Vector of SNP names for QTL loci (if they are included
#'   in the dataset). If specified, their positions will be highlighted.
#' @import GenABEL
#' @import ggplot2
#' @export
#'


FullGwasSummary <- function(scan.gwaa.object, corrected = FALSE, bonf = F, qtl.ids = NULL) {

  cumutemp <- CumuPos(scan.gwaa.object)

  if(bonf == F) bonf = 0.05/nrow(cumutemp)

  print(summary(scan.gwaa.object))
  #print(lambda(scan.gwaa.object))

  if(corrected == FALSE){

    if(is.null(qtl.ids)){

      p3 <- FullGwasPlot (scan.gwaa.object, corrected = FALSE) + labs(title = "Uncorrected GWAS")
      p4 <- FullPpPlot (scan.gwaa.object,  corrected = FALSE) + labs(title = paste0("lambda = ", lambda(scan.gwaa.object)))

    } else {

      p3 <- FullGwasPlot (scan.gwaa.object, corrected = FALSE, qtl.ids = qtl.ids) + labs(title = "Uncorrected GWAS")
      p4 <- FullPpPlot (scan.gwaa.object,  corrected = FALSE, qtl.ids = qtl.ids) + labs(title = paste0("lambda = ", lambda(scan.gwaa.object)))


    }

  } else {

    if(is.null(qtl.ids)){
      p3 <- FullGwasPlot (scan.gwaa.object, corrected = TRUE) + labs(title = "Corrected GWAS")
      p4 <- FullPpPlot (scan.gwaa.object,  corrected = TRUE) + labs(title = "lambda = 1")
    } else {
      p3 <- FullGwasPlot (scan.gwaa.object, corrected = TRUE, qtl.ids = qtl.ids) + labs(title = "Corrected GWAS")
      p4 <- FullPpPlot (scan.gwaa.object,  corrected = TRUE, qtl.ids = qtl.ids) + labs(title = "lambda = 1")

    }
  }

  multiplot(p3, p4, cols = 2, layout = matrix(c(1,1,2), 1, 3))

}
