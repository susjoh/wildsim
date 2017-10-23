#' FullGwasPlot: Create a Manhattan plot from a GenABEL scan.gwaa.object.
#' @param scan.gwaa.object GenABEL scan.gwaa.object object containing GWAS
#'   results.
#' @param corrected logical. If true, GWAS results will be reported after
#'   correction for genomic control.
#' @param bonf numeric. Threshold P-value for genome-wide significance. If not
#'   defined, will be calculated using a Bonferroni correction.
#' @param lambda numeric. Genomic control parameter, to be specified if a
#'   different value from the scan.gwaa.object is to be used.
#' @param qtl.ids vector. Vector of SNP names for QTL loci (if they are included
#'   in the dataset). If specified, their positions will be highlighted.
#' @import GenABEL
#' @import ggplot2
#' @export
#'



FullGwasPlot<-function(scan.gwaa.object, corrected = FALSE, bonf = NULL, lambda = NULL, qtl.ids = NULL) {

  require(ggplot2)

  cumu.object <- CumuPos(scan.gwaa.object)

  chrinfo <- NULL

  if(is.null(bonf)) bonf = 0.05/nrow(cumu.object)

  for(i in unique(cumu.object$Chromosome)){

    temp1 <- subset(cumu.object, Chromosome == i)

    temp2 <- data.frame(Chromosome = i,
                        Start = temp1[1,"Cumu"],
                        Stop = temp1[nrow(temp1),"Cumu"])

    chrinfo <- rbind(chrinfo, temp2)
  }

  chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)

  colourscale <- rep(c("red","blue"), times = length(unique(chrinfo$Chromosome)))
  colourscale <- colourscale[1:(length(colourscale)/2)]

  if(corrected){

    if(!is.null(lambda)){
      cumu.object$Pc1df <- 1-pchisq(cumu.object$chi2.1df/lambda, 1)
      warning("lambda adjusted in plot only")
    }

    p1 <- ggplot(cumu.object, aes(Cumu, -log10(Pc1df), col = factor(Chromosome))) +
      geom_point(alpha = 0.4) +
      geom_hline(yintercept = -log10(bonf),linetype = 2, alpha = 0.6, size = 1) +
      scale_colour_manual(values = colourscale) +
      theme(legend.position = "none") +
      theme(axis.text.x  = element_text (size = 16, vjust = 0),
            axis.text.y  = element_text (size = 14, hjust = 1.3),
            strip.text.x = element_text (size = 16, vjust = 0.7),
            axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
            axis.title.x = element_text (size = 16, vjust = 0.2),
            strip.background = element_blank()) +
      scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chromosome) +
      labs(x = "Chromosome", y = "-log10 P")

  } else {

    p1 <- ggplot(cumu.object, aes(Cumu,-log10(P1df), col = factor(Chromosome))) +
      geom_point(alpha = 0.4) +
      geom_hline(yintercept=-log10(bonf),linetype=2, alpha = 0.6, size = 1) +
      scale_colour_manual(values = colourscale) +
      theme(legend.position="none") +
      theme(axis.text.x  = element_text (size = 16, vjust = 0),
            axis.text.y  = element_text (size = 14, hjust = 1.3),
            strip.text.x = element_text (size = 16, vjust = 0.7),
            axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
            axis.title.x = element_text (size = 16, vjust = 0.2),
            strip.background = element_blank()) +
      scale_x_continuous(breaks=chrinfo$Mid,labels=chrinfo$Chromosome) +
      labs(x="Chromosome",y="-log10 P")
  }

  if(!is.null(qtl.ids)){
    cumu.object$SNP.Name <- row.names(cumu.object)
    qtl.info <- subset(cumu.object, SNP.Name %in% qtl.ids)

    if(corrected){
      p1 +
        geom_point(data = qtl.info, aes(Cumu, -log10(Pc1df)), col = "black", shape = 15, size = 2, alpha = 0.5) +
        geom_vline(xintercept = qtl.info$Cumu, alpha = 0.2)
    } else {
      p1 +
        geom_point(data = qtl.info, aes(Cumu, -log10(P1df)), col = "black", shape = 15, size = 2, alpha = 0.5) +
        geom_vline(xintercept = qtl.info$Cumu, alpha = 0.2)

    }
  } else {
    p1

  }


}
