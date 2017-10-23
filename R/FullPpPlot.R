#' FullPpPlot: Create a P-P plot from a GenABEL scan.gwaa.object.
#' @param scan.gwaa.object GenABEL scan.gwaa.object object containing GWAS
#'   results.
#' @param corrected logical. If true, GWAS results will be reported after
#'   correction for genomic control.
#' @param lambda numeric. Genomic control parameter, to be specified if a
#'   different value from the scan.gwaa.object is to be used.
#' @param qtl.ids vector. Vector of SNP names for QTL loci (if they are included
#'   in the dataset). If specified, their positions will be highlighted.
#' @import GenABEL
#' @import ggplot2
#' @export
#'



FullPpPlot<-function(scan.gwaa.object, corrected = FALSE, lambda = NULL, qtl.ids = NULL) {

  require(ggplot2)

  cumu.object <- CumuPos(scan.gwaa.object)
  cumu.object$SNP.Name <- row.names(cumu.object)

  if(corrected){

    if(!is.null(lambda)){
      cumu.object$Pc1df <- 1-pchisq(cumu.object$chi2.1df/lambda, 1)
      warning("lambda adjusted in plot only")
    }

    gwa_res.null<-data.frame(obs=sort(-log10(cumu.object$Pc1df)),
                             exp=sort(-log10(seq(1/nrow(cumu.object),1,1/nrow(cumu.object)))))

  } else {

    gwa_res.null<-data.frame(obs=c(sort(-log10(cumu.object$P1df))),
                             exp=sort(-log10(seq(1/nrow(cumu.object),
                                                 1,
                                                 1/nrow(cumu.object)))))

  }

  p1 <- ggplot(gwa_res.null,aes(x=exp,y=obs)) +
    geom_point(colour = "red") +
    geom_abline(intercept=0,slope=1) +
    theme(axis.text.x  = element_text (size = 16, vjust = 0),
          axis.text.y  = element_text (size = 14, hjust = 1.3),
          strip.text.x = element_text (size = 16, vjust = 0.7),
          axis.title.y = element_text (size = 16, angle = 90, vjust = 0.2),
          axis.title.x = element_text (size = 16, vjust = 0.2),
          strip.background = element_blank()) +
    labs(x="Expected -log10 P",y="Observed -log10 P")

  if(!is.null(qtl.ids)){

    qtl.info <- subset(cumu.object, SNP.Name %in% qtl.ids)

    if(corrected == FALSE){
      qtl.info <- gwa_res.null[which(gwa_res.null$obs %in% -log10(qtl.info$P1df)),]
      p1 +
        geom_point(data = qtl.info, aes(x=exp,y=obs), col = "black", shape = 15, size = 2, alpha = 0.5)

    } else {

      qtl.info <- gwa_res.null[which(gwa_res.null$obs %in% -log10(qtl.info$Pc1df)),]
      p1 +
        geom_point(data = qtl.info, aes(x=exp,y=obs), col = "black", shape = 15, size = 2, alpha = 0.5)


    }

  } else {

    p1
  }

}
