#' createKinshipMatrix: create a kinship matrix using ibs in R
#' @param gwaa.data GenABEL gwaa.data object.
#' @param plot.mds Default TRUE. Creates a plot of classical multidimensional
#'   scaling (MDS) of the kinship matrix.
#' @export
#' 

createKinshipMatrix <- function(gwaa.data = NULL, plot.mds = T){
  
  data1.gkin <- ibs(gwaa.data[, autosomal(gwaa.data)], weight="freq")
  data1.dist <- as.dist(0.5 - data1.gkin)
  data1.mds <- cmdscale(data1.dist)
  
  if(plot.mds) plot(data1.mds)
  
  data1.gkin

}



