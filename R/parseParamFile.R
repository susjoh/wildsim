#' parseParamFile: Parse the parameter file.
#' @param paramfile Path to the parameter file
#' @param markerfile If provided, then will assume param file is in the same
#'   directory and parse the *.prm file.
#' @export
#'


parseParamFile <- function(paramfile, markerfile = NULL){

  if(!is.null(markerfile)){
    marker.dir <- gsub("\\\\", "/", markerfile)
    marker.dir <- strsplit(marker.dir, split = "/")[[1]]
    if(length(marker.dir > 1)) marker.dir <- paste(marker.dir[1:(length(marker.dir)-1)], collapse = "/")

    paramfile <- readLines(paste0(marker.dir, "/", grep("*.prm", dir(marker.dir), value = T)))

  } else {

    paramfile <- readLines(paramfile)

  }


  paramfile <- paramfile[grep("=", paramfile)]
  paramfile <- gsub(" ", "", paramfile)
  paramfile <- gsub("\t", "", paramfile)

  paramfile <- sapply(paramfile, function(x) strsplit(x, split = "//")[[1]][1], USE.NAMES = F)
  paramfile <- sapply(paramfile, function(x) strsplit(x, split = ",pop=")[[1]][1], USE.NAMES = F)
  paramfile <- gsub("[n", ".n", paramfile, fixed = T)
  paramfile <- gsub(";", "", paramfile, fixed = T)

  data.frame(Parameter = sapply(paramfile, function(x) strsplit(x, split = "=")[[1]][1], USE.NAMES = F),
             Value     = sapply(paramfile, function(x) strsplit(x, split = "=")[[1]][2], USE.NAMES = F))

}
