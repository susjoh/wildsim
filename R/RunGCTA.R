#' RunGCTA: Run GCTA
#' @param input.string Command to be passed to GCTA. Do not include "gcta64" at
#'   the beginning of the beginning of the string e.g. RunGCTA("--help")
#' @export
#'

RunGCTA <- function(input.string){

  if(Sys.info()["sysname"] == "Windows") {

    gcta.path <- paste0(.libPaths()[1], "/wildsim/bin/windows64/gcta64.exe")

  } else {

    if(Sys.info()["sysname"] == "Linux"){
      gcta.path <- paste0(.libPaths()[1], "/wildsim/bin/linux/gcta64")
    } else {
      gcta.path <- paste0(.libPaths()[1], "/wildsim/bin/macosx/gcta64")
    }


  }

  system(paste(gcta.path, input.string))

}
