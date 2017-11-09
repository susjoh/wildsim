#' RunQMSim: Run QMSim
#' @param input.string Command to be passed to QMSim Do not include "gcta64" at
#'   the beginning of the beginning of the string e.g. RunGCTA("--help")
#' @export
#'

RunQMSim <- function(input.string){

  if(Sys.info()["sysname"] == "Windows") {

    qmsim.path <- paste0(.libPaths()[1], "/wildsim/bin/windows64/QMSim.exe")

  } else {

    if(Sys.info()["sysname"] == "Linux"){
      qmsim.path <- paste0(.libPaths()[1], "/wildsim/bin/linux/QMSim")
    } else {
      qmsim.path <- paste0(.libPaths()[1], "/wildsim/bin/macosx/QMSim")
    }


  }

  system(paste(qmsim.path, input.string))

}
