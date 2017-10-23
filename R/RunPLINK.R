#' RunPLINK: Run PLINK
#' @param input.string Command to be passed to PLINK. Do not include "plink" at
#'   the beginning of the beginning of the string e.g. RunPLINK("--help")
#' @export
#'

RunPLINK <- function(input.string){

  if(Sys.info()["sysname"] == "Windows") {

    plink.path <- paste0(.libPaths()[1], "/wildsim/bin/windows64/plink.exe")

  } else {

    if(Sys.info()["sysname"] == "Linux"){
      plink.path <- paste0(.libPaths()[1], "/wildsim/bin/linux/plink")
    } else {
      plink.path <- paste0(.libPaths()[1], "/wildsim/bin/macos/plink")
    }


  }

  system(paste(plink.path, input.string))

}
