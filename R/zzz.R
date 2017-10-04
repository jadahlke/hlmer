## Global variables to be called from data frames with known variables names:
globalVariables(c("."))


## Messages to be displayed when the user loads psychmeta:
.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    packageStartupMessage("This is ", paste(pkgname, version))
    packageStartupMessage("Please report any bugs to improve functionality. \n")
}
