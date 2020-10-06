#' jheem: A package
#'
#'@docType package
#'@name jheem
#'@useDynLib jheem, .registration = TRUE
#'@import Rcpp
NULL

.onUnload <- function (libpath) { library.dynam.unload("jheem", libpath) }
#if you get an error about loading jheem.dll, need to close Rstudio to release the file lock, then re-open
# (?possibly delete Documents/R/win-library/3.5/jheem/libs/x64/jheem.dll)
