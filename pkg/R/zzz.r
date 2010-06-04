.First.lib <- function(libname, pkgname)
{
  library.dynam("CompRandFld", package = pkgname, lib.loc = libname)
  invisible()
}


