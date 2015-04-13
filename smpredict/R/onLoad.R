.onLoad <- function(libname, pkgname) {
  .jpackage(pkgname, lib.loc = libname)
  packageStartupMessage('Loading smpredict (small molecule property prediction package)')
  packageStartupMessage('written by Daniel Murrell (dsmurrell@gmail.com)')
  packageStartupMessage('For a complete list of package functions, use ls("package:smpredict")')
}
