getbasistype <- function(basisfd) {
  #  Extracts the type of basis, permitting variants in spelling

  #  Last modified 6 Feb 2001
  
  if (!(inherits(basisfd, "basis.fd"))) stop(
    "First argument is not a basis object.")

  type <- basisfd$type

  if        (type == "Fourier" ||
             type == "fourier" ||
             type == "Fou"     ||
             type == "fou") {
                return("fourier")
  } else if (type == "bspline" ||
             type == "Bspline" ||
             type == "Bsp"     ||
             type == "bsp") {
                return("bspline")
  } else if (type == "poly"    ||
             type == "pol"     ||
             type == "polynomial") {
                return("poly")
  } else if (type == "exp"     ||
             type == "expon"   ||
             type == "exponential") {
                return("expon")
  } else if (type == "polygonal" ||
             type == "polyg"     ||
             type == "polygon") {
                return("polyg")
  } else if (type == "power"   ||
             type == "pow") {
                return("power")
  } else if (type == "const"   ||
             type == "con"     ||
             type == "constant") {
                return("const")
  } else {
                return("unknown")
  }
}
