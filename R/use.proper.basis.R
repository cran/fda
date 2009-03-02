use.proper.basis <- function(type) {
  #  recognizes type of basis by use of several variant spellings
  if(type == "bspline" ||
          type == "Bspline" ||
          type == "spline"  ||
          type == "Bsp"     ||
          type == "bsp") {
                return("bspline")
        }
  else if(type == "con"      ||
          type == "const"    ||
          type == "constant") {
                return("const")
        }
  else if(type == "exp"    ||
          type == "expon"  ||
          type == "exponential") {
                return("expon")
        }
  else if(type == "Fourier" ||
     type == "fourier" ||
     type == "Fou"     ||
     type == "fou") {
                return("fourier")
        }
  else if(type == "mon" ||
          type == "monom"  ||
          type == "monomial") {
                return("monom")
        }
  else if(type == "polyg"    ||
          type == "polygon"  ||
          type == "polygonal") {
                return("polyg")
        }
  else if(type == "poly"    ||
          type == "pol"     ||
          type == "polynom" ||
          type == "polynomial") {
                return("polynom")
        }
  else if(type == "pow"    ||
          type == "power") {
                return("power")
        }
  else {
                return("unknown")
        }
}
