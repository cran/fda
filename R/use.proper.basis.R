use.proper.basis <- function(type) {
  #  recognizes type of basis by use of several variant spellings
  if(type == 'Fourier' ||
     type == 'fourier' ||
     type == 'Fou'     ||
     type == 'fou') {
                return('fourier')
        }
  else if(type == 'bspline' ||
          type == 'Bspline' ||
          type == 'Bsp'     ||
          type == 'bsp') {
                return('bspline')
        }
  else if(type == 'poly' ||
          type == 'pol'  ||
          type == 'polynomial') {
                return('poly')
        }
  else if(type == 'polyg'    ||
          type == 'polygon'  ||
          type == 'polygonal') {
                return('polyg')
        }
  else if(type == 'exp'    ||
          type == 'expon'  ||
          type == 'exponential') {
                return('expon')
        }
  else if(type == 'con'   ||
          type == 'const' ||
          type == 'const') {
                return('const')
        }
  else {
                return('unknown')
        }
}
