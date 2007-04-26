quadset <- function(nquad=5, basisobj=NULL, breaks){
#function [basisobj, quadpts, quadwts] = quadset(nquad, basisobj, breaks)

# last modified 2007 April 25 by Spencer Graves
#%  Matlab version last modified 14 August 2006
##
## 1.  Check nquad
##  
  {
    if(nquad<5){
      warning("nquad must be at least 5;  increase to this minimum.")
      if((nquad%%2)!=1){
        warning("nquad must be an odd integer;  increased to enforce this.")
        nquad <- 1+2*ceil(nquad/2)
      }
    }
  }
##  
## 2.  check basisobj
##  
  if(!is.null(basisobj) && !is.basis(basisobj))
    stop('basisobj is not a basis object.')
##  
## 3.  check breaks
##
  if(missing(breaks)) {
    if(is.null(basisobj) || !is.basis(basisobj))
      stop("Either 'breaks' or 'basisobj' must be provided.")
#
#   type     = getbasistype(basisobj);
    type <- basisobj$type 
    if(type != 'bspline')
      stop(
        "'breaks' not supplied and 'basisobj' is not a spline basis.")
#
#    rangeval = getbasisrange(basisobj);
    rangeval = basisobj$range #????? 
#    params   = getbasispar(basisobj);
    params   = basisobj$params
#    knots    = [rangeval(1), params, rangeval(2)];
    knots <- c(rangeval[1], params, rangeval[2]) 
    breaks   = unique(knots)
  }
##
## 4.  quadpts and quadwts
##
  nbreaks = length(breaks);
#
  db <- diff(breaks)
  nquad1 <- nquad-1
  nbreaks1 <- nbreaks-1
  quadpts. <- array(NA, dim=c(nbreaks1, nquad))
  quadpts.[, 1] <- breaks[-nbreaks]
  db. <- db/nquad1
  for(i in 2:nquad)
    quadpts.[, i] <- (quadpts.[, i-1]+db.)
  quadpts <- as.vector(t(quadpts.))
#
  quadwts. <- outer(c(1, 4, rep(c(2, 4), (nquad1-2)/2), 1),
                   db/(nquad1*3) )
  quadwts <- as.vector(quadwts.)
#  quadpts1 = linspace(breaks(1),breaks(2),nquad)';
# quadwts1 = ones(nquad,1);
# quadwts1(2:2:nquad-1) = 4;
# quadwts1(3:2:nquad-2) = 2;
# quadwts1 = ((breaks(2)-breaks(1))/(nquad-1)).*quadwts1/3;
# quadvals = [quadpts1, quadwts1];
# for( i in 3:nbreaks )
#    quadptsi = linspace(breaks(i-1),breaks(i),nquad)';
#   quadwtsi = ones(nquad,1);
#   quadwtsi(3:2:nquad-2) = 2;
#   quadwtsi(2:2:nquad-1) = 4;
#    quadwtsi = ((breaks(i)-breaks(i-1))/(nquad-1)).*quadwtsi/3;
#    quadvals = [quadvals;[quadptsi, quadwtsi]];
#end
  quadvals <- cbind(quadpts=quadpts, quadwts=quadwts)
#quadpts = quadvals(:,1);
#quadwts = quadvals(:,2);

#basisobj = putquadvals(basisobj, quadvals);
  if(is.null(basisobj))return(quadvals)
#  
  basisobj$quadvals <- quadvals
  values <- vector('list', 2) 
  for( ivalue in 1:2){
#    basisvalues    = eval_basis(quadpts, basisobj, ivalue-1);
#    values{ivalue} = basisvalues;
    values[[ivalue]] <- eval.basis(quadpts, basisobj, ivalue-1)
  }
# basisobj = putvalues(basisobj, values);
  basisobj$values <- values
  basisobj
}

