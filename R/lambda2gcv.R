lambda2gcv <- function(log10lambda, argvals, y, fdParobj,
                       wtvec=rep(1,length(y)))
{
#  LAMBDA2GCV smooths data using smooth_basis and returns the
#  GCV criterion that results.  DFFACTOR is an optional multiplier of df
#  in the definition of GCV

fdParobj$lambda = 10^log10lambda

#smoothlist <- smooth.basis(argvals, y, fdParobj, wtvec, dffactor)
smoothlist <- smooth.basis(argvals, y, fdParobj, wtvec)

gcv <- smoothlist$gcv

gcv

}
