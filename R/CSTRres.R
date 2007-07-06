CSTRres <- function(kref=NULL, EoverR=NULL, a=NULL, b=NULL,
               datstruct, fitstruct, CSTRbasis, lambda,
               gradwrd=FALSE){
#%  2007.06.02 by Spencer Graves
##
## 1.  Construct 'parvec'
##
#  parvec <- vector("list", 4)
#  names(parvec) <- c("kref", "EoverR", "a", "b")
#  parvec$kref <- kref
#  parvec$Eover <- EoverR
#  parvec$a <- a
#  parvec$b <- b
  parvec <- list(kref=kref, EoverR=EoverR, a=a, b=b)
  pv <- unlist(parvec)
##
## 2.  Call CSTRfn
##
  cstr. <- CSTRfn(parvec=pv, datstruct=datstruct, fitstruct=fitstruct,
         CSTRbasis=CSTRbasis, lambda=lambda, gradwrd=gradwrd)
  Res <- as.vector(cstr.$res)
  if(length(d.r <- dim(Res))>1){
    ys <- dimnames(Res)[[2]]
    if(!is.null(ys)){
      resNames <- t(outer(ys, 1:d.r[1], paste, sep=""))
      Res <- as.vector(Res)
      names(Res) <- resNames
    }
  }
  if(gradwrd)
    attr(Res, "gradient") <- cstr.$Dres
#  
  Res
}

CSTRres0 <- function(kref=NULL, EoverR=NULL, a=NULL, b=NULL,
                     gradwrd=FALSE){
#%  2007.06.02 by Spencer Graves
#cat("CSTRres0: kref = ",kref, "; EoverR = ",EoverR,"; a = ",a, "; b = ",b,"\n")
##
## 1.  Construct 'parvec'
##
#  parvec <- vector("list", 4)
#  names(parvec) <- c("kref", "EoverR", "a", "b")
#  parvec$kref <- kref
#  parvec$Eover <- EoverR
#  parvec$a <- a
#  parvec$b <- b
  parvec <- list(kref=kref, EoverR=EoverR, a=a, b=b)
  pv <- unlist(parvec)
##
## 2.  'get' the other arguments of CSTRfn
##
  datstr <- get(".datstruct")
  fitstr <- get(".fitstruct")
  CSTRb <- get(".CSTRbasis")
  lam <- get(".lambda")
##
## 3.  Call CSTRfn
##
  cstr. <- CSTRfn(parvec=pv, datstruct=datstr, fitstruct=fitstr,
         CSTRbasis=CSTRb, lambda=lam, gradwrd=gradwrd)
  Res <- as.vector(cstr.$res)
  if(length(d.r <- dim(Res))>1){
    ys <- dimnames(Res)[[2]]
    if(!is.null(ys)){
      resNames <- t(outer(ys, 1:d.r[1], paste, sep=""))
      Res <- as.vector(Res)
      names(Res) <- resNames
    }
  }
  if(gradwrd)
    attr(Res, "gradient") <- cstr.$Dres
#  
  Res
}
