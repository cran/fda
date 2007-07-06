basisobj1 = create.bspline.basis(c(0,1),5)
basisobj2 = create.bspline.basis(c(0,1),6)
coef1 = matrix(rnorm(5),5,1))
coef2 = matrix(rnorm(6),6,1))
fdnames = list(3,1)
fdnames[[1]] = "Time"
fdnames[[2]] = "Reps"
fdnames1 = fdnames
fdnames1[[3]] = "f1"
fdobj1 = fd(coef1, basisobj1, fdnames1)
fdnames2 = fdnames
fdnames2[[3]] = "f2"
fdobj2 = fd(coef2, basisobj2, fdnames2)

figure(1)
subplot(2,1,1)
plot(fdobj1)
subplot(2,1,2)
plot(fdobj2)

fdprodobj = fdobj1*fdobj2

figure(2)
subplot(1,1,1)
plot(fdprodobj)

prodbasisobj = getbasis(fdprodobj)
prodbasisobj

figure(3)
plot(getbasis(fdprodobj))

coef1  = matrix(rnorm(5,N)
fdobj1 = fd(coef1, basisobj1, fdnames1)

fac = 2*ones(10,1)
fac = matrix(rnorm(10,1)

fdprodobj = fdobj1.*fac

figure(1)
subplot(1,1,1)
plot(fdobj1)

figure(2)
subplot(1,1,1)
plot(fdprodobj)

N = 20

conbas = create.constant.basis(c(0,1))
confd  = fd(matrix(1,1,N),conbas)

basisobj2 = create.bspline.basis(c(0,1),5)
basisobj3 = create.bspline.basis(c(0,1),6)

coef2 = matrix(rnorm(5*N),5,N)
xobj2 = fd(coef2, basisobj2, fdnames1)
coef3 = matrix(rnorm(6*N),6,N)
xobj3 = fd(coef3, basisobj3, fdnames2)

yfd = confd + 2*xobj2 - 2*xobj3

ybasis = getbasis(yfd)

plot(ybasis)

xfdlist = list(1,3)
xfdlist{1} = confd
xfdlist{2} = xobj2
xfdlist{3} = xobj3

betalist = list(1,3)

betalist{1} = fdPar(conbas)
betalist{2} = fdPar(conbas)
betalist{3} = fdPar(conbas)

yfdPar = fdPar(yfd)

fRegressCell = fRegress(yfdPar, xfdlist, betalist)

betaestlist = fRegressCell{4}

for i=1:3
    subplot(3,1,i)
    plot(getfd(betaestlist{i}))
end

betabasis = create.bspline.basis(c(0,1),4)
betafd1 = fd(1, conbas)
betalist{1} = fdPar(betafd1)
for j=2:3
    betafdj = fd(matrix(rnorm(4,1),betabasis)
    betalist{j} = fdPar(betafdj)
end

yfd0 = getfd(betalist{1}) + ...
       getfd(betalist{2}).*xfdlist{2} + ...
       getfd(betalist{3}).*xfdlist{3}
ybasis = getbasis(yfd0)
ynbasis = getnbasis(ybasis)

sigma = 0.1
efd = fd(matrix(rnorm(ynbasis,N).*sigma, ybasis)
yfd = yfd0 + efd

yfdPar = fdPar(yfd)

fRegressCell = fRegress(yfdPar, xfdlist, betalist)

betaestlist = fRegressCell{4}

tfine = linspace(0,1,101)"
for j=1:3
    subplot(3,1,j)
    plot(tfine, eval.fd(tfine,getfd(betaestlist{j})), "-", ...
         tfine, eval.fd(tfine,getfd(   betalist{j})), "--")
end
