
#  all 170688 half-hourly demand levels (3556) days

#  ymat is a 48 by 3556 matrix

ymat = SAelectdemand$y

#  plot two days of demand

plot(1:96, ymat[1:96], type="p")

#  plot the first week, Sunday through Saturday

matplot(1:48, ymat[,1:7], type="b")

nweek = 12
n2 = 0
par(ask=TRUE)
for (week in 1:nweek) {
  n1 = n2 + 1
  n2 = n2 + 7
  matplot(1:48, ymat[,n1:n2], type="b", xlim=c(0,48), ylim=c(800,1900), main=paste("week",week))
  lines(c(24,24),c(800,1900), lty=2)
}

tobs = matrix((1:48) - 0.5,48,1)

#  data smoothing

rng = c(0,48)
nbasis = 51
ebasis = create.bspline.basis(rng, nbasis)

lambda = 1e-4
efdPar = fdPar(ebasis, 2, lambda)

yfd = smooth.basis(tobs, ymat, efdPar)

loglam = seq(-4,0,0.5)
nloglam = length(loglam)
gcv = rep(0,nloglam)
df  = rep(0,nloglam)
for (ilam in 1:nloglam) {
  lambda = 10^loglam[ilam]
  efdPar = fdPar(ebasis, 2, lambda)
  fdsmooth  = smooth.basis(tobs, ymat, efdPar)
  df[ilam]  = fdsmooth$df
  gcv[ilam] = mean(fdsmooth$gcv)
}

# loglam   df   gcv
# [1,]   -4.0 47.2 465.5
# [2,]   -3.5 46.3 330.0
# [3,]   -3.0 44.7 221.4
# [4,]   -2.5 42.1 168.7
# [5,]   -2.0 38.4 151.3
# [6,]   -1.5 33.4 153.7
# [7,]   -1.0 27.8 178.4
# [8,]   -0.5 22.3 240.2
# [9,]    0.0 17.5 365.2

lambda = 1e-2
efdPar = fdPar(ebasis, 2, lambda)
yfd = smooth.basis(tobs, ymat, efdPar)$fd

plotfit.fd(ymat, tobs, yfd)

xmat = tempairport$y

xmatMon = mondaytempairport$y
xmatTue = tuesdaytempairport$y
xmatWed = wednesdaytempairport$y
xmatThu = thursdaytempairport$y
xmatFri = fridaytempairport$y
xmatSat = saturdaytempairport$y
xmatSun = sundaytempairport$y

#  best to analyze each day separately over 508 weeks

#  smooth the Sunday demand data

loglam = seq(-4,0,0.5)
nloglam = length(loglam)
gcv = rep(0,nloglam)
df  = rep(0,nloglam)
for (ilam in 1:nloglam) {
  lambda = 10^loglam[ilam]
  efdPar = fdPar(ebasis, 2, lambda)
  fdsmooth  = smooth.basis(tobs, sundaydemand$y, efdPar)
  df[ilam]  = fdsmooth$df
  gcv[ilam] = mean(fdsmooth$gcv)
}

round(cbind(loglam, df, gcv),1)

# loglam   df   gcv
# [1,]   -4.0 47.2 402.3
# [2,]   -3.5 46.3 278.8
# [3,]   -3.0 44.7 183.8
# [4,]   -2.5 42.1 142.8
# [5,]   -2.0 38.4 129.4  **
# [6,]   -1.5 33.4 130.4
# [7,]   -1.0 27.8 152.4
# [8,]   -0.5 22.3 210.4
# [9,]    0.0 17.5 324.9

lambda = 10^loglam[5]
efdPar = fdPar(ebasis, 2, lambda)
sunelectfd = smooth.basis(tobs, sundaydemand$y, efdPar)$fd
  
plotfit.fd(sundaydemand$y, tobs, sunelectfd)

#  smooth the Sunday temperature data

loglam = seq(-2,2,0.5)
nloglam = length(loglam)
gcv = rep(0,nloglam)
df  = rep(0,nloglam)
for (ilam in 1:nloglam) {
  lambda = 10^loglam[ilam]
  efdPar = fdPar(ebasis, 2, lambda)
  fdsmooth  = smooth.basis(tobs, sundaytempairport$y, efdPar)
  df[ilam]  = fdsmooth$df
  gcv[ilam] = mean(fdsmooth$gcv)
}

round(cbind(loglam, df, gcv),2)

# loglam    df  gcv
# [1,]   -2.0 38.39 0.58
# [2,]   -1.5 33.45 0.44
# [3,]   -1.0 27.83 0.37
# [4,]   -0.5 22.31 0.35 **
# [5,]    0.0 17.47 0.36
# [6,]    0.5 13.54 0.39
# [7,]    1.0 10.47 0.46
# [8,]    1.5  8.13 0.58
# [9,]    2.0  6.35 0.78

lambda = 10^loglam[4]
efdPar = fdPar(ebasis, 2, lambda)
suntempfd = smooth.basis(tobs, sundaytempairport$y, efdPar)$fd

plotfit.fd(sundaytempairport$y, tobs, suntempfd)

#  severe outliers and missing data cases: drop these:

c(2,6,64,65,66)

suntempmat = sundaytempairport$y
suntempmat = suntempmat[,-c(2,6,64,65,66)]

sunelectmat = sundaydemand$y
sunelectmat = sunelectmat[,-c(2,6,64,65,66)]

#  there are now 503 sunday records

lambda = 10^loglam[4]
efdPar = fdPar(ebasis, 2, lambda)
suntempfd = smooth.basis(tobs, suntempmat, efdPar)$fd

lambda = 10^loglam[5]
efdPar = fdPar(ebasis, 2, lambda)
sunelectfd = smooth.basis(tobs, sunelectmat, efdPar)$fd

#  beta basis

rng = c(0,48)
nbasis = 51
bbasis = create.bspline.basis(rng, nbasis)

lambda = 10^-0.5
bfdPar = fdPar(bbasis, 2, lambda)

cbasis  = create.constant.basis(rng)
constfd = fd(matrix(1,1,503),cbasis)

xfdlist = vector("list",2)
xfdlist[[1]] = constfd
xfdlist[[2]] = suntempfd

betalist = vector("list",2)
betalist[[1]] = bfdPar
betalist[[2]] = bfdPar

fRegressResult = fRegress(sunelectfd, xfdlist, betalist)

betaestfd1 = fRegressResult$betaestlist[[1]]
betaestfd2 = fRegressResult$betaestlist[[2]]

par(mfrow=c(2,1), ask=FALSE)
plot(betaestfd1)
title("Intercept")
plot(betaestfd2)
title("Regression on temperature")

yhatfd = fRegressResult$yhatfd$fd

tfine = seq(0,48,len=201)
par(mfrow=c(1,1), ask=TRUE)
for (i in 1:503) {
  yveci = eval.fd(tfine, sunelectfd[i])
  yhati = eval.fd(tfine, yhatfd[i])
  matplot(tfine, cbind(yveci,yhati), type="l", lty=1, col=c(1,2), 
          ylim=c(800,1800), main=paste("Sunday",i))
}
