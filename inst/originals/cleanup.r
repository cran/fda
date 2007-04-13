#  Set up the data objects used in the examples

#  ------------  gait  --------------------------

hip  <- matrix(scan("hip.txt",  0), 20, 39)
knee <- matrix(scan("knee.txt", 0), 20, 39)

#  set up a three-dimensional array of function values
gait <- array(0, c(20, 39, 2))
dimnames(gait) <- list(NULL, NULL, c("Hip Angle", "Knee Angle"))
gait[,,1] <- hip
gait[,,2] <- knee

save(gait, file="gait.rda", compress=TRUE)

#  ------------  goods index  --------------------------

temp    <- matrix(scan("nondurprod.txt",0), 18, 81)
tempmat <- temp[2:13,]
tempmat[12,81] <- 0
nondurables <- matrix(tempmat, 12*81, 1)
nondurables <- nondurables[1:971]
ndur <- 971

#  for completeness, make dec 99 equal to dec 98, jan 00 equal to jan 99

nondurables <- c(nondurables,nondurables[961])
nondurables <- c(nondurables,nondurables[962])
ndur <- 973

save(nondurables, file="nondurables.rda", compress=TRUE)

#  ------------  Berkeley growth  --------------------------

ncasem <- 39
ncasef <- 54
nage   <- 31

hgtm <- t(matrix(scan("hgtm.txt", 0), ncasem, byrow=TRUE))
hgtf <- t(matrix(scan("hgtf.txt", 0), ncasef, byrow=TRUE))


age <- c( seq(1, 2, 0.25), seq(3, 8, 1), seq(8.5, 18, 0.5))

growth <- list(hgtm=hgtm, hgtf=hgtf, age=age)

save(growth, file="growth.rda", compress=TRUE)


#  ------------  handwriting  --------------------------

temp <- array(scan("fdareg.txt",0), c(20,2,1401))

#  set up a three-dimensional array

handwrit <- array(0, c(1401, 20, 2))
handwrit[,,1] <- t(temp[,1,])/1000
handwrit[,,2] <- t(temp[,2,])/1000
dimnames(handwrit) <- list(NULL, NULL, c("X", "Y") )

save(handwrit, file="handwrit.rda", compress=TRUE)

#  ------------  lip  --------------------------

lip <- matrix(scan("lip.txt", 0), 51, 20)

save(lip, file="lip.rda", compress=TRUE)

#  ------------  melanoma  --------------------------

tempmat <- t(matrix(scan("melanoma.txt", 0), 3, 37))
colnames(tempmat) <- c("index","year","melanoma")

melanoma <- as.data.frame(tempmat)

save(melanoma, file="melanoma.rda", compress=TRUE)


#  ------------  pinch  --------------------------

pinch   <- matrix(scan("pinch.txt",0), 151, 20, byrow=TRUE)
save(pinch, file="pinch.rda", compress=TRUE)

#  ------------  refinery  --------------------------

refinery <- t(matrix(scan("refinery.txt", 0), 3))
#tval <- refinery[,1]  #  observation time
#uval <- refinery[,2]  #  reflux flow
#yval <- refinery[,3]  #  tray 47 level
colnames(refinery) <- c("tval", "uval", "yval")

refinery <- as.data.frame(refinery)
#  center the data on mean values prior to change
refinery <- transform(refinery, 
  yval = yval - mean(yval[1:60]), 
  uval = uval - mean(uval[1:60])
)

save(refinery, file="refinery.rda", compress=TRUE)

#  ------------  daily weather  --------------------------

tempav <- matrix(scan("dailtemp.txt",0), 365, 35)
precav <- matrix(scan("dailprec.txt",0), 365, 35)

#  define 11-character names for stations

place <- c(
"Arvida     ", "Bagottville", "Calgary    ", "Charlottvl ", "Churchill  ", "Dawson     ",
"Edmonton   ", "Fredericton", "Halifax    ", "Inuvik     ", "Iqaluit    ", "Kamloops   ",
"London     ", "Montreal   ", "Ottawa     ", "Pr. Albert ", "Pr. George ", "Pr. Rupert ",
"Quebec     ", "Regina     ", "Resolute   ", "Scheffervll", "Sherbrooke ", "St. Johns  ",
"Sydney     ", "The Pas    ", "Thunderbay ", "Toronto    ", "Uranium Cty", "Vancouver  ",
"Victoria   ", "Whitehorse ", "Winnipeg   ", "Yarmouth   ", "Yellowknife")

dimnames(tempav) <- list(NULL,place)
dimnames(precav) <- list(NULL,place)

#  set up indices that order the stations from east to west to north

geogindex <- c(24,  9, 25, 34,  4,  8, 22,  1,  2, 19, 23, 14, 15, 28, 13, 
               27, 33, 26,  5, 20, 16, 29,  7,  3, 12, 30, 31, 17, 18, 32, 
                6, 35, 11, 10, 21)

#  put the stations in geographical order, from east to west to north
#  rather in the original alphatical order.

library(gdata)
Place <- trim(place)
CanadianWeather <- daily

library(gdata)
(place <- trim(daily$place))

CanadianWeather$place <- place
dimnames(CanadianWeather$tempav)[[2]] <- place
dimnames(CanadianWeather$precav)[[2]] <- place

daysPerMonth <- rep(31, 12)
names(daysPerMonth) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
daysPerMonth[c("Sep", "Apr", "Jun", "Nov")] <- 30
daysPerMonth["Feb"] <- 28

CanadianWeather$month365 <- rep(names(daysPerMonth), daysPerMonth)
CanadianWeather$dayOfYear <- 1:365
Mon. <- with(CanadianWeather,
     tapply(dayOfYear, month365, mean))
oM <- order(Mon.)
CanadianWeather$Month <- Mon.[oM]

CanadianWeather$monthlyTemp <- (with(CanadianWeather,
     apply(tempav, 2, function(x)tapply(x, month365, mean)) )[oM,])

CanadianWeather$monthlyPrecip <- (with(CanadianWeather,
     apply(precav, 2, function(x)tapply(x, month365, mean)) )[oM, ])
CanadianWeather$geogindex <- geogindex
save(CanadianWeather, file="CanadianWeather.rda", compress=TRUE)
str(CanadianWeather)

tempav <- tempav[,geogindex]
precav <- precav[,geogindex]
place  <- place[geogindex]

daily <- list(place=place, tempav=tempav, precav=precav)
save(daily, file="daily.rda", compress=TRUE)
