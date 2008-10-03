#  Sets up the daily weather data as a named list containing:
#  weather$temp:  a 365 by 35 matrix containing average temperature 
#         in degrees Celsius for
#  weather$prec:  a 365 by 35 matrix containing average precipitation 
#         in millimetres for
#  weather$time:  a vector of 365 mid-day times in days
#  weather$names: a vector of station names

#  The list is saved with file name "weatherdata" and
#  can be loaded with the command
#  load("weatherdata")

#  These commands are executed in folder fdaR/demo

#  Last modified 28 February 2007

#  ------------------------  input the data  -----------------------

tempav <- matrix(scan("../data/dailtemp.txt",0), 365, 35)
precav <- matrix(scan("../data/dailprec.txt",0), 365, 35)

#  set up the times of observation at noon

daytime   <- (1:365)-0.5

#  define 11-character names for stations

station <- c(
"Arvida     ", "Bagottville", "Calgary    ", "Charlottvl ", "Churchill  ", "Dawson     ",
"Edmonton   ", "Fredericton", "Halifax    ", "Inuvik     ", "Iqaluit    ", "Kamloops   ",
"London     ", "Montreal   ", "Ottawa     ", "Pr. Albert ", "Pr. George ", "Pr. Rupert ",
"Quebec     ", "Regina     ", "Resolute   ", "Scheffervll", "Sherbrooke ", "St. Johns  ",
"Sydney     ", "The Pas    ", "Thunderbay ", "Toronto    ", "Uranium Cty", "Vancouver  ",
"Victoria   ", "Whitehorse ", "Winnipeg   ", "Yarmouth   ", "Yellowknife")

#  set up indices that order the stations from east to west to north

geogindex <- c(24,  9, 25, 34,  4,  8, 22,  1,  2, 19, 23, 14, 15, 28, 13,
               27, 33, 26,  5, 20, 16, 29,  7,  3, 12, 30, 31, 17, 18, 32,
                6, 35, 11, 10, 21)

#  put the stations in geographical order, from east to west to north
#  rather in the original alphatical order.

tempav  <- tempav[,geogindex]
precav  <- precav[,geogindex]
station <- station[geogindex]

weatherdata <- list(tempav  = tempav,  precav = precav,
                    daytime = daytime, station = station)

save(weatherdata, file = "weatherdata")

