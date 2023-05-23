\name{matplot}
\alias{matplot}
\alias{matplot.default}
\alias{matplot.Date}
\alias{matplot.POSIXct}
\title{
  Plot Columns of Matrices
}
\description{
  Plot the columns of one matrix against the columns of another.
}
\usage{
matplot(x, ...)
\method{matplot}{default}(x, y, type = "p", lty = 1:5, lwd = 1,
    lend = par("lend"), pch = NULL, col = 1:6, cex = NULL, bg = NA,
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, ..., add = FALSE,
    verbose = getOption("verbose"))
\method{matplot}{Date}(x, y, type = "p", lty = 1:5, lwd = 1,
    lend = par("lend"), pch = NULL, col = 1:6, cex = NULL, bg = NA,
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, ..., add = FALSE,
    verbose = getOption("verbose"))
\method{matplot}{POSIXct}(x, y, type = "p", lty = 1:5, lwd = 1,
    lend = par("lend"), pch = NULL, col = 1:6, cex = NULL, bg = NA,
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, ..., add = FALSE,
    verbose = getOption("verbose"))
}
\arguments{
  \item{x, y}{
    vectors or matrices of data for plotting.  The number of rows should
    match.  If one of them are missing, the other is taken as 'y' and an
    'x' vector of '1:n' is used.  Missing values ('NA's) are allowed.
  }
  \item{type}{
    character string (length 1 vector) or vector of 1-character strings
    indicating the type of plot for each column of 'y', see 'plot' for
    all possible 'type's.  The first character of 'type' defines the
    first plot, the second character the second, etc.  Characters in
    'type' are cycled through; e.g., '"pl"' alternately plots points and
    lines.
  }
  \item{lty, lwd, lend}{
    vector of line types, widths, and end styles.  The first element is
    for the first column, the second element for the second column,
    etc., even if lines are not plotted for all columns. Line types will
    be used cyclically until all plots are drawn.
  }
  \item{pch}{
    character string or vector of 1-characters or integers for plotting
    characters, see 'points'.  The first character is the
    plotting-character for the first plot, the second for the second,
    etc.  The default is the digits (1 through 9, 0) then the lowercase
    and uppercase letters.
  }
  \item{col}{
    vector of colors.  Colors are used cyclically.
  }
  \item{cex}{
    vector of character expansion sizes, used cyclically.  This works as
    a multiple of 'par("cex")'.  'NULL' is equivalent to '1.0'.
  }
  \item{bg}{
    vector of background (fill) colors for the open plot symbols given
    by 'pch=21:25' as in 'points'.  The default 'NA' corresponds to the
    one of the underlying function 'plot.xy'.
  }
  \item{xlab, ylab}{
    titles for x and y axes, as in 'plot'.
  }
  \item{xlim, ylim}{
    ranges of x and y axes, as in 'plot'.
  }
  \item{\dots}{
    Graphical parameters (see 'par') and any further arguments of
    'plot', typically 'plot.default', may also be supplied as arguments
    to this function.  Hence, the high-level graphics control arguments
    described under 'par' and the arguments to 'title' may be supplied
    to this function.
  }
  \item{add}{
    logical.  If 'TRUE', plots are added to current one, using 'points'
    and 'lines'.
  }
  \item{verbose}{
    logical.  If 'TRUE', write one line of what is done.
  }
}
\details{
  Note that for multivariate data, a suitable array must first be
  defined using the \code{par} function.

  \code{matplot.default} calls \code{\link[graphics]{matplot}}.  The
  other methods are needed, because the default methods ignore the
  \code{Date} or \code{POSIXct} character of \code{x}, labeling the
  horizontal axis as numbers, thereby placing it on the user to
  translate the numbers of days or seconds since the start of the epoch
  into dates (and possibly times for \code{POSIXct} \code{x}).
}
\section{Side Effects}{
  a plot of the functional observations
}
\seealso{
  \code{\link[graphics]{matplot}}, \code{\link{plot}},
  \code{\link{points}}, \code{\link{lines}}, \code{\link{matrix}},
  \code{\link{par}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
##
## matplot.Date, matplot.POSIXct
##
# Date
invasion1 <- as.Date('1775-09-04')
invasion2 <- as.Date('1812-07-12')
earlyUS.Canada <- c(invasion1, invasion2)
Y <- matrix(1:4, 2, 2)
matplot(earlyUS.Canada, Y)

# POSIXct
AmRev.ct <- as.POSIXct1970(c('1776-07-04', '1789-04-30'))
matplot(AmRev.ct, Y)

##
## matplot.default (copied from matplot{graphics})
##
matplot((-4:5)^2, main = "Quadratic") # almost identical to plot(*)
sines <- outer(1:20, 1:4, function(x, y) sin(x / 20 * pi * y))
matplot(sines, pch = 1:4, type = "o", col = rainbow(ncol(sines)))
matplot(sines, type = "b", pch = 21:23, col = 2:5, bg = 2:5,
             main = "matplot(...., pch = 21:23, bg = 2:5)")

x <- 0:50/50
matplot(x, outer(x, 1:8, function(x, k) sin(k*pi * x)),
             ylim = c(-2,2), type = "plobcsSh",
             main= "matplot(,type = \"plobcsSh\" )")
## pch & type =  vector of 1-chars :
matplot(x, outer(x, 1:4, function(x, k) sin(k*pi * x)),
             pch = letters[1:4], type = c("b","p","o"))

lends <- c("round","butt","square")
matplot(matrix(1:12, 4), type="c", lty=1, lwd=10, lend=lends)
text(cbind(2.5, 2*c(1,3,5)-.4), lends, col= 1:3, cex = 1.5)

table(iris$Species) # is data.frame with 'Species' factor
iS <- iris$Species == "setosa"
iV <- iris$Species == "versicolor"
op <- par(bg = "bisque")
matplot(c(1, 8), c(0, 4.5), type= "n", xlab = "Length", ylab = "Width",
             main = "Petal and Sepal Dimensions in Iris Blossoms")
matpoints(iris[iS,c(1,3)], iris[iS,c(2,4)], pch = "sS", col = c(2,4))
matpoints(iris[iV,c(1,3)], iris[iV,c(2,4)], pch = "vV", col = c(2,4))
legend(1, 4, c("    Setosa Petals", "    Setosa Sepals",
                    "Versicolor Petals", "Versicolor Sepals"),
            pch = "sSvV", col = rep(c(2,4), 2))

nam.var <- colnames(iris)[-5]
nam.spec <- as.character(iris[1+50*0:2, "Species"])
iris.S <- array(NA, dim = c(50,4,3),
                     dimnames = list(NULL, nam.var, nam.spec))
for(i in 1:3) iris.S[,,i] <- data.matrix(iris[1:50+50*(i-1), -5])

matplot(iris.S[,"Petal.Length",], iris.S[,"Petal.Width",], pch="SCV",
             col = rainbow(3, start = .8, end = .1),
             sub = paste(c("S", "C", "V"), dimnames(iris.S)[[3]],
                         sep = "=", collapse= ",  "),
             main = "Fisher's Iris Data")
par(op)
par(oldpar)
}
\keyword{hplot}
