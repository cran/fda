\name{readHMD}
\alias{readHMD}
\title{
  Download data from the Human Mortality Database (HMD)
}
\description{
  Obtain data from the Human Mortality Database (HMD) maintained at the
  University of California at Berkeley and the Max Planck Institute for
  Demographic Research.  To use \code{readHMD}, you must first obtain a
  \code{username} and \code{password} from HMD and supply them as
  arguments to this function.
}
\usage{
readHMD(username, password,
        country=c(AUS='Australia', AUT='Austria', BLR='Belarus',
                  BEL='Belgium', BGR='Bulgaria', CAN='Canada',
                  CHL='Chile', CZE='Czech Republic', DNK='Denmark',
                  EST='Estonia', FIN='Finland', FRA='France',
                  DEU='Germany', HUN='Hungary', ISL='Iceland',
                  IRL='Ireland', ISR='Israel', ITA='Italy',
                  JPN='Japan', LVA='Latvia', LTU='Lithuania',
                  LUX='Luxemburg', NDL='Netherlands', NZL='New Zealand',
                  NOR='Norway', POL='Poland', PRT='Portugal',
                  RUS='Russia', SVK='Slovakia', SVN='Slovenia',
                  ESP='Spain', SWE='Sweden', CHE='Switzerland',
                  TWN='Taiwan', GBR='U.K.', USA='U.S.A.',
                  UKR='Ukraine'),
        sex=c('m', 'f', 'b'), HMDurl='http://www.mortality.org/hmd',
        dataType = 'lt',
        ltCol=c('m', 'q', 'a', 'l', 'd', 'L', 'T', 'e'),
        cohper = c(coh='cohort', per='periodic'),
        ageInterval=c(1, 5), yearInterval=c(1, 5, 10),
        url, ...)
}
\arguments{
  \item{username, password}{
    user name and password for HMD.
  }
  \item{country}{
    country for which data are required.  Must be one of the available
    options given with the "usage", specified either by the name of the
    3-letter code.
  }
  \item{sex}{
    gender for which data are required (male, female or both), specified
    via a single letter.
  }
  \item{HMDurl}{
    character string giving the URL of HMD.
  }
  \item{dataType}{
    data type:  'lt' for 'life table'.  If a different data type is
    required, you should go to the HMD web site, find the URL you want,
    and specify it directly via the \code{url} argument.
  }
  \item{ltCol}{
    life table column:

    \itemize{
      \item{m}{
	Central death rate between ages x and x+n
      }
      \item{q}{
	Probability of death between ages x and x+n
      }
      \item{a}{
	Average length of survival between ages x and x+n for persons
	dying in the interval
      }
      \item{l}{
	Number of survivors at exact age x, assuming l(0) = 100,000
      }
      \item{d}{
	Number of deaths between ages x and x+n
      }
      \item{L}{
	Number of person-years lived between ages x and x+n
      }
      \item{T}{
	Number of person-years remaining after exact age x
      }
      \item{e}{
	Life expectancy at exact age x (in years)
      }
    }
  }
  \item{cohper}{
    periodic or chorhort data?
  }
  \item{ageInterval}{
    width of the age intervals:  either 1 or 5 years.
  }
  \item{yearInterval}{
    width of the intervals for which results are desired:  either, 1, 5
    or 10 years.
  }
  \item{url}{
    The Universal Resource Locator for the desired table.  If provided,
    all the other arguments except \code{username} and \code{password}
    will be ignored.  If missing, this will be constructed from the
    other arguments.
  }
  \item{...}{
    options arguments for \code{\link[RCurl]{getURL}}
  }
}
\details{
  To use this function, you need a username and password with HMD.  As
  of 2012.07.28, the HMD access is free for individuals, but you must
  register with them.  To start the registration process, go to
  \url{http://www.mortality.org} and click "New User" near the top in
  the left margin.  Then click "New User" again to see the User
  Agreement.

  In all published work and presentations, please acknowledge the HMD as
  either the source or the intermediary of the data with the date on
  which you extracted the data.  See
  \url{http://www.mortality.org/mp/auth.pl} for recommended citation
  formats.

  If you agree to comply with their User Agreement, click "I agree" at
  the bottom of that page.  This will take you to another page where you
  can enter requested information about you including an email address
  to which your password will be mailed.
}
\value{
  If url is provided or dataType != 'lt' or \code{\link{read.table}}
  fails, this returns a list with the following components:

  \item{URL}{
    the URL that was used, being the argument \code{url} if provided or
    the url that was constructed after noting that it was missing.
  }
  \item{getURL}{
    object returned by \code{\link[RCurl]{getURL}}.
  }
  \item{readLines}{
    object returned by readLines(getURL(...)).
  }
  \item{read.table}{
    object returned by read.table(getURL(...)).
  }

  Otherwise, this returns a list with the following components:

  \item{x}{ numeric vector of ages }
  \item{y}{ matrix of the data for each age-year combination }
  \item{time}{ numeric vector of the years }
  \item{xname}{ "age" }
  \item{yname}{
    c(m='Mortality rate',
    q='Mortality probability',
    a='Survival time for mortalities',
    l='Number of survivors',
    d='Number of deaths',
    L='Person-years in interval',
    T='Person-years remaining',
    e='Life expectancy')[ltCol]

    = data type from the life table
  }
  \item{class}{ c('fts', 'fds') }
}
\source{
  Human Mortality Database.  University of California, Berkeley (USA),
  and Max Planck Institute for Demographic Research (Germany). Available
  at \url{www.mortality.org} and \url{www.humanmortality.de}.
}
\author{
  Han Lin Shang, Monash University in Australia, and Spencer Graves
}
\seealso{
  \code{\link[RCurl]{getURL}}
  \code{\link[demography]{hmd.mx}}
}
\examples{
\dontrun{
user <- "your.email@wherever.org" # as registered with HMD.
pw <- "1234567890"
# password as provided in an email from HMD possibly changed by you.

country <- 'Sweden'
SwedeMatAll <- readHMD(user, pw, 'Sweden', sex='f')
}
}
% docclass is function
\keyword{IO}
