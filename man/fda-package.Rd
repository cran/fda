\docType{package}
\name{fda}
\alias{fda}
\title{Functions for statistical analyses of functions}

\description{
The data analyses that we use for data matrices or dataFrames, such as means, covariances, linear regressions, principal and canonical component analyses, can also be applied to samples of functions or curves.  This package provides functional versions of these analyses along with plotting and other assessment tools.  

The functions themselves are often a consequence of smoothing discrete data values over domains like time, space, and other continuous variables.  Data smoothing tools are also provided.

But there are transformations of functions that have no meaning for data matrices for which rows may be re-ordered.  Derivatives and integrals are often used to set up dynamic models, and can also play a constructive role in the data smoothing process.  Methods for these and other functional operations are also available in this package.

There are now many texts and papers on functional data analysis.  The two resources provided by the buildeers of the fda package are:

James Ramsay and Bernard Silverman (2005) Functional Data Analysis. New York, Springer.

James Ramsay, Giles Hooker and Spencer Graves (2009) Functional Data Analysis with R and Matlab.  New York: Springer.

Another relevant package is the package Data2LD that offers ann introduction and methods for constructing linear and nonlinear differential equation models, along with the text: 

James Ramsay and Giles Hooker, Dynamic Data Analysis, New York: Springer.
}