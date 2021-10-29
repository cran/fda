\name{NEWS}
\alias{NEWS}
\title{
Recent changes to the fda package
}
\description{
    Changes in version fda_5.5.0 2021-10-28:
    
    Many data smoothing situations require that the smooth
    curves satisfy some constraints. 
    
    Take function \code{smooth.monotone.R} for example. Its curves are either strictly increasing or strictly decreasing, even though the data are not.  This is the case in modelling human growth, where we can reasonably assume that daily or monthly measurements will reflect a trend that increases everywhere.  
    
    Function \code{smooth.morph.R}, which plays an important role in curve registration, adds the additional constraint that the domain limits mapped exactly into the range limits.
    
    In this version two new constrained curves are introduced.  Nonsingular multinomial probability vectors contain nonzero probabilities that sum to zero.  A simple transformation of these probabilities, $S = -log(P)$, converts probabilities into what is often called surprisal.  Surprisal is a measure of information where the unit of measurement is the M-bit, where $M$ is the length of the multinomial vector.  Information measured in this way can be added and subtracted, and fixed differences mean the same thing anywhere along the surprisal continuum, which is positive with an origin at 0.  Probability 1 corresponds to surprisal 0, and a very small probability produces a very large positive surprisal.  Probabilities 0.05 and 0.01 correspond to 2-bit surprisals 4.3 and 6.1, respectively.
    
    Probability curves result if the probabilities change with over continuous scale, often called a latent variable in statistics.  The corresponding surprisal curves satisfy the constraint at any index value $log(sum(M^S)) = 0.$  The unbounded nature of surprisal curves plus their metric property render them much easier to work with computationally, as well having the metric property.
    
    Functions smooth.surp.R and error sum of squares fit function surp.fit.R are added in this version in order to support a package \code{TestGardener} that analyzes choice or psychometric data.
    
    Function \code{smooth.morph.R} is also now extended by function \code{smooth.morph2.R} in order to map the limits of a domain into different limits for the range.
    
}

