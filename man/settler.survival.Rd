\name{settler.survival}
\alias{settler.survival}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Settler-recruit survival function
}
\description{
}
\usage{
settler.survival(S,a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
Number of settlers (scalar)
}
  \item{a}{
Parameter defining the strength of density-dependence (scalar)
}
}
\details{
This function is not needed by the user. It is used by the sim.metapopgen function to calculate density-dependent recruitment. The density-dependent function is \code{y = (1 / (1 + a * S )}
}
\value{
A probability of survival 
}
\references{
}
\author{
}
\note{
}
\seealso{
}
\examples{

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
