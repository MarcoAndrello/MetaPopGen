\name{freq.all}
\alias{freq.all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating allele frequencies 
}
\description{
This function calculates allele frequencies from a vector of numbers of individuals of each genotype at one locus.
}
\usage{
freq.all(N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Vector (1-D array) giving the number of individuals of each genotype, ordered in the standard way (see description in \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}})
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of allele frequencies. 
}
\references{
}
\author{Marco Andrello}
\note{
}
\seealso{
}
\examples{
N <- c(30,50,20) # Means that there are 30 A1A1, 50 A1A2 AND 20 A2A2
freq.all(N)       # The frequency of A1 is 0.55 and the frequency of A2 is 0.45 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
