\name{fst}
\alias{fst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic function to calculate Fst 
}
\description{
  Fst for a single locus
}
\usage{
fst(N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
A 2-D array. It gives the number of individuals of each genotype (first dimension) in each group (second dimension)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function calculates the genetic differentiation index Fst among the groups contained in N as ratio of heterozygosities: \code{(Ht - Hs) / Ht}. \code{Ht} is the expected heterozygosity over the population formed by pooling the groups together. \code{Hs} is the average expected heterozygosity over the groups. The average is weighted by the number of individuals in each group.
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
# Data: three groups with two alleles
# Genotypic composition of first, second and third groups
N1 <- c(18,299,277)
N2 <- c(284,14,145)
N3 <- c(115,228,254)
N <- array(c(N1,N2,N3),dim=c(3,3),dimnames=list(genotype=c("A1A1","A1A2","A2A2"),group=c(1:3)))
fst(N)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
