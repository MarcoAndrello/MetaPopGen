\name{genotype.index}
\alias{genotype.index.multilocus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotype indices
}
\description{
Functions giving genotype indices as a function of the number of alleles (for a single locus) or the number of allele at each locus (for multiple loci).
}
\usage{
genotype.index(l)
genotype.index.multilocus(allele_vec,r,method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{
Number of alleles
}
  \item{allele_vec}{
Vector giving the number of alleles at each locus 
}
  \item{r}{
Recombination rate. A single value between 0 and 0.5. For more than two loci, recombination rate must be set to 0.5.
}
  \item{method}{
A string defining the method used to store the genotypes. One of either \code{"gamete"} or \code{"genotype"}. For two loci and \code{r} < 0.5, method must be \code{"gamete"}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
For \code{genotype.index}, a square matrix whose [i,j] element gives the index of the genotype formed by the i and j alleles.\cr
\cr
For \code{genotype.index.multilocus}, the value depends on the argument \code{method}. If \code{method = "gamete"}, it is a square matrix whose [i,j] element gives the index of the multilocus genotype formed by the i and j gametes. If \code{method = "genotype"}, it is an n-dimensional array, where n is the number of loci, giving the index of the multilocus genotype formed by the ith single-locus genotype at the first locus, the jth single-locus genotype at the second locus, etc.
}
\references{
}
\author{
Marco Andrello
}
\note{
}
\seealso{
}
\examples{
# Monolocus case:
genotype.index(2)

# Multilocus case:
allele_vec <- c(2,2)
genotype.index.multilocus(allele_vec,r=0.5,method="gamete")
genotype.index.multilocus(allele_vec,r=0.5,method="genotype")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
