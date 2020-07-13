\name{freq_genotypes}
\alias{freq_genotypes}
\alias{freq_alleles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotype and allele frequencies 
}
\description{
Functions to calculate the frequencies of single-locus genotypes and single-locus alleles from a vector of numbers of individuals at each multilocus genotype.
}
\usage{
freq_genotypes(N, init.par)
freq_alleles(N, init.par, fG=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Vector (1-D array) giving the number of individuals of each multilocus genotype. Typically this is obtained from the output of a multilocus simulation by subsetting a single deme / age class / time step. See the example.
}
  \item{init.par}{
A list of arguments necessary to initialize the simulations and obtained with the function initialize.monoecious.multilocus. It must be the same list that has been used to simulate \code{N}
}
  \item{fG}{
The genotype frequencies that will be used to the calculate allele frequencies: the output of \code{freq_genotypes}. If \code{NULL}, the function will call \code{freq_genotypes}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Lists with two components. For \code{freq_genotypes}, the first component (\code{$counts}) is itself a list giving the genotype counts for each single-locus genotype at each locus. The second component (\code{$frequencies}) is like the first one, but giving frequencies instead of counts. For \code{freq_alleles}, the first component gives the allele counts and the second component gives the allele frequencies
}
\references{
}
\author{Marco Andrello}
\note{
}
\seealso{
\code{\link{freq.all}} to calculate allele frequencies in single-locus simulations
}
\examples{
\dontrun{
data(MultilocusMonoecious)
N <- sim.metapopgen.monoecious.multilocus(init.par=init.par,
                                          sigma=sigma,
                                          phi_F=phi_F, phi_M=phi_M,
                                          delta.prop=delta.prop, delta.ad=delta.ad,
                                          T_max=T_max)
# Subsetting to get results from deme 1, age class 1 at time 3
N_deme1_age1_time3 <- N[,1,1,3]

# Genotype frequencies
f_genotypes <- freq_genotypes(N_deme1_age1_time3, init.par)
f_genotypes$counts
f_genotypes$frequencies

# Allele frequencies
fA <- freq_alleles(N_deme1_age1_time3, init.par)
# Equivalent but a bit faster:
fA <- freq_alleles(N_deme1_age1_time3, init.par, fG)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
