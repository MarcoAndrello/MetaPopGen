\name{het_obs}
\alias{het_obs}
\alias{het_exp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heterozygosities 
}
\description{
Functions to calculate the observed and expected heterozygosities per locus. The input is the vector of numbers of individuals at each multilocus genotype.
}
\usage{
het_obs(N, init.par, fG = NULL)
het_exp(N, init.par, fA = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Vector (1-D array) giving the number of individuals of each multilocus genotype. Typically this is obtained from the output of a multilocus simulation by subsetting a single deme / age class / time step. See the example.
}
  \item{init.par}{
A list of arguments necessary to initialize the simulations and obtained with the function initialize.monoecious.multilocus. It must be the same list that has been used to simulate \code{N}.
}
  \item{fG}{
The genotype frequencies that will be used to the calculate allele frequencies: the output of \code{\link{freq_genotypes}}. If \code{NULL}, the function will call \code{\link{freq_genotypes}}.
}
  \item{fA}{
The allele frequencies that will be used to the calculate the expected heterozygosity: the output of \code{\link{freq_alleles}}. If \code{NULL}, the function will call \code{\link{freq_alleles}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
For \code{het_obs}, a vector of observed heterozygosities (ratio of number of heterozygotes to total number of individuals) per locus.\cr\cr
For \code{het_exp}, a vector of expected heterozygosities per locus, calculated as 1 - sum(p^2), where p is the allele frequency and the sum is taken over all alleles at a locus.
}
\references{
}
\author{Marco Andrello}
\note{
}
\seealso{
\code{\link{freq_genotypes}}
\code{\link{freq_alleles}}
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

# Observed heterozygosity
Ho <- het_obs(N_deme1_age1_time3, init.par)

# Expected heterozygosity
He <- het_exp(N_deme1_age1_time3, init.par)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
