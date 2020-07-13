\name{freq_genotypes}
\alias{freq_genotypes}
\alias{freq_alleles}
\alias{het_obs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotype and allele frequencies, heterozygosities 
}
\description{
Functions to calculate the frequencies of single-locus genotypes, single-locus alleles and heterozygosities from a vector of numbers of individuals at each multilocus genotype.
}
\usage{
freq_genotypes(N, init.par)
freq_alleles(N, init.par, fG = NULL)
het_obs(N, init.par, fG = NULL)
het_exp(N, init.par, fA = NULL)
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
  \item{fA}{
The allele frequencies that will be used to the calculate the expected heterozygosity: the output of \code{freq_alleles}. If \code{NULL}, the function will call \code{freq_alleles}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
For \code{freq_genotypes}, a list with two components. The first component (\code{$counts}) is itself a list giving the genotype counts for each single-locus genotype at each locus. The second component (\code{$frequencies}) is like the first one, but giving frequencies instead of counts.\cr\cr
For \code{freq_alleles}: as for \code{freq_genotypes}, but giving allele counts and allele frequencies instead of genotype counts and genotype frequencies.\cr\cr
For \code{het_obs}, a vector of observed heterozygosities (ratio of number of heterozygotes to total number of individuals) per locus.\cr\cr
For \code{het_exp}, a vector of expected heterozygosities per locus, calculated as 1 - sum(p^2), where p is the allele frequency and the sum is taken over all alleles at a locus.
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
fG <- freq_genotypes(N_deme1_age1_time3, init.par)
fG$counts
fG$frequencies

# Allele frequencies
fA <- freq_alleles(N_deme1_age1_time3, init.par)
# Equivalent but a bit faster:
fA <- freq_alleles(N_deme1_age1_time3, init.par, fG)

# Observed heterozygosity
Ho <- het_obs(N_deme1_age1_time3, init.par, fG)

# Expected heterozygosity
He <- het_exp(N_deme1_age1_time3, init.par, fA)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
