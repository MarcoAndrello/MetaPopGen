\name{freq_genotypes}
\alias{freq_genotypes}
\alias{freq_alleles}
\alias{freq_gametes}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotype, allele and gamete frequencies
}
\description{
Functions to calculate the frequencies of single-locus genotypes, the frequencies of single-locus alleles and the frequencies of multi-locus gametes. The input is the vector of numbers of individuals at each multilocus genotype.
}
\usage{
freq_genotypes(N, init.par)
freq_alleles(N, init.par, fG = NULL)
freq_gametes(N, init.par)
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
The genotype frequencies that will be used to the calculate allele frequencies: the output of \code{freq_genotypes}. If \code{NULL}, the function will call \code{freq_genotypes}.
}
}
\details{
When the recombination rate \code{r} = 0.5, the loci are independent and MetaPopGen uses a memory-efficient storing method for genotypes that does not consider phasing ("locus-based" method, see Andrello et al. under review). For example, the two genotypes A1B1/A2B2 and A1B2/A2B1 are pooled into the same genotype A1B1/A2B2 because they are equivalent in terms of multi-locus gamete production. The computation of gamete frequencies is affected by this loss of information, as the A1B1/A2B2 genotype could be formed i) by gametes A1B1 and A2B2 (with probability 0.5) or ii) by gametes A1B2 and A2B1 (with probability 0.5). The gamete counts obtained by \code{freq_gametes} are calculated considering these probabilities, thus might not be integers. Conversely, when \code{r} < 0.5 (linked loci), the loci are non-independent (linked) and MetaPopGen uses a storing method for genotypes that considers phasing ("gamete-based" method). In this case, the gamete counts obtained by \code{freq_gametes} are calculated by counting the gametes forming the phased genotypes and are always integers.  
}
\value{
For \code{freq_genotypes}, a list with two components. The first component (\code{$counts}) is itself a list giving the genotype counts for each single-locus genotype at each locus. The second component (\code{$frequencies}) is like the first one, but giving frequencies instead of counts.\cr\cr
For \code{freq_alleles}: as for \code{freq_genotypes}, but giving allele counts and allele frequencies instead of genotype counts and genotype frequencies.\cr\cr
For \code{freq_gametes}, a list with two components. The first component (\code{$counts}) is a vector of gamete counts and the second component (\code{$frequencies}) is a vector of gamete frequencies.\cr\cr
}
\references{
Andrello M, Noirot C, Debarre F and Manel S. (under review). METAPOPGEN 2.0: a multi-locus genetic simulator to model populations of large size. Molecular Ecology Resources **Issue**, pp-pp.
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

# Gamete frequencies
fgam <- freq_alleles(N_deme1_age1_time3, init.par)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
