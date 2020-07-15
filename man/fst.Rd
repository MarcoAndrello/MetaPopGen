\name{fst}
\alias{fst}
\alias{fst_multilocus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
F-statistics 
}
\description{
  F_ST for a single locus (\code{fst}) or for multiple loci (\code{fst_multilocus}).
}
\usage{
fst(N)
fst(N, init.par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
A 2-D array. It gives the number of individuals of each genotype (first dimension) in each group (second dimension). Typically this is obtained from the output of a single-locus (for \code{fst}) or multilocus (for \code{fst_multilocus}) simulation by subsetting a single deme / age class / time step. See the example.
}
  \item{init.par}{
A list of arguments necessary to initialize the simulations and obtained with the function initialize.monoecious.multilocus. It must be the same list that has been used to simulate \code{N}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The functions calculate the genetic differentiation index F_ST among the groups contained in \code{N} as ratio of heterozygosities: \code{(Ht - Hs) / Ht}. \code{Ht} is the expected heterozygosity over the population formed by pooling the groups together. \code{Hs} is the average expected heterozygosity over the groups. The average is weighted by the number of individuals in each group.\cr\cr
The multilocus function \code{fst_multilocus} returns F_ST for each locus and mean F_ST over all loci, calculated by taking the means of \code{Ht} and \code{Hs} over all loci.
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
# Generic single-locus version
# Three groups with two alleles
N1 <- c(18, 299, 277)
N2 <- c(284, 14, 145)
N3 <- c(115, 228, 254)
N <- array(c(N1,N2,N3),dim=c(3,3),dimnames=list(genotype=c("A1A1","A1A2","A2A2"),group=c(1:3)))
fst(N)

# Multi-locus version
data(MultilocusMonoecious)
N <- sim.metapopgen.monoecious.multilocus(init.par=init.par,
                                          sigma=sigma,
                                          phi_F=phi_F, phi_M=phi_M,
                                          delta.prop=delta.prop, delta.ad=delta.ad,
                                          T_max=T_max)
# Subsetting to get results from all demes, age class 1 at time 3
N_age1_time3 <- N[,,1,3]
fst_multilocus(N_age1_time3, init.par)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
