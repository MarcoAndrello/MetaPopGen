\name{ld}
\alias{ld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linkage disequilibrium
}
\description{
  Linkage disequilibrium for two loci with two alleles.
}
\usage{
ld(N, init.par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Vector (1-D array) giving the number of individuals of each multilocus genotype. Typically this is obtained from the output of a multilocus simulation by subsetting a single deme / age class / time step. See the example.
}
  \item{init.par}{
A list of arguments necessary to initialize the simulations and obtained with the function initialize.monoecious.multilocus. It must be the same list that has been used to simulate \code{N}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Linkage disequilibrium (gametic-phase imbalance) calculated as f(A1B1) - f(A1)*f(B1), where f(A1B1) is the frequency of gamete A1B1 and f(A1) and f(B1) are the frequencies of alleles A1 and B1. Possible only for results of multilocus simulations using two loci with two alleles each.
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
data(MultilocusMonoecious)
N <- sim.metapopgen.monoecious.multilocus(init.par=init.par,
                                          sigma=sigma,
                                          phi_F=phi_F, phi_M=phi_M,
                                          delta.prop=delta.prop, delta.ad=delta.ad,
                                          T_max=T_max)
# Subsetting to get results from deme 1, age class 1 at time 3
N_deme1_age1_time3 <- N[,1,1,3]
ld(N_deme1_age1_time3, init.par)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
