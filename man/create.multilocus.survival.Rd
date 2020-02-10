\name{create.multilocus.survival}
\alias{create.multilocus.survival}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create multilocus survival
}
\description{
  This function returns survival probabilities for every genotype, deme, age-class and time as a function of the genetic architecture (total number of loci, number of selected loci, number of alleles per locus) and an environmental variable (selective variable) changing over demes. See Andrello et al. (submitted) for a description of the function.
}
\usage{
create.multilocus.survival(allele_vec, n.sel, x, xi, surv.max, z, T_max, omega=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allele_vec}{
    Vector giving the number of alleles at each locus. Currently accepts only biallelic loci under selection, while the neutral loci can have more than two alleles.
  }
  \item{n.sel}{
    Number of loci under selection. The first \code{n.sel} loci in \code{allele_vec} are considered under selection. Must be lower than the length of \code{allele_vec}. 
  }
  \item{x}{
  Vector of values of the environmental variable in each deme. The environmental variable acts as a selection pressure on survival through the \code{n.sel} loci under selection. The length of \code{x} must equal the number of demes \code{n}
}
  \item{xi}{
  Vector of optimal values of the environmental variable for each phenotype. The optimal values is the one for which survival of that phenotype is the highest. The elements of the \code{xi} vector are usually monotonically increasing or decreasing. The length of the \code{xi} vector is equal to the number of survival phenotypes, which is \code{K = 2*n.sel+1}
}
  \item{surv.max}{
  Maximal survival probability, comprised between 0 and 1. This is the survival probability of a phenotype when the envronmental value in a deme equals its optimal value.   
}
  \item{z}{
    Number of age-classes.
}
  \item{T_max}{
    Number of years of simulations.
}
  
    \item{omega}{
Parameter that can be interpreted as an inverse of selection strength of the environment on individuals. Lower values result in stronger selection. it must be strictly positive.  
  }
}
\value{
Survival probabilities to be used as input in \code{\link[MetaPopGen:sim.metapopgen.monoecious.multilocus]{sim.metapopgen.monoecious.multilocus}} and \code{\link[MetaPopGen:sim.metapopgen.dioecious.multilocus]{sim.metapopgen.dioecious.multilocus}}. 4-D array of length c(m,n,z,T_max). \code{sigma[j,i,x,t]} is the survival probability for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. 
}
\examples{
allele_vec <- c(2,2,2,2)
init.par <- initialize.multilocus(allele_vec=allele_vec,
                                  r=0.5,
                                  mu=c(0.01,0.01),
                                  n=3,
                                  z=1,
                                  kappa0=100,
                                  sexuality="monoecious")
x <- c(35,39,38) # Environmental values (e.g. water salinity per deme)
xi <- seq(36,39,0.5) # Phenotype optima (e.g. optimal water salinity per phenotype)
sigma <- create.multilocus.survival(allele_vec,
                                    init.par,
                                    n.sel=3,
                                    x,
                                    xi,
                                    surv.max=0.8,
                                    omega=1,
                                    T_max=5)
}
\author{
Marco Andrello
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
