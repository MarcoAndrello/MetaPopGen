\name{create.multilocus.rate}
\alias{create.multilocus.rate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create multilocus vital rates
}
\description{
  This function returns vital rates (survival probabilities or fecundities) for every genotype, deme, age-class and time as a function of the genetic architecture (total number of loci, number of selected loci, number of alleles per locus) and an environmental variable (selective variable) changing over demes. See Andrello et al. (under review) for a description of the function.
}
\usage{
create.multilocus.rate(allele_vec, init.par, n.sel, x, xi,
                       rate.max, omega = 1, T_max)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allele_vec}{
    Vector giving the number of alleles at each locus. Currently accepts only biallelic loci under selection, while the neutral loci can have more than two alleles.
  }
    \item{init.par}{
        A list of arguments necessary to initialize the simulations and obtained with the function initialize.monoecious.multilocus. See the example.
        }
  \item{n.sel}{
    Number of loci under selection. The first \code{n.sel} loci in \code{allele_vec} are considered under selection. Must be lower than the length of \code{allele_vec}. 
  }
  \item{x}{
  Vector of values of the environmental variable in each deme. The environmental variable acts as a selection pressure on survival through the \code{n.sel} loci under selection. The length of \code{x} must equal the number of demes \code{n}.
}
  \item{xi}{
  Vector of optimal values of the environmental variable for each genotype "class". The optimal values is the one for which survival of that genotype "class" is the highest. The elements of the \code{xi} vector are usually monotonically increasing or decreasing. The length of the \code{xi} vector is equal to the number of genotype "classes", which is \code{2*n.sel+1}. See details for the definition of "classes".
}
  \item{rate.max}{
  Maximal vital rate. This is the vital rate of a genotype "class" when the environmental value in a deme equals the  optimal environmental value of that "class".   
}
    \item{omega}{
Parameter that can be interpreted as an inverse of selection strength of the environment on individuals. Lower values result in stronger selection. it must be strictly positive. }
  \item{T_max}{
    Number of years of simulations.
}
}

\details{
This function is used to set survival probabilities and fecundities as functions of environmental conditions and the genotype of biallelic loci (e.g. SNPs). This can be understood as if the genotype j determines the optimal environmental condition \eqn{\xi_j} for the organism. It is assumed that, at each locus, one allele reduces the optimal condition (the - allele) and the other increases it (the + allele). \eqn{\xi_j} is thus a function of the number of + alleles, and the vital rate (survival or fecundity) is an exponential function of the difference between the environmental condition of deme i, \eqn{x}, and the optimal environmental condition for the genotype \eqn{\xi_j}:

\deqn{rate.max*exp[ (-(x_i - \xi_j )^2) / (2\omega^2) ]}

If the number of biallelic loci under selection is n.sel, then there are \code{2*n.sel+1} combinations of - and + alleles ("genotype "classes"). This means thatthe number of + alleles in a multilocus genotype can range from 0 to is \code{2*n.sel}. It is thus necessary to define \code{2*n.sel+1} values for the optimal environmental conditions, each corresponding to a combination of - and + alleles (a genotype "class"). See the example.


}
\value{
Vital rates (survival probabilities or fecundities) to be used as input in \code{\link[MetaPopGen:sim.metapopgen.monoecious.multilocus]{sim.metapopgen.monoecious.multilocus}} and \code{\link[MetaPopGen:sim.metapopgen.dioecious.multilocus]{sim.metapopgen.dioecious.multilocus}}. 4-D array of length c(m,n,z,T_max). For example, \code{sigma[j,i,x,t]} is the survival probability for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. 
}



\examples{
# Defining four biallelic loci, of which the first three are under selection
allele_vec <- c(2,2,2,2)
n.sel <- 3
# Initializing the simulations: three demes. ?initialize.multilocus for details
init.par <- initialize.multilocus(allele_vec=allele_vec,
                                  r=0.5,
                                  mu=c(0.01,0.01),
                                  n=3,
                                  z=1,
                                  kappa0=100,
                                  sexuality="monoecious")
# Defining the environmental conditions in the three demes (e.g. values of water salinity per deme)         
x <- c(35,39,38)
# Defining the optimal conditions for each "genotype class"
# Because there are 3 loci under selection, the number of combinations of - and + alleles is 2*3+1 = 7 classes
# # We need to define 7 values:
xi <- c(36, 36.5, 37, 37.5, 38, 38.5, 39) # (e.g. optimal water salinity per phenotype)
# Creating survival dependent on water salinity and genotype
sigma <- create.multilocus.rate(allele_vec,
                                    init.par,
                                    n.sel=n.sel,
                                    x,
                                    xi,
                                    rate.max=0.8,
                                    omega=1,
                                    T_max=5)
}
\author{
Marco Andrello
}
\references{
Andrello M, Noirot C, Debarre F and Manel S. (under review). METAPOPGEN 2.0: a multi-locus genetic simulator to model populations of large size. Molecular Ecology Resources **Issue**, pp-pp.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
