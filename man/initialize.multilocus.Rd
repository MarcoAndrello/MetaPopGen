\name{initialize.multilocus}
\alias{initialize.multilocus}
\title{
  Initialize multilocus simulations
}
\description{
Function used to initialize multilocus simulations.
}
\usage{
initialize.multilocus(allele_vec,
                      r,
                      mu=rep(0,length(allele_vec)),
                      n,
                      z,
                      kappa0,
                      init.state="fixed",
                      sexuality=NULL)
}

\arguments{
  \item{allele_vec}{
Vector giving the number of alleles at each locus.
}
  \item{r}{
Recombination rate. A single value between 0 and 0.5. For more than two loci, recombination rate must be set to 0.5.
}
  \item{mu}{
Vector giving the mutation rate for each locus, i.e. the probability that an allele mutate to a different allele. 
}
  \item{n}{
Number of demes. 
}
  \item{z}{
Number of age classes. 
}
  \item{kappa0}{
Initial number of individuals per deme. This is used to assign individuals to genotype and age-classes in each deme according to the criterion defined in \code{"init.state"}. Because the assignment involves rounding, the realized number of individuals might be different than the one set by the user.
}
  \item{init.state}{
String defining the genetic composition of the initial population. If \code{"fixed"}, all genotypes and age-classes will have the same number of individuals. If \code{"random"}, numbers of individuals will be assigned randomly.
}
  \item{sexuality}{
One of \code{"monoecious"} or \code{"dioecious"}, defining the sexuality of the modelled species.
}
}

\details{
}
\value{
A list with the following items:
\item{r}{
The recombination rate input by the user
}
\item{mu}{
The vector of mutation rates input by the user
}
\item{m}{
The number of multilocus genotypes
}
\item{n}{
The number of demes input by the user
} 
\item{z}{
The number of age-classes input by the user
}
\item{kappa0}{
 The initial number of individuals per deme given by the user
 }
\item{index_matr}{
The index of the multilocus genotypes given by the function \code{\link{genotype.index.multilocus}}
 }
\item{meiosis_matrix}{
A matrix giving the probabilities of producing multilocus gametes (rows) for each multilocus genotype (columns)
}
\item{mat_geno_to_index_mapping}{
A square matrix giving the genotype index for each combination of female (rows) and male (columns) multilocus gametes
}
\item{N1}{
The initial composition of the population as a 3D array giving the number of individuals of each genotype (first dimension) in each deme (second dimension) in each age-class (third dimension)
}

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
\dontrun{
allele_vec <- c(2,2)
init.par <- initialize.multilocus(allele_vec=allele_vec,
                                  r=0.5,
                                  mu=c(0.01,0.01),
                                  n=3,
                                  z=1,
                                  kappa0=100,
                                  sexuality="monoecious")

}
}

