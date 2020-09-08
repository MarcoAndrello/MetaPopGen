\name{sim.metapopgen.monoecious}
\alias{sim.metapopgen.monoecious}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate one-locus metapopulation genetics\cr(monoecious life cycles)
}
\description{
Simulates the metapopulation genetics. The function will calculate automatically from its arguments the number of alleles \code{l}, the number of genotypes \code{m}, the number of demes \code{n} and the number of age-classes \code{z}. It is therefore important that all the arguments have the correct dimensions and size. See the example \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}} for a detailed description of each argument, its biological meaning, required dimensions and dimension sizes.
}
\usage{
res <- sim.metapopgen.monoecious(input.type,
                                  demographic.data,
                                  N1,
                                  sigma,
                                  phi_F, phi_M,
                                  fec.distr_F = "poisson", fec.distr_M = "poisson",
                                  mu,
                                  delta,
                                  recr.dd="settlers", kappa0,
                                  T_max,
                                  save.res = F, save.res.T = seq(1,T_max),
                                  verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.type}{
Character vector defining whether the input data are giving as arrays (\code{"array"}) or as a data frame (\code{"data.frame"})
}
  \item{demographic.data}{
If the input is given as a data.frame (\code{input.type = "data.frame"}), a data frame giving 1) survival probabilities and fecundities for each genotype, deme, age-class and time and 2) number of individuals for each genotype, deme and age-class at the beginning of the simulation
}
  \item{N1}{
Number of individuals at the beginning of the simulation. 3-D array giving the number of individuals of each genotype in each deme and each age-class at the beginning of the simulation. The first dimension is the genotype, the second dimension is the deme and the third dimension is the age-class. To be used when \code{input.type = "array"}.
}
  \item{sigma}{
Survival probabilities. 4-D array of length c(m,n,z,T_max). \code{sigma[j,i,x,t]} is the survival probability for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. Alternatively, \code{sigma} can be 
a 3-D array of length c(m,n,z): this will be interpreted as survival probabilities constant in time. All the elements of \code{sigma} are bounded between [0,1]. To be used when \code{input.type = "array"}.
}
  \item{phi_F}{
Mean female fecundity. 4-D array of length c(m,n,z,T_max). \code{phi_F[j,i,x,t]} is the mean female fecundity for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. Alternatively, \code{phi_F} can be 
a 3-D array of length c(m,n,z): this will be interpreted as female fecundities constant in time. The female fecundity is the per capita number of female gametes (eggs for animals, embrio sacs for plants). All the elements of \code{phi_F} must be natural numbers. To be used when \code{input.type = "array"}.
}
  \item{phi_M}{
Mean male fecundity. 4-D array of length c(m,n,z,T_max). \code{phi_M[j,i,x,t]} is the mean male fecundity for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. Alternatively, \code{phi_M} can be 
a 3-D array of length c(m,n,z): this will be interpreted as male fecundities constant in time. The male fecundity is the per capita number of male gametes (sperms for animals, pollen grains for plants). All the elements of \code{phi_F} must be natural numbers. To be used when \code{input.type = "array"}.
}
    \item{fec.distr_F, fec.distr_M}{
The distribution used to set fecundities: one of \code{"fixed"} or \code{"poisson"} (default). The \code{"fixed"} option sets the fecundity of each individual equal to the mean (\code{phi_F} and \code{phi_M}). Option \code{"poisson"} draws individual fecundities form a Poisson distribution with means \code{phi_F} and \code{phi_M}.
  }
  
  \item{mu}{
Mutation rates. Square matrix of size \code{l*l}. \code{mu[i,j]} is the probability that allele \code{j} mutates into allele \code{i}. The elements of the diagonal of \code{mu} are the probabilities that the allele does not mutate, i.e. it must be \code{colSums(mu)} = 1
}
  \item{delta}{
Dispersal probabilities. 3-D array of length c(n,n,T_max). \code{delta[i,j,t]} is the dispersal probability from deme \code{j} to deme \code{i} at time \code{t}. Alternatively, \code{delta} can be a 2-D array of length c(n,n): this will be interpreted as dispersal constant in time.
}
  \item{recr.dd}{
Defines how population regulation is performed. Currently inactive.
}
  \item{kappa0}{
  2-D array of length c(n,T_max). Deme carrying capacity used to perform population regulation. \code{kappa0} can also be a 1-D array of length n: this will be interpreted as a carrying capacity constant in time.
}
  \item{T_max}{
  Number of years of simulations.
}
  \item{save.res}{
Boolean. If FALSE (default), the result of the simulation is stored in memory and returned to the user. If TRUE, the function does not return anything, but save the result in different files on disk. See "value" below.
}
  \item{save.res.T}{
Vector. If save.res=TRUE, this argument is used to define the time steps for which the state of the metapopulation will be saved on disk. It is used to save disk space in the case of long simulations with many demes or genotypes. It deafults to \code{seq(1,T_max)}, that is, the results are saved every time step.
}
  \item{verbose}{
Boolean. If TRUE, the function will print much information, mainly used for debugging purposes. Default to FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If \code{save.res} is FALSE: a 4-D array (ex: \code{N[k,i,x,t]}) of length c(m,n,z,tmax) giving the number of individuals of genotype \code{k}, age \code{x} in deme \code{i} at time \code{t}. \cr
If \code{save.res} is TRUE: the function does not return anything. Results are saved in files on disk. The function will create a folder named with the date and time of your simulation and create one RData file per year of simulation (but see argument \code{save.res.T} for defining the time steps for which the results are saved). Each file will contain the variable N, a 3-D array (\code{N[k,i,x]}) of length c(m,n,z) giving the number of individuals of genotype \code{k}, age \code{x} in deme \code{i} at that time.
}
\references{
}
\author{
}
\note{
}
\seealso{
}
\examples{
data(FourPopMonoecious)
N <- sim.metapopgen.monoecious(input.type="array",N1=N1, sigma=sigma, phi_F=phi_F, phi_M=phi_M,
  mu=mu, delta=delta, kappa0=kappa0, T_max=T_max)

# This will save results on disk:
sim.metapopgen.monoecious(input.type="array",N1=N1, sigma=sigma, phi_F=phi_F, phi_M=phi_M,
  mu=mu, delta=delta, kappa0=kappa0, T_max=T_max,save.res=T)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
