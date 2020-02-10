\name{sim.metapopgen.monoecious.multilocus}
\alias{sim.metapopgen.monoecious.multilocus}
\alias{sim.metapopgen.dioecious.multilocus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Main function to simulate metapopulation genetics\cr(multiple loci) 
}
\description{
  Simulates the metapopulation genetics.
}
\usage{
sim.metapopgen.monoecious.multilocus(init.par,
                                    sigma,
                                    phi_F, phi_M,
                                    delta, delta.ad,
                                    recr.dd="settlers", kappa0,
                                    T_max,
                                    save.res=FALSE, save.res.T=seq(1:T_max),
                                    verbose=FALSE)
                                    
sim.metapopgen.dioecious.multilocus(init.par,
                                    sigma_F, sigma_M,
                                    phi_F, phi_M,
                                    delta, delta.ad,
                                    recr.dd="settlers", kappa0,
                                    T_max,
                                    save.res=FALSE, save.res.T=seq(1:T_max),
                                    verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init.par}{
        A list of arguments necessary to initialize the simulations and obtained with the function initialize.monoecious.multilocus. See the example.
        }
  \item{sigma, sigma_F, sigma_M}{
    Survival probabilities. 4-D array of length c(m,n,z,T_max). \code{sigma[j,i,x,t]} is the survival probability for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. Alternatively, \code{sigma} can be a 3-D array of length c(m,n,z): this will be interpreted as survival probabilities constant in time. All the elements of \code{sigma} are bounded between [0,1].
  }
  \item{phi_F, phi_M}{
    Mean female fecundity. 4-D array of length c(m,n,z,T_max). \code{phi_F[j,i,x,t]} is the female fecundity for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. Alternatively, \code{phi_F, phi_M} can be 
    a 3-D array of length c(m,n,z): this will be interpreted as female or male fecundities constant in time. Female fecundity is the per capita number of female gametes (eggs for animals, embrio sacs for plants). Male fecundity is the per capita number of male gametes (sperms for animals, pollen grains for plants). All the elements of \code{phi_F} and \code{phi_M} must be positive integers. Variation in fecundity among individuals follows a Poisson distribution.
  }
  \item{delta}{
    Dispersal probabilities for propagules. Square matrix of size \code{n*n}. \code{delta[i,j]} is the dispersal probability from deme \code{j} to deme \code{i} (see also Details). Alternatively, \code{delta} can be a 3-D array of length c(n,n,T_max): this will be interpreted as dispersal varying in time.
  }
  \item{delta.ad}{
    Dispersal probabilities for adults. Square matrix of size \code{n*n}. \code{delta.add[i,j]} is the dispersal probability from deme \code{j} to deme \code{i} (see also Details). Alternatively, \code{delta.ad} can be a 3-D array of length c(n,n,z): this will be interpreted as dispersal varying with age; or \code{delta.ad} can be a 4-D array of length c(n,n,z,T_max), which will be interpreted as dispersal varying with age and in time. To set dispersal probabilities constant with age and varying in time, use the 4-D format.
  }
  \item{recr.dd}{
    Defines how settler survival probabilities depend on density. If "settlers" (default), settler survival probabilities depends on the number of settlers, through the function \code{\link[MetaPopGen:settler.survival]{settler.survival}}. If "adults", settler survival probabilities depends on the number of adults (all individuals already present in the deme): only up to \code{kappa0[i,t]}  individuals are allowed in a deme (see argument \code{kappa0}).
  }
  \item{kappa0}{
    2-D array of length c(n,T_max). Carrying capacity used to calculate the density-dependent settler survival probabilities. If the argument \code{recr.dd = "settlers"}, \code{kappa0[i,t]} is the settler carrying capacity for deme i at time t (maximum number of settlers allowed to recruit), and is used in the function \code{\link[MetaPopGen:settler.survival]{settler.survival}}. If \code{recr.dd = "adults"}, \code{kappa0[i,t]} is the maximum number of individuals allowed for deme i at time t. \code{kappa0} can also be a 1-D array of length n: this will be interpreted as a carrying capacity constant in time.
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
  If TRUE, the function will print much information, mainly used for debugging purposes. Default to FALSE.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
Element \code{i,j} of propagule and adult dispersal matrices always indicates dispersal from deme \code{j} to deme \code{i}. Note that this may be not equal to the dispersal probability from deme \code{i} to deme \code{j}. As a consequence, the sum of the elements of each column is bounded to one, while the sum of the elements of rows is unbounded.
}
\value{
  If \code{save.res} is FALSE: a named list with one element, \code{N} a 4-D array (ex: \code{N[k,i,x,t])} of length c(m,n,z,tmax) giving the number individuals of genotype \code{k}, age \code{x} in deme \code{i} at time \code{t}. \cr
  \cr
  If \code{save.res} is TRUE: the function does not return anything. Results are saved in files on disk. The function will create a folder named with the date and time of your simulation and create one RData file per year of simulation (but see argument \code{save.res.T} for defining the time steps for which the results are saved). Each file will contain the variable N, 3-D arrays (\code{N[k,i,x]}) of length c(m,n,z) giving the number of individuals of genotype \code{k}, age \code{x} in deme \code{i} at that time.
}
\references{
}
\author{
Christelle Noirot\cr
Marco Andrello
}
\note{
}
\seealso{
}
\examples{
\dontrun{
data("MultilocusMonoecious")
N <- sim.metapopgen.monoecious.multilocus(init.par=init.par,
                                          sigma=sigma,
                                          phi_F=phi_F, phi_M=phi_M,
                                          delta=delta, delta.ad=delta.ad,
                                          T_max=T_max)

}
}

