\name{sim.metapopgen.monoecious.multilocus}
\alias{sim.metapopgen.monoecious.multilocus}
\alias{sim.metapopgen.dioecious.multilocus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate metapopulation genetics (multiple loci) 
}
\description{
  Simulates the metapopulation genetics.
}
\usage{
sim.metapopgen.monoecious.multilocus(init.par,
                                    sigma,
                                    phi_F, phi_M,
                                    fec.distr_F, fec.distr_M,
                                    migration="forward",
                                    delta.prop, delta.ad, migr,
                                    recr.dd="settlers", kappa0,
                                    T_max,
                                    save.res=FALSE, save.res.T=seq(1:T_max),
                                    output.var = "N",
                                    verbose=FALSE)
                                    
sim.metapopgen.dioecious.multilocus(init.par,
                                    sigma_F, sigma_M,
                                    phi_F, phi_M,
                                    delta.prop, delta.ad,
                                    recr.dd="settlers", kappa0,
                                    T_max,
                                    save.res=FALSE, save.res.T=seq(1:T_max),
                                    verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init.par}{
        A list of arguments necessary to initialize the simulations and obtained with the function \code{\link{initialize.multilocus}}. See the example.
        }
  \item{sigma, sigma_F, sigma_M}{
    Survival probabilities. 4-D array of length c(m,n,z,T_max). \code{sigma[j,i,x,t]} is the survival probability for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. Alternatively, \code{sigma} can be a 3-D array of length c(m,n,z): this will be interpreted as survival probabilities constant in time. All the elements of \code{sigma} are bounded between [0,1].
  }
  \item{phi_F, phi_M}{
    Mean female and male fecundity. 4-D array of length c(m,n,z,T_max). \code{phi_F[j,i,x,t]} is the female fecundity for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. Alternatively, \code{phi_F} and \code{phi_M} can be a 3-D array of length c(m,n,z): this will be interpreted as female or male fecundities constant in time. Female fecundity is the per capita number of female gametes (eggs for animals, embrio sacs for plants). Male fecundity is the per capita number of male gametes (sperms for animals, pollen grains for plants). All the elements of \code{phi_F} and \code{phi_M} must be positive integers.
  }
    \item{fec.distr_F, fec.distr_M}{
The distribution used to set fecundities: one of \code{"fixed"} or \code{"poisson"}. The \code{"fixed"} option sets the fecundity of each individual equal to the mean (\code{phi_F} and \code{phi_M}). Option \code{"poisson"} draws individual fecundities form a Poisson distribution with means \code{phi_F} and \code{phi_M}.
  }
      \item{migration}{
The type of migration used, \code{"forward"} (the default) or \code{"backward"}.\cr\cr
With forward migration, dispersal of propagules is modelled using the \code{delta.prop} argument (see below) while argument \code{migr} is disregarded.\cr\cr
With backward migration, migration is modelled using argument \code{migr} as in a finite island model (Latter, 1973). In this model, the gene pool of the jth deme contains a proportion \code{1-migr} of gametes derived at random from adults resident in the jth deme in the previous generation, and a proportion \code{migr/(n-1)} of gametes derived at random from adults resident in each of the remaining \code{n-1} demes. Variables \code{L} and \code{S} are not calculated with backward migration (see argument \code{output.var}), because the number of individuals recruiting the deme is calculated directly from the gametes of the previous generations.\cr
The function will return a warning if \code{(1-migr)*kappa0[j,t+1]} or \code{migr/(n-1)*kappa0[j,t+1]} are not integers, where \code{kappa0[j,t+1]} is the number of individuals that recruits in the generation that is being formed. Backward migration is currently limited to only one age class and monoecious life cycles.\cr\cr
\code{migration} does not affect adult dispersal.
  }
  
  \item{delta.prop}{
    Dispersal probabilities for propagules. Square matrix of size \code{n*n}. \code{delta.prop[i,j]} is the dispersal probability from deme \code{j} to deme \code{i} (see also Details). Alternatively, \code{delta.prop} can be a 3-D array of length c(n,n,T_max): this will be interpreted as dispersal varying in time.
  }
  \item{delta.ad}{
    Dispersal probabilities for adults. Square matrix of size \code{n*n}. \code{delta.ad[i,j]} is the dispersal probability from deme \code{j} to deme \code{i} (see also Details). Alternatively, \code{delta.ad} can be a 3-D array of length c(n,n,z): this will be interpreted as dispersal varying with age; or \code{delta.ad} can be a 4-D array of length c(n,n,z,T_max), which will be interpreted as dispersal varying with age and in time. To set dispersal probabilities constant with age and varying in time, use the 4-D format.
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
    \item{output.var}{
A vector containing the names of the variables output by the simulation: \code{"N"}, \code{"Nprime"}, \code{"Nprimeprime"}, \code{"L"} and \code{"S"}. For example, to output genotype numbers in adults and propagules, set \code{output.var = c("N","L")}. Default to \code{"N"}. 
  }

  \item{verbose}{
  If TRUE, the function will print additional information, mainly used for debugging purposes. Default to FALSE.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
Element \code{i,j} of propagule and adult dispersal matrices always indicates dispersal from deme \code{j} to deme \code{i}. Note that this may be not equal to the dispersal probability from deme \code{i} to deme \code{j}. As a consequence, the sum of the elements of each column is bounded to one, while the sum of the elements of rows is unbounded.
}

\value{
\itemize{
  \item {If \code{save.res = FALSE}:\cr\cr
  a named list where each element is one of the tracked variables (see argument \code{output.var}). Each of these elements is an array giving the number of individuals of each genotype and deme. Specifically,
  \describe{
  \item{\code{N}}{
a 4-D array (ex: \code{N[k,i,x,t])} of length c(m,n,z,T_max) giving the number of individuals of genotype \code{k}, age \code{x} in deme \code{i} at time \code{t}. This is the number of individuals at the beginning of the life cycle, just after recruitment and before adult survival.
}
\item{\code{Nprime}}{
a 4-D array (ex: \code{Nprime[k,i,x,t])} of length c(m,n,z,T_max) giving the number of individuals of genotype \code{k}, age \code{x} in deme \code{i} at time \code{t} after adult survival and before adult dispersal.
}
\item{\code{Nprimeprime}}{
a 4-D array (ex: \code{Nprimeprime[k,i,x,t])} of length c(m,n,z,T_max) giving the number of individuals of genotype \code{k}, age \code{x} in deme \code{i} at time \code{t} after adult dispersal.
}
\item{\code{L}}{
a 3-D array (ex: \code{L[k,i,t])} of length c(m,n,T_max) giving the number of propagules of genotype \code{k}, in deme \code{i} at time \code{t} before propagule dispersal.
}
\item{\code{S}}{
a 3-D array (ex: \code{S[k,i,t])} of length c(m,n,T_max) giving the number of propagules of genotype \code{k}, in deme \code{i} at time \code{t} after propagule dispersal.
}
}
}

  \item {If \code{save.res = TRUE}:\cr\cr
  the function does not return anything. Results are saved in files on disk. The function will create a folder named with the date and time of your simulation and create one RData file per year of simulation (but see argument \code{save.res.T} for defining the time steps for which the results are saved) and per chosen output variables (\code{N}, \code{Nprime}, \code{Nprimeprime}, \code{L} and \code{S}). The files will be named with the name of the output variable and the time step of simulation, e.g. N1.RData, S1.RData, etc.\cr\cr
  The variables contained in each file are as described above for the case (\code{save.res = FALSE}), but they do not have the time dimension. Specifically, \code{N}, \code{Nprime} and \code{Nprimeprime} have length c(m,n,z), and \code{L} and \code{S} have length c(m,n).
  }
  }
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
\code{\link{initialize.multilocus}}
\code{\link{create.dispersal.coord}}
\code{\link{create.dispersal.IM}}
}
\examples{
\dontrun{
data("MultilocusMonoecious")
res <- sim.metapopgen.monoecious.multilocus(init.par=init.par,
                                            sigma=sigma,
                                            phi_F=phi_F, phi_M=phi_M,
                                            fec.distr_F="fixed", fec.distr_M="fixed",
                                            delta.prop=delta.prop, delta.ad=delta.ad,
                                            T_max=T_max)

}
}

