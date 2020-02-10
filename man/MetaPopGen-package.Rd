\name{MetaPopGen-package}
\alias{MetaPopGen-package}
\alias{MetaPopGen}
\docType{package}
\title{
Simulation of metapopulation genetics
}
\description{
Functions and code to simulate metapopulation genetics. Includes genetic drift, mutation, selection and migration. The metapopulation has age-structure. Migration is limited to the first age class (like seeds in plants and larvae in marine organisms). Functions are available to simulate one-locus genetics or for multi-locus simulations.
}
\details{
\tabular{ll}{
Package: \tab MetaPopGen\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2020-02-07\cr
License: \tab GPL-2\cr
}
This package simulates metapopulation genetics accounting for genetic drift, migration, selection and mutation. 
}
\author{
Marco Andrello \email{marco.andrello@gmail.com} \cr
Christelle Noirot
}
\references{
Andrello, M., and S. Manel. 2015. MetaPopGen: an R package to simulate population genetics in large size metapopulations. Molecular Ecology Resources 15:1153-1162.
}

\keyword{ package }
\seealso{
The best way to get started with MetaPopGen is to follow the tutorials:
https://github.com/MarcoAndrello/MetaPopGen-2.0/blob/master/Tutorial1.pdf
https://github.com/MarcoAndrello/MetaPopGen-2.0/blob/master/Tutorial2.pdf

The main functions for simulating one-locus metapopulation genetics are \code{\link[MetaPopGen:sim.metapopgen.monoecious]{sim.metapopgen.monoecious}} and \code{\link[MetaPopGen:sim.metapopgen.dioecious]{sim.metapopgen.dioecious}}.\cr
There are several example datasets, such as \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}}.\cr
For multilocus simulations, see functions \code{\link[MetaPopGen:sim.metapopgen.monoecious.multilocus]{sim.metapopgen.monoecious.multilocus}} and \code{\link[MetaPopGen:sim.metapopgen.dioecious.multilocus]{sim.metapopgen.dioecious.multilocus}}
}
\examples{

}
