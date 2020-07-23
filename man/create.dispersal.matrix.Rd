\name{create.dispersal.coord}
\alias{create.dispersal.coord}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create dispersal matrix
}
\description{
Functions to create dispersal matrices

}
\usage{
create.dispersal.coord(x,a,longlat = FALSE)
create.dispersal.IM(n,m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 A matrix of 2D (n-D) points with row denoting points, first column x/longitude, second column y/latitude (or a Spatial object that has a coordinates method).
}
  \item{a}{
Parameter of the Laplacian dispersal function. \code{1/a} is the mean dispersal distance.
}
\item{longlat}{
It defines how coordinates in \code{x} are used to calculate spatial distances. If \code{longlat = FALSE}, Euclidean distance; if TRUE, Great Circle distance.
}
\item{n}{
Number of demes
}
\item{m}{
Migration probability
}
}
\details{
\code{create.dispersal.coord} creates a dispersal matrix based on spatial distance between demes, using a Laplacian dispersal function:\cr \code{p(i,j) = 0.5 * a * exp(-a*dist(i,j))},\cr where \code{dist(i,j)} is the spatial distance between deme \code{i} and deme \code{j} and \code{1/a} is the mean dispersal distance.\cr\cr
\code{create.dispersal.IM} creates a dispersal matrix based on an island model. The probability to disperse to any \code{n-1} non-natal patch is \code{m/(n-1)} while the probability to stay in the natal patch is \code{1-m}.

}
\value{
A square symmetric matrix whose [i,j] element gives the probability of dispersal from deme j to deme i.
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
lon <- c(12.3, 13.3, 14.2, 16.9, 16.1, 17.5, 18.1, 19.7, 20.3)
lat <- c(45.4, 46.5, 48.9, 47.1, 45.8, 46.1, 47.0, 46.8, 47.2)
x <- cbind(lon,lat)
create.dispersal.coord(x,0.01,longlat=T)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
