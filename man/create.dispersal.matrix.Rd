\name{create.dispersal.matrix}
\alias{create.dispersal.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create dispersal matrix
}
\description{
Function to create a simple dispersal matrix based on spatial distance between demes.
The function returns a matrix of dispersal probabilities between demes. It first calculates the spatial distance between demes, and then the dispersal probability using a Laplacian dispersal function:
\code{p(i,j) = 0.5 * a * exp(-a*dist(i,j))}, where \code{dist(i,j)} is the spatial distance between deme \code{i} and deme \code{j} and \code{1/a} is the mean dispersal distance.
}
\usage{
create.dispersal.matrix(x,a,longlat = FALSE)
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
create.dispersal.matrix(x,0.01,longlat=T)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
