# Create dispersal matrix from deme coordinates
create.dispersal.coord = function(x,a,longlat = FALSE) {
  require(sp)
  dist.mat <- spDists(x,longlat=longlat)
  y   <- 0.5*a*exp(-a*abs(dist.mat))
  y
}
create.dispersal.IM = function(n,m) {
  y <- matrix(m/(n-1),nrow=n,ncol=n)
  diag(y) <- (1-m)
  y
}
