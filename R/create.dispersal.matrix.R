# Create dispersal matrix from deme coordinates
create.dispersal.matrix = function(x,a,longlat = FALSE) {
  require(sp)
  dist.mat <- spDists(x,longlat=longlat)
  y   <- 0.5*a*exp(-a*abs(dist.mat))
  y
}