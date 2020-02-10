#calcul the allele frequency in a population N 

freq.all = function(N){
  m <- length(N)              # Number of genotypes
  l <- (sqrt(1+8*m)-1)/2      # Number of alleles
  p <- array(0,dim=l)
  ntot=sum(N)
  pos.gen <- genotype.index(l)
  for (i in 1 : l) {
    p[i] <- 0.5 * ( sum(N[pos.gen[i,]]) + N[pos.gen[i,i]] ) / ntot
  }
  return(p)
}