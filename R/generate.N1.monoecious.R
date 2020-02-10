generate.N1.monoecious<- function(l,n,z,n.deme=100) {
  
  # l = number of alleles
  # n = number of demes
  # z = number of age.classes
  m = l * (l + 1) / 2 # Number of genotypes
  
  N1 <- array(100,dim=c(m,n,z))
#   
#   for (deme in 1 : n){
#     for (age in 1 : z) {
#       for (genotype in 1 : m) {
#         N1[genotype,deme,age] <- round(runif(1)*(n.deme*2/m))      
#       }
#     }
#   }
  
  return(N1)
  
}


