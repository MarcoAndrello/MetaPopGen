# Generic FST function
fst <- function(N) {
  
  # Define basic variables
  m <- dim(N)[1]            # Number of genotypes
  l <- (sqrt(1+8*m)-1)/2    # Number of alleles
  n <- dim(N)[2]            # Number of groups
  
  # Genotype arrays and number of individuals for each deme
  N_i <- apply(N,c(1,2),sum)
  n_i <- colSums(N)
  
  # Genotype arrays and number of individuals for total population
  N_T <- rowSums(N)
  n_T  <- sum(N)
    
  # Calculate allele frequencies in each group and in the total population
  p_i <- array(NA,dim=c(n,l),dimnames=list(group=c(1:n),allele=c(1:l)))
  for (i in 1 : n) {
    p_i[i,] <- freq.all(N[,i])    
  }
  p_T	<- freq.all(N_T)
  
  # Calculate expected heterozygosities in each group and in the total population
  H_S_i <- array(NA,dim=n,dimnames=list(group=c(1:n)))
  for (i in 1 : n) {
    H_S_i[i] <- 1 - sum(p_i[i,]^2)
  }
  H_t	<- 1 - sum(p_T^2)
  
  # Calculate fst as ratio of heterozygosities
  fst <- ( H_t - sum(n_i * H_S_i,na.rm=T)/n_T)  / H_t
  return(fst)
  
}