# Generic FST function for single locus and multi-locus simulations

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

fst_multilocus <- function(N, init.par) {
  if(length(dim(N))!=2) stop("N must have two dimensions: genotype and group")
  nLoc <- length(init.par$allele_vec) # Number of loci
  num_groups <- dim(N)[2]             # Number of groups
  locus_names <- paste0("Locus",LETTERS[1:nLoc])
  
  # Number of individuals for each group
  n_group <- colSums(N) # number of individuals for each group
  
  # Genotype arrays for the total population
  N_T <- rowSums(N)
  
  # H_S per locus per group
  H_S_group <- array(NA,c(nLoc,num_groups))
  dimnames(H_S_group) <- list(Locus = locus_names,
                              Deme = colnames(N))
  for(i.group in 1 : num_groups) {
    H_S_group[,i.group] <- het_exp(N[,i.group], init.par)
  }
  
  # H_S per locus as weighted means of the H_S per group
  H_S <- rep(0,nLoc)
  for(i.locus in 1 : nLoc) {
    H_S[i.locus] <- weighted.mean(H_S_group[i.locus,],n_group)
  }
  
  # H_T per locus
  H_T <- het_exp(N_T, init.par)
  
  # Fst per locus
  fst <- rep(NA,nLoc)
  for(i.locus in 1 : nLoc) {
    fst[i.locus] <- (H_T[i.locus] - H_S[i.locus]) / H_T[i.locus]
  }
  
  # Mean over loci
  fst_mean <- (mean(H_T) - mean(H_S)) / mean(H_T)
  
  # Out
  out <- c(fst,fst_mean)
  names(out) <- c(locus_names,"Mean")
  out
}

