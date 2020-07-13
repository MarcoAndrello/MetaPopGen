# Observed and expected heterozygosities

# Observed heterozygosity
het_obs <- function(N, init.par, fG = NULL) {
  nLoc <- length(allele_vec)
  if(is.null(fG)) {
    fG <- freq_genotypes(N,init.par)
  }
  het_obs <- rep(NA,nLoc)
  names(het_obs) <- names(fG$counts)
  for(i.locus in 1 : nLoc) {
    names_genotypes <- names(fG$counts[[i.locus]])
    
    # Use package "genetics" to find heterozygous and homozygous genotypes
    if(any(sapply(names_genotypes,nchar)>4)) stop("All single-locus genotype names must be 4 character length")
    names_genotypes_genetics <- genetics::genotype(names_genotypes,sep=2)
    id_heterozygotes <- which(genetics::heterozygote(names_genotypes_genetics))
    
    # Observed heterozygosity
    het_obs[i.locus] <- sum(fG$counts[[i.locus]][id_heterozygotes]) / sum(fG$counts[[i.locus]])
  }
  het_obs
}

# Expected heterozygosity
het_exp <- function(N, init.par, fA = NULL) {
  nLoc <- length(allele_vec)
  if(is.null(fA)) {
    fA <- freq_alleles(N,init.par)
  }
  het_exp <- rep(NA,nLoc)
  names(het_exp) <- names(fA$counts)
  for (i.locus in 1 : nLoc) {
    het_exp[i.locus] <- 1 - sum((fA$frequencies[[i.locus]])^2)
  }
  het_exp
}

