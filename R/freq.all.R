# Calculates genotype and allele frequencies 

# freq.all
# change_allele_order
# freq_genotypes

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

change_allele_order <- function(corr_locus){
  new.corr_locus <- corr_locus
  for (i.genotype in 1 : length(corr_locus)) {
    all.female <- as.numeric(substr(corr_locus[i.genotype],2,2))
    all.male <- as.numeric(substr(corr_locus[i.genotype],4,4))
    if (all.female > all.male) {
      all.female.ext <- substr(corr_locus[i.genotype],1,2)
      all.male.ext <- substr(corr_locus[i.genotype],3,4)
      new.corr_locus[i.genotype] <- paste0(all.male.ext, all.female.ext)
    }
  }
  new.corr_locus
}


freq_genotypes <- function(N, init.par) {
  # Check if N has names...
  if (is.null(names(N))) stop("N must be a named vector. Multilocus simulations with MetaPopGen return named objects. You must have done something wrong!")
  
  if (init.par$method == "locus") {
    # Check no more than 26 loci...
    if (length(init.par$allele_vec) > 26) stop("Not yet implemented for more than 26 loci")
    nLoc <- length(allele_vec)
    # Set genotype_counts vectors, one per locus
    genotype_counts <- list()
    for (i.locus in 1 : nLoc){
      genotype_counts[[i.locus]] <- array(0,dim=dim(init.par$index_matr)[i.locus])
      names(genotype_counts[[i.locus]]) <- dimnames(init.par$index_matr)[[i.locus]]
    }
    names(genotype_counts) <- paste0("Locus",LETTERS[1:nLoc])
    
    # loop on multilocus genotypes
    for (i.genotype in 1 : init.par$m){
      index_loci <- which(init.par$index_matr == i.genotype, arr.ind = T) # vector of indices for each locus
      for (i.locus in 1 : nLoc) {
        # Check if index_loci is always a 1 by nLoc matrix
        # print(dim(index_loci))
        id_genotype <- index_loci[1,i.locus]
        genotype_counts[[i.locus]][id_genotype] <- genotype_counts[[i.locus]][id_genotype] + N[i.genotype]
      }
    }
    
  } else {
    nLoc <- 2 # Valid only for two loci
    # Names of single-locus genotypes
    names_genotypes <- def_genotype.name.locus(init.par$allele_vec)
    
    # Names of multilocus genotypes
    names_multilocus_genotypes <- names(N)
    
    # Correspondence to single-locus genotypes
    if(any(allele_vec>9)) stop("More than 10 alleles per locus is not implemented yet")
    
    multi_to_single_geno <- list()
    # First locus
    corr_locus <- paste0(substr(names_multilocus_genotypes,1,2), substr(names_multilocus_genotypes,6,7))
    new.corr_locus <- change_allele_order(corr_locus)
    multi_to_single_geno[[1]] <- match(new.corr_locus, names_genotypes[[1]])
    # Second locus
    corr_locus <- paste0(substr(names_multilocus_genotypes,3,4), substr(names_multilocus_genotypes,8,9))
    new.corr_locus <- change_allele_order(corr_locus)
    multi_to_single_geno[[2]] <- match(new.corr_locus, names_genotypes[[2]])
    names(multi_to_single_geno) <- c("LocusA", "LocusB")
    
    # Set genotype_counts vectors, one per locus
    genotype_counts <- list()
    for (i.locus in 1 : nLoc){
      genotype_counts[[i.locus]] <- array(0,dim=length(names_genotypes[[i.locus]]))
      names(genotype_counts[[i.locus]]) <- names_genotypes[[i.locus]]
    }
    names(genotype_counts) <- paste0("Locus",LETTERS[1:nLoc])
    
    # Loop on multilocus genotypes
    for (i.genotype in 1 : init.par$m){
      for (i.locus in 1 : nLoc) {
        id_single_locus <- multi_to_single_geno[[i.locus]][i.genotype]
        genotype_counts[[i.locus]][id_single_locus] <- genotype_counts[[i.locus]][id_single_locus] + N[i.genotype]
      }
    }
  }
  
  # Genotype frequencies
  genotype_frequencies <- list()
  for (i.locus in 1 : nLoc) {
    genotype_frequencies[[i.locus]] <- genotype_counts[[i.locus]] / sum(genotype_counts[[i.locus]])
  }
  names(genotype_frequencies) <- names(genotype_counts)
  
  # Output
  out <- list(counts = genotype_counts,
              frequencies = genotype_frequencies)
  out
  
}

freq_alleles <- function(N,init.par,fG=NULL) {
  # Get genotype frequencies
  if(is.null(fG)) {
    fG <- freq_genotypes(N,init.par)
  }
  
  # Allele counts
  allele_counts <- list()
  nLoc <- length(init.par$allele_vec)
  for (i.locus in 1 : nLoc) {
    allele_counts[[i.locus]] <- rep(0, allele_vec[i.locus])
    names_genotypes <- names(fG$counts[[i.locus]])
    
    # Use package "genetics" to find allele counts per genotype
    if(any(sapply(names_genotypes,nchar)>4)) stop("All single-locus genotype names must be 4 character length")
    names_genotypes_genetics <- genetics::genotype(names_genotypes,sep=2)
    allele_counts_per_genotype <- genetics::allele.count(names_genotypes_genetics)
    
    # Order alleles from 1 to n
    allele_counts_per_genotype <- allele_counts_per_genotype[,order(colnames(allele_counts_per_genotype))] 
    
    # Add genotype name (not needed, just for clarity in debugging)
    # rownames(allele_counts_per_genotype) <- names_genotypes
    
    num_alleles_per_locus <- allele_vec[i.locus]
    for (i.allele in 1 : num_alleles_per_locus) {
      allele_counts[[i.locus]][i.allele] <- sum(allele_counts_per_genotype[,i.allele] * fG$counts[[i.locus]])
    }
    names(allele_counts[[i.locus]]) <- colnames(allele_counts_per_genotype)
  }
  names(allele_counts) <- names(fG$counts)
  
  # Allele frequencies
  allele_frequencies <- list()
  for (i.locus in 1 : nLoc) {
    allele_frequencies[[i.locus]] <- allele_counts[[i.locus]] / sum(allele_counts[[i.locus]])
  }
  names(allele_frequencies) <- names(allele_counts)
  
  # Output
  out <- list(counts = allele_counts,
              frequencies = allele_frequencies)
  out
}

# Norig <- N
# N <- Norig[,1,1,3]
##


