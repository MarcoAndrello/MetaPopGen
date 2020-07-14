# Gamete frequencies
freq_gametes <- function(N, init.par) {
  
  if (!init.par$method %in% c("gamete","locus")) stop("init.par$method must be either \"gamete\" or \"locus\"")
  
  # Set gamete_counts vector and gamete_names  
  num_gametes <- prod(init.par$allele_vec)
  gamete_names <- rownames(init.par$meiosis_matrix)
  gamete_counts <- rep(0,num_gametes)
  names(gamete_counts) <- gamete_names
  
  if (init.par$method == "locus") {
    # This is for locus based
    # Setting mutation to 0 and recombination to 0.5 will give a "recombination" matrix giving the parental gametes
    # As the method is "locus", phasing is ignored
    RECOMB <- create.meiosis.matrix(init.par$index_matr, init.par$allele_vec, r=0.5, mu = rep(0,length(init.par$allele_vec)), method="locus")
    # Multiply by 2 because each genotype is made of 2 gametes (diploidy)
    RECOMB <- RECOMB * 2
    
    # Counts gametes
    for (i.genotype in 1 : init.par$m) {
      gamete_counts <- gamete_counts + N[i.genotype] * RECOMB[,i.genotype]
    }

    
  } else {
    # This is for gamete based
    
    # For each genotype, get names and the "position" of the two gametes forming the genotype
    if( any(grep("/",names(N)) != c(1:init.par$m)) ) stop("Unable to find the separating slash in genotype names in N")
    gametes_per_genotype <- strsplit(names(N),"/") # Names of the two gametes forming the genotype
    names(gametes_per_genotype) <- names(N)
    corr_genotypes_gametes <- lapply(gametes_per_genotype,match,gamete_names) # Position of the two gametes forming the genotype, relative to the order of the gametes in gamete_names
    
    # Counts gametes
    # Add N to the first gamete. Add N to the second gamete.
    # Need to do this separately for the first and second gamete
    # Otherwise it does not count well the homozygous genotypes
    for (i.genotype in 1 : init.par$m) {
      gamete_counts[corr_genotypes_gametes[[i.genotype]][1]] <- gamete_counts[corr_genotypes_gametes[[i.genotype]][1]] + N[i.genotype]
      gamete_counts[corr_genotypes_gametes[[i.genotype]][2]] <- gamete_counts[corr_genotypes_gametes[[i.genotype]][2]] + N[i.genotype]
    }
  }
  
  # Gamete frequencies
  gamete_frequencies <- gamete_counts / sum(gamete_counts)
  
  # Out
  out <- list(counts = gamete_counts,
              frequencies = gamete_frequencies)
  out
  
}
