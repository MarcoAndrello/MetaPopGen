fst_multilocus_WC84 <- function(N, init.par) {
    if(length(dim(N))!=2) stop("N must have two dimensions: genotype and group")
    nLoc <- length(init.par$allele_vec) # Number of loci
    num_groups <- dim(N)[2]             # Number of groups
    locus_names <- paste0("Locus",LETTERS[1:nLoc])
    
    # Number of individuals for each group
    n_group <- colSums(N) # number of individuals for each group
    
    # Calculate allele frequency per "sample" (group)
    p_tilde_i <- array(NA,c(nLoc,num_groups))
    for (i.group in 1 : num_groups) {
        allele_freq_i <- freq_alleles(N[,i.group], init.par)$frequencies
        for (i.locus in 1 : nLoc) {
            p_tilde_i[i.locus,i.group] <- allele_freq_i[[i.locus]][1] # biallelic loci only
        }
    }
    
    # Observed heterozygosity
    h_tilde_i <- apply(N,2,het_obs,init.par)
    
    # Average sample size
    n_bar <- mean(n_group)
    
    # nc
    n_c <- (num_groups * n_bar - sum(n_group^2) / (num_groups * n_bar)) / (num_groups - 1)
    
    # Average sample frequency of the allele
    # Sample variance of allele frequencies
    # Average heterozygote frequency
    p_bar <- array(NA, nLoc)
    s_square <- array(NA, nLoc)
    h_bar <- array(NA, nLoc)
    for (i.locus in 1 : nLoc) {
        p_bar[i.locus] <- sum(n_group * p_tilde_i[i.locus,]) / (num_groups * n_bar)
        s_square[i.locus] <- sum(n_group * (p_tilde_i[i.locus,] - p_bar[i.locus])^2) / ((num_groups - 1)*n_bar)
        h_bar[i.locus] <- sum(n_group * h_tilde_i[i.locus,]) / (num_groups * n_bar)
    }
    
    # Components of the variance and fstWC84 per locus
    a_comp <- b_comp <- c_comp <- fst_WC84 <- array(NA,nLoc)
    for (i.locus in 1 : nLoc) {
        a_comp[i.locus] <- (n_bar/n_c) * (s_square[i.locus] - 1/(n_bar-1) * (p_bar[i.locus]*(1-p_bar[i.locus]) -
                                                                                 ((num_groups-1)/num_groups)*s_square[i.locus] - (h_bar[i.locus]/4)))
        b_comp[i.locus] <- (n_bar/(n_bar-1)) * (p_bar[i.locus]*(1-p_bar[i.locus]) - (num_groups-1)/num_groups*s_square[i.locus] - 
                                                    (2*n_bar-1)/(4*n_bar)*h_bar[i.locus])
        c_comp[i.locus] <- 0.5 * h_bar[i.locus]
        
        fst_WC84[i.locus] <- a_comp[i.locus] / (a_comp[i.locus] + b_comp[i.locus] + c_comp[i.locus])
    }
    
    return(fst_WC84)
}
