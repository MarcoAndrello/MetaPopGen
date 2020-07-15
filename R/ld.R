# Linkage disequilibrium
ld <- function(N, init.par){
  if(length(init.par$allele_vec)!=2) stop(paste0("Detected ",length(init.par$allele_vec)," loci. Only available for two loci"))
  if(any(init.par$allele_vec!=c(2,2))) stop("Detected more than two alleles per locus. only available for two alleles per locus")
  
  fgam <- freq_gametes(N, init.par)$frequencies
  ld <- fgam[1]*fgam[4] - fgam[2]*fgam[3]
  names(ld) <- "D"
  # # Alternative calculation
  # fA <- freq_alleles(N, init.par)$frequencies
  # fgam[1] - fA$LocusA[1]*fA$LocusB[1]
  ld
}
