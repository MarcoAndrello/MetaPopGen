genotype.index = function(l){
  y <- array(NA,dim=c(l,l))
  diagy <- array(NA,dim=l)
  diagy[1] <- 1
  for (i in 2 : l) {
    diagy[i] <- diagy[i-1] + (l-(i-2))
  }
  diag(y)<-diagy
  for (i in 1 : (l-1)) {
    y[i,i:l] <- seq(y[i,i],(y[i,i]+l-i))
    y[i:l,i] <- seq(y[i,i],(y[i,i]+l-i))
  }
  
  # Set dimension names
  dn <- vector("character",l)
  for (i in 1 : l){
    dn[i] <- paste0("A",i)
  }
  dimnames(y) <- list(dn,dn)
  
  return(y)
}