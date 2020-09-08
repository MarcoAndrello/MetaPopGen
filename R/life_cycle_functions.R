# Life-cycle functions
# Marco Andrello
# 06/02/2020

# Survival
# Adult dispersal
# Propagule dispersal
# Union gametes
# Reproduction
# Recruitment with backward migration
# Recruitment
# Regulate pop size
# Recruit settlers
# Recruit adults
# Recruitment with backward migration one locus
# Reproduction monoecious one locus
# Reproduction dioecious one locus


# Survival
surv <-
  function(sigma,N) {
    rbinom(1,N,sigma)
  }

# Adult dispersal
# N.deme.age = Nprime[k,i,x]
# delta.ad.deme = delta.ad[,i,x,t]
disp.ad <-
  function(N.deme.age,delta.ad.deme){
    delta_lost <- max(0,1 - sum(delta.ad.deme)) # The maximum function is needed to avoid errors due to precision
    delta_add <- c(delta.ad.deme,delta_lost)
    as.vector(rmultinom(1,N.deme.age,delta_add))
  }

# Dispersal
disp <-
  function(L, delta){
    delta_lost <- max(0,1 - sum(delta)) # The maximum function is needed to avoid errors due to precision
    delta_add <- c(delta,delta_lost)
    S <- rmultinom(1,L,delta_add)
    return(S)
  }


# Union gametes
union.gametes <- function(G_M, G_F, l) {
  if (sum(G_F) <= sum(G_M)) {
    G_max <- G_M
    G_min <- G_F
  } else {
    G_max <- G_F
    G_min <- G_M
  }
  mat_geno <- array(0,dim=c(l,l))
  Gprime_max <- G_max
  for (j in 1 : l) {
    in_dist <- Gprime_max 
    odds    <- array(1,dim=l)
    ndraws  <- G_min[j]
    err1 <- try(rMWNCHypergeo(1,in_dist,ndraws,odds),silent=T)
    if (class(err1)=="try-error") {
      err2 <- try(as.numeric(rmultinom(1,ndraws,in_dist)),silent=T)
      if (class(err2)=="try-error") {
        # Use multivariate normal
        prob      <- in_dist/sum(in_dist)
        mu.mvr    <- ndraws * prob  			    # Vector of means of the multivariate normal distribution
        var.mvr   <- ndraws * prob * (1-prob)	    # Vector of variances of the multivariate normal distribution
        sigma.mvr <- diag(var.mvr, l)			    # Variance-covariance matrix of the multivariate normal distribution
        for (i.mvr in 1 : l){
          for (j.mvr in 1 : l) {
            if (i.mvr == j.mvr) next
            sigma.mvr[i.mvr,j.mvr] <- -ndraws * prob[i.mvr] * prob[j.mvr]
          }
        }
        extr <- as.vector(round(loinorm(1,mu.mvr,sigma.mvr)))
      } else {
        extr<- err2                                  # Use multinomial
      }
    } else {
      extr <- err1
    }
    mat_geno[j,] <- extr
    Gprime_max <- Gprime_max - extr
    Gprime_max[which(Gprime_max<0)] <- 0 # Needed bcs multinomial and multivariate normal could extract more gametes than there are
  }
  mat_geno
}

# Reproduction multilocus
repr <- function(Nprime_F, Nprime_M=NULL, phi_F, phi_M, l, m, z, Proba,
                 mat_geno_to_index_mapping, fec.distr_F, fec.distr_M, migration) {
  # Monoecious: N_prime_M will be NULL: set it to Nprime_F 
  if(is.null(Nprime_M)) {
    Nprime_M <- Nprime_F
    sexuality <- "monoecious"
  } else {
    sexuality <- "dioecious"
  }
  # Force age dimension when there is only one age class
  # This was added on 29/10/2019 to fix case of n demes, 1 age class.
  # Think of how to fix the general case (1|n) demes with (1|n) age-class
  if (is.null(dim(Nprime_F))) dim(Nprime_F) <- c(length(Nprime_F),1)
  if (is.null(dim(Nprime_M))) dim(Nprime_M) <- c(length(Nprime_M),1)
  if (is.null(dim(phi_F))) dim(phi_F) <- c(length(phi_F),1)
  if (is.null(dim(phi_M))) dim(phi_M) <- c(length(phi_M),1)
  
  # Calculate number of female gametes for each gametype
  fecx <- array(0,dim=c(m,z))	# Number of female gametes produced by all the individuals of each genotype in each age class
  fec  <- array(0,dim=m)	    # Number of female gametes produced by all the individuals of each genotype
  for (k in 1 : m) {
    for (x in 1 : z) {
      if (fec.distr_F == "poisson") {
        fecx[k,x] <- sum(as.numeric(rpois(Nprime_F[k,x],phi_F[k,x]))) # This is the contribution of variation in reproductive success among individuals to genetic drift
      } else {
        fecx[k,x] <- sum(as.numeric(Nprime_F[k,x]*phi_F[k,x])) 
      }
    }
    fec[k] <- sum(fecx[k,])
  }
  G_F <-Type_gamete(fec,Proba)
  
  # Calculate number of male gametes for each gametype
  fecx <- array(0,dim=c(m,z))	# Number of male gametes produced by all the individuals of each genotype in each age class
  fec  <- array(0,dim=m)		# Number of male gametes produced by all the individuals of each genotype
  for (k in 1 : m) {
    for (x in 1 : z) {
      if (fec.distr_M == "poisson") {
        fecx[k,x] <- sum(as.numeric(rpois(Nprime_M[k,x],phi_M[k,x]))) # This is the contribution of variation in reproductive success among individuals to genetic drift
      } else {
        fecx[k,x] <- sum(as.numeric(Nprime_M[k,x],phi_M[k,x]))
      }
    }
    fec[k] <- sum(fecx[k,])
  }
  G_M <- Type_gamete(fec,Proba)
  
  if(migration == "backward") return(list(G_F=G_F, G_M=G_M))
  
  # Union of gametes to form zygotes
  mat_geno <- union.gametes(G_M, G_F, l)
  
  # Calculate genotype numbers in propagules L by adding the combinations of gametotypes
  # to see which gametotype is which: rownames(mat_geno) <- colnames(mat_geno) <- rownames(Proba)
  L <- rep(0,m)
  for (i.L in 1 : m) {
    L[i.L] <- L[i.L] + sum( mat_geno[which(mat_geno_to_index_mapping==i.L)] )
  }
  if (sexuality == "monoecious") {
    return(L)
  } else {
    # Sex differentiation is independent of genotype, with equal proportion of male and females produced
    LL <- array(NA,dim=c(m,2))
    for (i in 1 : m){
      LL[i,1] <- rbinom(1,L[i],0.5)
      LL[i,2] <- L[i] - LL[i,1] 
    }
    return(LL)
  }
}

# Recruitment with backward migration
recr.backward.migration <- function(migr, l, m, n, z, j.local.deme, kappa0,
                                    sexuality, mat_geno_to_index_mapping, t, 
                                    gamete_pool) {
  # G_F is genotype * deme
  # G_M is genotype * deme
  if (sexuality == "dioecious") {
    stop("Backward migration not implemented yet for dioecious life cycles")
  }
  if (z == 1) {
    # does not need to consider adults
    # can recruit *gametes*
    Ntot_recruits <- kappa0
    
    # # Build local gamete pool: "exact migration" method
    # G_F_local <- array(0,l)
    # for (j in 1 : n) {
    #   if (j == j.local.deme) {
    #     num_gametes <- (1-migr)*kappa0 
    #   } else {
    #     num_gametes <- migr/(n-1)*kappa0
    #   }
    #   if (num_gametes != round(num_gametes)) warning("num female gametes at deme ",j.local.deme," from deme ",j," at time ",t," is not an integer")
    #   G_F_local <- G_F_local + as.vector(rmultinom(1,num_gametes,G_F[,j]))
    #   # cat("deme",j.local.deme,"from deme",j,"done. Local female gametes:",sum(G_F_local),"\n")
    # }
    # G_M_local <- array(0,l)
    # for (j in 1 : n) {
    #   if (j == j.local.deme) {
    #     num_gametes <- (1-migr)*kappa0
    #     } else {
    #       num_gametes <- migr/(n-1)*kappa0
    #     }
    #   if (num_gametes != round(num_gametes)) warning("num male gametes at deme ",j.local.deme," from deme ",j," at time ",t," is not an integer")
    #   G_M_local <- G_M_local + as.vector(rmultinom(1,num_gametes,G_M[,j]))
    #   # cat("deme",j.local.deme,"from deme",j,"done. Local male gametes:",sum(G_M_local),"\n")
    # }
    
    # Build local gamete pool: "sampling migration" method
    # Calculates immigration probabilities for this deme
    immigr.prob <- rep(migr/(n-1),n)
    immigr.prob[j.local.deme] <- (1-migr)
    # Calculates gamete frequencies "seen" from this deme and sample gametes from them
    # freq_G_F <- prop.table(G_F,2)
    # G_F_j <- freq_G_F %*% immigr.prob
    # G_F_local <- as.vector(rmultinom(1, Ntot_recruits, G_F_j))
    # freq_G_M <- prop.table(G_M,2)
    # G_M_j <- freq_G_M %*% immigr.prob
    # G_M_local <- as.vector(rmultinom(1, Ntot_recruits, G_M_j))
    freq_G <- prop.table(t(gamete_pool),2)
    G_j <- freq_G %*% immigr.prob
    G_F_local <- as.vector(rmultinom(1, Ntot_recruits, G_j))
    G_M_local <- as.vector(rmultinom(1, Ntot_recruits, G_j))
    
    # Then, union gametes to form zygotes
    mat_geno <- union.gametes(G_M_local, G_F_local, l)
    
    # Calculate genotype numbers in new recruits (Naged) by adding the combinations of gametotypes
    # to see which gametotype is which: rownames(mat_geno) <- colnames(mat_geno) <- rownames(Proba)
    Naged <- rep(0,m)
    for (i.L in 1 : m) {
      Naged[i.L] <- Naged[i.L] + sum( mat_geno[which(mat_geno_to_index_mapping==i.L)] )
    }
    if (sexuality == "monoecious") {
      return(Naged)
    } else {
      stop("Backward migration not implemented for dioecious life cycles")
      # See in function repr when you want to implement it
    }
    # dimnames(Naged) <- dimnames(N)
  } else {
   stop("Backward migration not implemented yet for number of age classes > 1")
  }
}


# Recruitment
recr <- function(N, N_F, N_M, S, S_F, S_M, m, z, kappa0, recr.dd, sexuality) {
  if (sexuality == "dioecious") {
    N <- abind(N_F,N_M,along=3)
    # dimnames(N)[[3]] <- c("females","males")
    # names(dimnames(N)) <- c("genotype","age","sex")
    
    S <- abind(S_F,S_M,along=3)
    # dimnames(S) <- dimnames(S_F)
    # dimnames(S)[[3]] <- c("females","males")
    # names(dimnames(S)) <- c("genotype","age","sex")
  }
  # Only offspring compete
  if (z == 1) {
    # does not need to consider adults
    # can recruit settlers
    Ntot_recruits <- kappa0
    sett_reg <- recruit.settlers(S, Ntot_recruits)
    # Build aged population
    Naged <- sett_reg
    # dimnames(Naged) <- dimnames(N)
  } else {
    # Considering only adults from age class 2
    if (sexuality == "dioecious") {
      Nadlt <- N[,1:(z-1),,drop=F] 
    } else {
      Nadlt <- N[,1:(z-1),drop=F]
    }
    adlt_reg <- recruit.adults(Nadlt, kappa0)
    Ntot_recruits <- kappa0 - sum(adlt_reg)
    sett_reg <- recruit.settlers(S, Ntot_recruits)
    # Build aged population
    Naged <- abind(sett_reg,adlt_reg,along=2)
    # dimnames(Naged) <- dimnames(N)
  }
  return(Naged)
}


regulate.pop.size <- function(x,k) {
  xvec <- as.vector(x)
  # Perform a multivariate hypergeometric draw
  err1 <- try(rMWNCHypergeo(1,
                            xvec,
                            k,
                            odds=rep(1,length(xvec))),
              silent = T)
  if (class(err1)=="try-error") {
    err2 <- try(as.numeric(rmultinom(1,
                                     k,
                                     xvec)),
                silent=T)
    if (class(err2)=="try-error") { 
      # use multivariate normal
    }
    #print("Using multinomial")
    yvec <- err2
  } else {
    #print("Using multivariate hypergeometric")
    yvec <- err1
  }
  y <- array(yvec,dim=dim(x))
  return(y)
}


recruit.settlers <- function(S, Ntot_recruits) {
  if (Ntot_recruits == 0) {
    #print("cannot recruit settlers")
    sett_reg <- array(0,
                      dim=dim(S))
    return(sett_reg)
  }
  if (sum(S) > Ntot_recruits) {
    #print("need to regulate settlers")
    sett_reg <- regulate.pop.size(S, Ntot_recruits)
    return(sett_reg)
  }
  #print("does not need to regulate settlers")
  sett_reg <- S
  return(sett_reg)
}


recruit.adults <- function(Nadlt, kappa0) {
  if (sum(Nadlt) > kappa0) {
    #print("# need to regulate adults")
    adlt_reg <- regulate.pop.size(Nadlt, kappa0)
    # dimnames(adlt_reg) <- dimnames(Nadlt)
    return(adlt_reg)
  }
  #print("# does not need to regulate adults")
  adlt_reg <- Nadlt
  return(adlt_reg)
}



recr.backward.migration.onelocus <- function(G_F, G_M, migr, l, m, n, z, j.local.deme, kappa0, sexuality, t) {
  # G_F is genotype * deme
  # G_M is genotype * deme
  if (sexuality == "dioecious") {
    stop("Backward migration not implemented yet for dioecious life cycles")
  }
  # z == 1
  # does not need to consider adults
  # can recruit *gametes*
  Ntot_recruits <- kappa0
  
  # Build local gamete pool: "sampling migration" method
  # Calculates immigration probabilities for this deme
  immigr.prob <- rep(migr/(n-1),n)
  immigr.prob[j.local.deme] <- (1-migr)
  # Calculates gamete frequencies "seen" from this deme and sample gametes from them
  freq_G_F <- prop.table(G_F,2)
  G_F_j <- freq_G_F %*% immigr.prob
  G_F_local <- as.vector(rmultinom(1, Ntot_recruits, G_F_j))
  freq_G_M <- prop.table(G_M,2)
  G_M_j <- freq_G_M %*% immigr.prob
  G_M_local <- as.vector(rmultinom(1, Ntot_recruits, G_M_j))
  
  
  # Then, union gametes
  # Union of gametes to form zygotes
  
  if (sum(G_F_local) <= sum(G_M_local)) {
    
    mat_geno<- array(0,dim=c(l,l))
    Gprime_M<- G_M_local
    for (j in 1 : l) {
      in_dist<- Gprime_M 
      odds<- array(1,dim=l)
      ndraws<- G_F_local[j]
      err<- try(rMWNCHypergeo(1,in_dist,ndraws,odds),silent=T)
      if (class(err)=="try-error") {
        extr<- as.numeric(rmultinom(1,ndraws,in_dist))
      } else {
        extr<- err
      }
      mat_geno[j,]<- extr
      Gprime_M<- Gprime_M - extr
      Gprime_M[which(Gprime_M<0)] <- 0 # Needed bcs multinomial and multivariate normal could extract more gametes than there are
    }
    mat_geno_l<- mat_geno
    mat_geno_l[upper.tri(mat_geno_l)]<- 0
    mat_geno_u<- mat_geno
    mat_geno_u[lower.tri(mat_geno_u,diag=T)]<- 0
    mat_geno_f<- mat_geno_l + t(mat_geno_u)
    L<- mat_geno_f[lower.tri(mat_geno_f,diag=T)]
    
  } else {
    
    mat_geno<- array(0,dim=c(l,l))
    Gprime_F<- G_F_local
    for (j in 1 : l) {
      in_dist<- Gprime_F 
      odds<- array(1,dim=l)
      ndraws<- G_M_local[j]
      err<- try(rMWNCHypergeo(1,in_dist,ndraws,odds),silent=T)
      if (class(err)=="try-error") {
        extr<- as.numeric(rmultinom(1,ndraws,in_dist))
      } else {
        extr<- err
      }
      mat_geno[j,]<- extr
      Gprime_F<- Gprime_F - extr
      Gprime_F[which(Gprime_F<0)] <- 0 # Needed bcs multinomial and multivariate normal could extract more gametes than there are
    }
    mat_geno_l<- mat_geno
    mat_geno_l[upper.tri(mat_geno_l)]<- 0
    mat_geno_u<- mat_geno
    mat_geno_u[lower.tri(mat_geno_u,diag=T)]<- 0
    mat_geno_f<- mat_geno_l + t(mat_geno_u)
    Naged <- mat_geno_f[lower.tri(mat_geno_f,diag=T)]
  }
}


# Reproduction monoecious one locus
repr.monoecious.onelocus <-
  function(Nprime, phi_F, phi_M,
           mu, i, l, m, n, z,
           fec.distr_F, fec.distr_M, migration, verbose) {
    
    # Adjust dimensions if there is only one age-class and/or one deme
    #if (length(dim(phi_F))==2) dim(phi_F)[3]<-1
    #if (length(dim(phi_M))==2) dim(phi_M)[3]<-1
    dim(phi_F) <- c(m,n,z)
    dim(phi_M) <- c(m,n,z)
    # n and z are not passed to the function? But it works...
    
    if (verbose) print("Calculates total number of female gametes")
    fecx<- array(0,dim=c(m,z))	# Number of female gametes produced by all the individuals of each genotype in each age class
    fec <- array(0,dim=m)			  # Number of female gametes produced by all the individuals of each genotype
    for (k in 1 : m) {
      for (x in 1 : z) {
        if (fec.distr_F == "poisson") {
          fecx[k,x] <- sum(as.numeric(rpois(Nprime[k,i,x],phi_F[k,i,x])))	# This is the contribution of variation in reproductive success among individuals to genetic drift
        } else {
          fecx[k,x] <- sum(as.numeric(Nprime[k,i,x]*phi_F[k,i,x]))
        }
      }
      fec[k] <- sum(fecx[k,])
    }
    if (verbose) print("Calculate number of gametes for each allele")
    G_F <- array(0,dim=l)
    k <- 1
    for (j in 1 : l) {
      for (jj in j : l) {
        if (j == jj) {
          G_F[j] <- G_F[j] + fec[k]
        } else {
          meiosis_j<- rbinom(1,fec[k],0.5)# This is the contribution of Mendelian segregation to genetic drift
          meiosis_jj<- fec[k] - meiosis_j
          G_F[j]<- G_F[j] + meiosis_j
          G_F[jj]<- G_F[jj] + meiosis_jj
        }
        k <- k + 1
      }
    }
    
    
    if (verbose) print("Calculates total number of male gametes")
    fecx <- array(0,dim=c(m,z))	# Number of male gametes produced by all the individuals of each genotype in each age class
    fec <- array(0,dim=m)			# Number of male gametes produced by all the individuals of each genotype
    for (k in 1 : m) {
      for (x in 1 : z) {
        if (fec.distr_F == "poisson") {
          fecx[k,x] <- sum(as.numeric(rpois(Nprime[k,i,x],phi_M[k,i,x])))	# This is the contribution of variation in reproductive success among individuals to genetic drift
        } else {
          fecx[k,x] <- sum(as.numeric(Nprime[k,i,x]*phi_M[k,i,x]))
        }
      }
      fec[k] <- sum(fecx[k,])
    }
    if (verbose) print("Calculate number of gametes for each allele")
    G_M <- array(0,dim=l)
    k <- 1
    for (j in 1 : l) {
      for (jj in j : l) {
        if (j == jj) {
          G_M[j] <- G_M[j] + fec[k]
        } else {
          meiosis_j<- rbinom(1,fec[k],0.5)
          meiosis_jj<- fec[k] - meiosis_j
          G_M[j]<- G_M[j] + meiosis_j
          G_M[jj]<- G_M[jj] + meiosis_jj
        }
        k <- k + 1
      }
    }
    
    if (verbose) print("Mutation")
    Gprime_F <- array(0,dim=l)
    Gprime_M <- array(0,dim=l)
    for (j in 1 : l){
      Gprime_F <- Gprime_F + as.vector(rmultinom(1,G_F[j],mu[,j]))
      Gprime_M <- Gprime_M + as.vector(rmultinom(1,G_M[j],mu[,j]))
    }
    G_F <- Gprime_F
    G_M <- Gprime_M
    
    if(migration == "backward") return(list(G_F=G_F, G_M=G_M))
    
    if (verbose) print("Union of gametes to form zygotes")
    if (sum(G_F) <= sum(G_M)) {
      
      mat_geno<- array(0,dim=c(l,l))
      Gprime_M<- G_M
      for (j in 1 : l) {
        in_dist<- Gprime_M 
        odds<- array(1,dim=l)
        ndraws<- G_F[j]
        err<- try(rMWNCHypergeo(1,in_dist,ndraws,odds),silent=T)
        if (class(err)=="try-error") {
          extr<- as.numeric(rmultinom(1,ndraws,in_dist))
        } else {
          extr<- err
        }
        mat_geno[j,]<- extr
        Gprime_M<- Gprime_M - extr
      }
      mat_geno_l<- mat_geno
      mat_geno_l[upper.tri(mat_geno_l)]<- 0
      mat_geno_u<- mat_geno
      mat_geno_u[lower.tri(mat_geno_u,diag=T)]<- 0
      mat_geno_f<- mat_geno_l + t(mat_geno_u)
      L<- mat_geno_f[lower.tri(mat_geno_f,diag=T)]
      
    } else {
      
      mat_geno<- array(0,dim=c(l,l))
      Gprime_F<- G_F
      for (j in 1 : l) {
        in_dist<- Gprime_F 
        odds<- array(1,dim=l)
        ndraws<- G_M[j]
        err<- try(rMWNCHypergeo(1,in_dist,ndraws,odds),silent=T)
        if (class(err)=="try-error") {
          extr<- as.numeric(rmultinom(1,ndraws,in_dist))
        } else {
          extr<- err
        }
        mat_geno[j,]<- extr
        Gprime_F<- Gprime_F - extr
      }
      mat_geno_l<- mat_geno
      mat_geno_l[upper.tri(mat_geno_l)]<- 0
      mat_geno_u<- mat_geno
      mat_geno_u[lower.tri(mat_geno_u,diag=T)]<- 0
      mat_geno_f<- mat_geno_l + t(mat_geno_u)
      L<- mat_geno_f[lower.tri(mat_geno_f,diag=T)]
      
    }
    
    return(L)
  }


# Reproduction dioecious one locus
repr.dioecious.onelocus <-
  function(Nprime_M, Nprime_F, phi_F, phi_M,
           mu, i, l, m, n, z,
           fec.distr_F, fec.distr_M, migration, verbose) {
    
    # Adjust dimensions if there is only one age-class and/or one deme
    #if (length(dim(phi_F))==2) dim(phi_F)[3]<-1
    #if (length(dim(phi_M))==2) dim(phi_M)[3]<-1
    dim(phi_F) <- c(m,n,z)
    dim(phi_M) <- c(m,n,z)
    
    # Adjust dimensions if there is only one age-class
    if (length(dim(phi_F))==2) dim(phi_F)[3]<-1
    if (length(dim(phi_M))==2) dim(phi_M)[3]<-1
    
    
    
    ## Female gametes
    
    if (verbose) print("Calculates total number of female gametes")
    
    fecx <- array(0,dim=c(m,z))	# Number of female gametes produced by all the individuals of each genotype in each age class
    fec <- array(0,dim=m)			# Number of female gametes produced by all the individuals of each genotype
    
    for (k in 1 : m) {
      for (x in 1 : z) {
        if (fec.distr_F == "poisson") {
          fecx[k,x] <- sum(as.numeric(rpois(Nprime_F[k,i,x], phi_F[k,i,x])))	# This is the contribution of variation in reproductive success among individuals to genetic drift
        } else {
          fecx[k,x] <- sum(as.numeric(Nprime_F[k,i,x]*phi_F[k,i,x]))
        }
      }
      fec[k] <- sum(fecx[k,])
    }
    
    
    if (verbose) print("Calculates number of gametes for each allele")
    G_F <- array(0,dim=l)
    k <- 1
    for (j in 1 : l) {    # First allele
      for (jj in j : l) { # Second allele; looping in this way reproduces the genotype order: A1A1,A1A2,A1A3,A2A2,A2A3,A3A3
        if (j == jj) {
          G_F[j] <- G_F[j] + fec[k]
        } else {
          
          meiosis_j <- rbinom(1,fec[k],0.5) 	# This is the contribution of Mendelian segregation to genetic drift
          if (is.na(meiosis_j)) meiosis_j <- round(rnorm(1,(0.5*fec[k]),sqrt(0.25*fec[k]))) 	# If fec[k] is too large, uses the normal approximation to the binomial distribution
          
          meiosis_jj<- fec[k] - meiosis_j
          G_F[j]<- G_F[j] + meiosis_j
          G_F[jj]<- G_F[jj] + meiosis_jj
        }
        k <- k + 1
      }
    }
    
    ## Male gametes
    if (verbose) print("Calculates total number of gametes")
    
    fecx <- array(0,dim=c(m,z))	# Number of male gametes produced by all the individuals of each genotype in each age class
    fec <- array(0,dim=m)			# Number of male gametes produced by all the individuals of each genotype
    for (k in 1 : m) {
      for (x in 1 : z) {
        if (fec.distr_M == "poisson") {
          fecx[k,x] <- sum(as.numeric(rpois(Nprime_M[k,i,x], phi_M[k,i,x])))	# This is the contribution of variation in reproductive success among individuals to genetic drift
        } else {
          fecx[k,x] <- sum(as.numeric(Nprime_M[k,i,x]*phi_M[k,i,x]))
        }
      }
      fec[k] <- sum(fecx[k,])
    }
    
    if (verbose) print("Calculates number of gametes for each allele")
    G_M <- array(0,dim=l)
    k <- 1
    for (j in 1 : l) {
      for (jj in j : l) {
        if (j == jj) {
          G_M[j] <- G_M[j] + fec[k]
        } else {
          meiosis_j<- rbinom(1,fec[k],0.5)
          meiosis_jj<- fec[k] - meiosis_j
          G_M[j]<- G_M[j] + meiosis_j
          G_M[jj]<- G_M[jj] + meiosis_jj
        }
        k <- k + 1
      }
    }
    
    
    if (verbose) print("Mutation")
    Gprime_F <- array(0,dim=l)
    Gprime_M <- array(0,dim=l)
    
    # Female gametes  
    for (j in 1 : l){
      if (G_F[j] == 0) next
      
      
      err <- try(as.vector(rmultinom(1,G_F[j],mu[,j])),silent=T)
      
      if (class(err)=="try-error") {
        mu.mvr <- G_F[j] * mu[,j]				# Vector of means of the multivariate normal distribution
        var.mvr <- G_F[j]* mu[,j] * (1-mu[,j])	# Vector of variances of the multivariate normal distribution
        sigma.mvr <- diag(var.mvr, l)			# Variance-covariance matrix of the multivariate normal distribution
        
        for (i.mvr in 1 : l){
          for (j.mvr in 1 : l) {
            if (i.mvr == j.mvr) next
            sigma.mvr[i.mvr,j.mvr] <- -G_F[j] * mu[i.mvr,j] * mu[j.mvr,j]
          }
        }
        
        Gprime_F <- Gprime_F + as.vector(round(mvrnorm(1,mu.mvr,sigma.mvr)))
        
      } else {
        Gprime_F <- Gprime_F + err
      }
    }
    
    # Male gametes  
    for (j in 1 : l){
      if (G_M[j] == 0) next
      
      err <- try(as.vector(rmultinom(1,G_M[j],mu[,j])),silent=T)
      
      if (class(err)=="try-error") {
        mu.mvr <- G_M[j] * mu[,j]				# Vector of means of the multivariate normal distribution
        var.mvr <- G_M[j]* mu[,j] * (1-mu[,j])	# Vector of variances of the multivariate normal distribution
        sigma.mvr <- diag(var.mvr, l)			# Variance-covariance matrix of the multivariate normal distribution
        
        for (i.mvr in 1 : l){
          for (j.mvr in 1 : l) {
            if (i.mvr == j.mvr) next
            sigma.mvr[i.mvr,j.mvr] <- -G_M[j] * mu[i.mvr,j] * mu[j.mvr,j]
          }
        }
        
        Gprime_M <- Gprime_M + as.vector(round(mvrnorm(1,mu.mvr,sigma.mvr)))
        
      } else {
        Gprime_M <- Gprime_M + err
      }
      
    } 
    
    G_F <- Gprime_F
    G_M <- Gprime_M
    
    
    if (verbose) print("Union of gametes to form zygotes")
    if (sum(G_F) <= sum(G_M)) {
      
      mat_geno<- array(0,dim=c(l,l))
      Gprime_M<- G_M
      for (j in 1 : l) {
        in_dist<- Gprime_M 
        odds<- array(1,dim=l)
        ndraws<- G_F[j]
        err<- try(rMWNCHypergeo(1,in_dist,ndraws,odds),silent=T)
        if (class(err)=="try-error") {
          extr<- as.numeric(rmultinom(1,ndraws,in_dist))
        } else {
          extr<- err
        }
        mat_geno[j,]<- extr
        Gprime_M<- Gprime_M - extr
      }
      mat_geno_l<- mat_geno
      mat_geno_l[upper.tri(mat_geno_l)]<- 0
      mat_geno_u<- mat_geno
      mat_geno_u[lower.tri(mat_geno_u,diag=T)]<- 0
      mat_geno_f<- mat_geno_l + t(mat_geno_u)
      L <- mat_geno_f[lower.tri(mat_geno_f,diag=T)]
      
    } else {
      
      
      mat_geno <- array(0,dim=c(l,l))
      Gprime_F <- G_F
      for (j in 1 : l) {
        in_dist <- Gprime_F 
        odds <- array(1,dim=l)
        ndraws <- G_M[j]
        # If the mutlivariate hypergeometric does not work, use the multinomial
        # If the multinomial does not work, use the multivariate normal
        
        err1 <- try(rMWNCHypergeo(1,in_dist,ndraws,odds),silent=T)
        
        if (class(err1)=="try-error") {
          
          err2 <- try(as.numeric(rmultinom(1,ndraws,in_dist)),silent=T)
          
          if (class(err2)=="try-error") {
            
            # Use multivariate normal
            prob <- in_dist/sum(in_dist)
            
            mu.mvr <- ndraws * prob  			               # Vector of means of the multivariate normal distribution
            var.mvr <- ndraws * prob * (1-prob)	       # Vector of variances of the multivariate normal distribution
            sigma.mvr <- diag(var.mvr, l)			             # Variance-covariance matrix of the multivariate normal distribution
            
            for (i.mvr in 1 : l){
              for (j.mvr in 1 : l) {
                if (i.mvr == j.mvr) next
                sigma.mvr[i.mvr,j.mvr] <- -ndraws * prob[i.mvr] * prob[j.mvr]
              }
            }
            
            extr <- as.vector(round(mvrnorm(1,mu.mvr,sigma.mvr)))
            
          } else {
            extr<- err2                                     # Use multinomial
          }
          
        } else {
          extr <- err1                                       # Use multivariate hypergeometric
        }
        
        mat_geno[j,]<- extr
        Gprime_F<- Gprime_F - extr
      }
      mat_geno_l<- mat_geno
      mat_geno_l[upper.tri(mat_geno_l)]<- 0
      mat_geno_u<- mat_geno
      mat_geno_u[lower.tri(mat_geno_u,diag=T)]<- 0
      mat_geno_f<- mat_geno_l + t(mat_geno_u)
      L<- mat_geno_f[lower.tri(mat_geno_f,diag=T)]
      
    }
    
    if (verbose) print("Determine sex")
    # Sex differentiation is independent of genotype, with equal proportion of male and females produced
    LL <- array(NA,dim=c(m,2))
    for (i in 1 : m){
      LL[i,1] <- rbinom(1,L[i],0.5)
      LL[i,2] <- L[i] - LL[i,1] 
    }
    
    return(LL)
  }

