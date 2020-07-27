# Life-cycle functions
# Marco Andrello
# 06/02/2020

# Survival
# Adult dispersal
# Propagule dispersal
# Union gametes
# Reproduction
# Recruitment
# Regulate pop size
# Recruit settlers
# Recruit adults


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
    function(L,delta){
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
    }
    mat_geno
}

# Reproduction multilocus
repr <- function(Nprime_F, Nprime_M=NULL, phi_F, phi_M, l, m, z, Proba, mat_geno_to_index_mapping) {
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
            fecx[k,x] <- sum(as.numeric(rpois(Nprime_F[k,x],phi_F[k,x]))) # This is the contribution of variation in reproductive success among individuals to genetic drift
        }
        fec[k] <- sum(fecx[k,])
    }
    G_F <-Type_gamete(fec,Proba)
    
    # Calculate number of male gametes for each gametype
    fecx <- array(0,dim=c(m,z))	# Number of male gametes produced by all the individuals of each genotype in each age class
    fec  <- array(0,dim=m)		# Number of male gametes produced by all the individuals of each genotype
    for (k in 1 : m) {
        for (x in 1 : z) {
            fecx[k,x] <- sum(as.numeric(rpois(Nprime_M[k,x],phi_M[k,x]))) # This is the contribution of variation in reproductive success among individuals to genetic drift
        }
        fec[k] <- sum(fecx[k,])
    }
    G_M <- Type_gamete(fec,Proba)
    
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
