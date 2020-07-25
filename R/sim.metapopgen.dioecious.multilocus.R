sim.metapopgen.dioecious.multilocus <- function(init.par,
                                                sigma_F, sigma_M,
                                                phi_F, phi_M,
                                                delta.prop = NULL, delta.ad = NULL,
                                                recr.dd="settlers",
                                                T_max,
                                                save.res=F, save.res.T=seq(1:T_max),
                                                verbose=F) {
  
  # Reading basic variables
  m                         <- init.par$m                          # Number of genotypes
  meiosis_matrix            <- init.par$meiosis_matrix             # Meiosis matrix
  mat_geno_to_index_mapping <- init.par$mat_geno_to_index_mapping  # Mapping of combinations of gametotypes to genotypes
  if (is.null(mat_geno_to_index_mapping)) mat_geno_to_index_mapping <- init.par$index_matr # This is for gamete-based representation
  l                         <- dim(meiosis_matrix)[1]              # Number of multilocus gametotypes
  n                         <- init.par$n                          # Number of demes
  z                         <- init.par$z                          # Number of age-classes
  kappa0                    <- init.par$kappa0                     # Carrying capacity
  N1_F                      <- init.par$N1_F
  N1_M                      <- init.par$N1_M
  
  # If only one age-class and recr.dd=="adults", gives an error
  if (z == 1 & recr.dd == "adults") {
    stop("Detected only one age class (z=1) and recruitment probability dependent on adult density (recr.dd == 'adults'). This combination is not supported. Use recr.dd == 'settlers' instead.")
  }
  
  # .............................Check the existence of dispersal matrices
  if (is.null(delta.prop)) {
    print("Setting propagule dispersal probability (delta.prop)")
    delta.prop <- diag(1,n)
    delta.prop <- array(delta.prop,c(n,n,T_max))
  }
  if (is.null(delta.ad)) {
    print("Setting adult dispersal probability (delta.ad)")
    delta.ad <- diag(1,n)
    delta.ad <- array(delta.ad,c(n,n,z,T_max))
  }
  
  
  ##########################################################################
  # Check if input data are time-dependent or not; in case, augment them
  ##########################################################################
  
  # Female Survival
  if (is.na(dim(sigma_F)[4])) {
    print("Augmenting sigma_F for time dimension")
    sigma_F <- array(rep(sigma_F,T_max),c(m,n,z,T_max))
  }
  
  # Male Survival
  if (is.na(dim(sigma_M)[4])) {
    print("Augmenting sigma_M for time dimension")
    sigma_M <- array(rep(sigma_M,T_max),c(m,n,z,T_max))
  }
  
  # Female fecundity
  if (is.na(dim(phi_F)[4])) {
    print("Augmenting phi_F for time dimension")
    phi_F <- array(rep(phi_F,T_max),c(m,n,z,T_max))
  }
  
  # Male fecundity
  if (is.na(dim(phi_M)[4])) {
    print("Augmenting phi_M for time dimension")
    phi_M <- array(rep(phi_M,T_max),c(m,n,z,T_max))
  }
  
  # Dispersal
  if (is.na(dim(delta.prop)[3])) {
    print("Augmenting delta.prop for time dimension")
    delta.prop <- array(rep(delta.prop,T_max),c(n,n,T_max))
  }
  
  # Adult dispersal
  if (is.na(dim(delta.ad)[3])) {
    print("Augmenting delta.ad for age dimension")
    delta.ad <- array(rep(delta.ad,z),c(n,n,z))
  }
  
  if (is.na(dim(delta.ad)[4])) {
    print("Augmenting delta.ad for time dimension")
    delta.ad <- array(rep(delta.ad,T_max),c(n,n,z,T_max))
  }
  
  # Carrying capacity
  if (is.vector(kappa0)) {
    print("Augmenting kappa0 for time dimension")
    kappa0 <- array(rep(kappa0,T_max),c(n,T_max))
  }
  
  ##########################################################################
  # Initialize state variables
  ##########################################################################
  print("Initializing variables...")
  if (save.res){
    N_F <- N1_F
    N_M <- N1_M
    dimnames(N_F) <- dimnames(N_M) <- dimnames(N1_M)
    rm(N1_M,N1_F)    
  } else {
    N_F       <- array(NA,dim=c(m,n,z,T_max))
    N_F[,,,1] <- N1_F
    N_M       <- array(NA,dim=c(m,n,z,T_max))
    N_M[,,,1] <- N1_M
    dimnamesN1 <- dimnames(N1_M)
    rm(N1_M,N1_F)
    dimnamesN1$time <- c(1:T_max)
    dimnames(N_F) <- dimnames(N_M) <- dimnamesN1
    L_F       <- array(NA,dim=c(m,n,T_max))
    L_M       <- array(NA,dim=c(m,n,T_max))
    S_F       <- array(0,dim=c(m,n,T_max))
    S_M       <- array(0,dim=c(m,n,T_max))
  }
  
  ##########################################################################
  # Define functions
  ##########################################################################
  
  # Recruitment:
  # S_M[,i],S_F[,i],Nprime_M[,i,],Nprime_F[,i,],m,kappa0[i,t]
  # S, S_othersex.  Number of settlers of all genotypes. Dimension: m
  # N_M, N_F.       Number of adults of all genotyps and age-classes. Dimensions: n*z
  # m.              Number of genotypes
  # kappa0.         Carrying capacity. Scalar.
  recr <-
    function(S,S_othersex,N_M,N_F,m,kappa0,recr.dd) {
      switch(recr.dd,
             # Dependence on settler density
             settlers = {
               Ntot <- sum(S+S_othersex)
               sigma0 <- settler.survival(Ntot,kappa0)
             },
             # Dependence on adult density
             adults = {
               if (z==1) Ntot <- 0 else Ntot <- sum(N_M[,1:(z-1)] + N_F[,1:(z-1)])
               Stot <- sum(S+S_othersex)
               Recr <- kappa0 - Ntot
               if (Recr <= 0){
                 sigma0 <- 0
               } else {
                 sigma0 <- Recr / Stot
               }
               if (sigma0 > 1) sigma0 <- 1
             },
             # No-match: error
             stop("Unknown value for argument recr.dd. Valid values: 'settlers', 'adults'")
      )
      
      # Use recruitment Probability to calculate the number of recruits
      R <- array(0,dim=m)
      for (k in 1 : m) {
        R[k] <- rbinom(1,S[k],sigma0)
      }
      return(R)
    }
  
  ##########################################################################
  # Simulate metapopulation genetics
  ##########################################################################
  if (save.res){
    dir.res.name <- paste(getwd(),format(Sys.time(), "%Y-%b-%d-%H.%M.%S"),sep="/")
    dir.create(dir.res.name)
    if (1 %in% save.res.T) {
      file.name <- "N1.RData"
      save(N_M,N_F,file=paste(dir.res.name,file.name,sep="/"))
    }
  }
  
  print("Running simulation...")
  for (t in 1 : (T_max-1)) {
    if (t %% 10 == 0) print(t)
    # m <- dim(N_M)[1]                 # Number of genotypes
    # l <- length(Proba [,1])          # Nuber of alleles
    # n <- dim(N_M)[2]                 # Number of demes
    # z <- dim(N_M)[3]                 # Number of age-classes 
    
    # At each time-step, redefine variable Nprime
    # If save.res, redefine also larval and settlers numbers
    if (save.res) {
      Nprime_F  <- array(NA,dim=c(m,n,z))
      Nprime_M  <- array(NA,dim=c(m,n,z))
      Nprimeprime_F  <- array(0,dim=c(m,n,z))
      Nprimeprime_M  <- array(0,dim=c(m,n,z))
      L_F       <- array(NA,dim=c(m,n))
      L_M       <- array(NA,dim=c(m,n))
      S_F       <- array(0,dim=c(m,n))
      S_M       <- array(0,dim=c(m,n))
    } else {
      Nprime_F  <- array(NA,dim=c(m,n,z))
      Nprime_M  <- array(NA,dim=c(m,n,z))
      Nprimeprime_F  <- array(0,dim=c(m,n,z))
      Nprimeprime_M  <- array(0,dim=c(m,n,z))
    }
    
    # Survival
    if (verbose) cat("t =",t,"Apply survival function \n")
    # If there is only one age-class, we must force the third dimension. What if only one year?
    if (length(dim(sigma_M))==2) dim(sigma_M)[3] <- 1
    if (length(dim(sigma_F))==2) dim(sigma_F)[3] <- 1
    for (i in 1 : n) {
      for (x in 1 : z) {
        for (k in 1 : m) {
          if (save.res){
            Nprime_M[k,i,x] = surv(sigma_M[k,i,x,t],N_M[k,i,x])
            Nprime_F[k,i,x] = surv(sigma_F[k,i,x,t],N_F[k,i,x])
          } else {
            Nprime_M[k,i,x] = surv(sigma_M[k,i,x,t],N_M[k,i,x,t])
            Nprime_F[k,i,x] = surv(sigma_F[k,i,x,t],N_F[k,i,x,t])
          }
          if(is.na(Nprime_M[k,i,x])){
            Nprime_M[k,i,x]<-0
            cat("NA adjustment in deme",i,"age",x,"genotype",k,"males \n")
          }
          if(is.na(Nprime_F[k,i,x])){
            Nprime_F[k,i,x]<-0
            cat("NA adjustment in deme",i,"age",x,"genotype",k,"females \n")
          }
        }
      }
    }  
    
    # Adult dispersal
    if (verbose) cat("t =",t,"Apply adult dispersal function \n") 
    for (i in 1 : n) {
      for (x in 1 : z) {
        for (k in 1 : m) {
          y = disp.ad(Nprime_F[k,i,x],delta.ad[,i,x,t])
          Nprimeprime_F[k,,x] <- Nprimeprime_F[k,,x] + y[1:n]
          y = disp.ad(Nprime_M[k,i,x],delta.ad[,i,x,t])
          Nprimeprime_M[k,,x] <- Nprimeprime_M[k,,x] + y[1:n]
          
        }
      }
    } 
    
    # Reproduction
    if (verbose) cat("t =",t,"Apply reproduction function \n")
    # If there is only one age-class, we must force the third dimension
    if (length(dim(phi_F))==2) dim(phi_F)[3] <- 1
    if (length(dim(phi_M))==2) dim(phi_M)[3] <- 1
    
    for (i in 1 : n) {
      if (save.res) {
        if (sum(Nprimeprime_F[,i,],Nprimeprime_M[,i,])==0) { # To save computing time. In the older version it was: if (sum(Nprime_M[,i,])==0 | sum(Nprime_F[,i,])==0)
          L_M[,i] = 0
          L_F[,i] = 0
          next
        } else {
          LL <- repr(Nprimeprime_F[,i,], Nprimeprime_M[,i,], phi_F[,i,,t], phi_M[,i,,t], l, m, z, meiosis_matrix, mat_geno_to_index_mapping)
          L_M[,i] <- LL[,1]
          L_F[,i] <- LL[,2]
          rm(LL)
        }
      } else {
        if (sum(Nprimeprime_F[,i,],Nprimeprime_M[,i,])==0) { # To save computing time. In the older version it was: if (sum(Nprime_M[,i,])==0 | sum(Nprime_F[,i,])==0)
          L_M[,i,t] = 0
          L_F[,i,t] = 0
          next
        } else {
          LL <- repr(Nprimeprime_F[,i,], Nprimeprime_M[,i,], phi_F[,i,,t], phi_M[,i,,t], l, m, z, meiosis_matrix, mat_geno_to_index_mapping)
          L_M[,i,t] <- LL[,1]
          L_F[,i,t] <- LL[,2]
        }
      }
    }
    
    # Propagule dispersal
    if (verbose) cat("t =",t,"Apply propagule dispersal function \n")
    for (i in 1 : n) {
      for (k in 1 : m) {
        if (save.res) {
          y = disp(L_M[k,i],delta.prop[,i,t])
          S_M[k,] <- S_M[k,] + y[1:n]
          y = disp(L_F[k,i],delta.prop[,i,t])
          S_F[k,] <- S_F[k,] + y[1:n] 
        } else {
          y = disp(L_M[k,i,t],delta.prop[,i,t])
          S_M[k,,t] <- S_M[k,,t] + y[1:n]
          y = disp(L_F[k,i,t],delta.prop[,i,t])
          S_F[k,,t] <- S_F[k,,t] + y[1:n]
        }
      }
    }  
    
    # Recruitment
    if (verbose) cat("t =",t,"Apply recruitment function \n")
    for (i in 1 : n) {
      if (save.res) {
        N_M[,i,1] <- recr(S_M[,i],S_F[,i],Nprimeprime_M[,i,],Nprimeprime_F[,i,],m,kappa0[i,t],recr.dd)  # Pass the abundance of other classes too
        N_F[,i,1] <- recr(S_F[,i],S_M[,i],Nprimeprime_M[,i,],Nprimeprime_F[,i,],m,kappa0[i,t],recr.dd)
      } else {
        N_M[,i,1,t+1] <- recr(S_M[,i,t],S_F[,i,t],Nprimeprime_M[,i,],Nprimeprime_F[,i,],m,kappa0[i,t],recr.dd)
        N_F[,i,1,t+1] <- recr(S_F[,i,t],S_M[,i,t],Nprimeprime_M[,i,],Nprimeprime_F[,i,],m,kappa0[i,t],recr.dd)
      }
    }
    
    if (verbose) cat("t =",t,"Calculating N at t+1 \n")
    for (i in 1 : n){
      for (x in 1 : z) {
        if (x == 1) next
        for (k in 1 : m) {
          if (save.res) {
            N_M[k,i,x] <- Nprimeprime_M[k,i,x-1]
            N_F[k,i,x] <- Nprimeprime_F[k,i,x-1]
          } else {
            N_M[k,i,x,t+1] <- Nprimeprime_M[k,i,x-1]
            N_F[k,i,x,t+1] <- Nprimeprime_F[k,i,x-1]
          }
        }
      }
    }
    
    # Save results if save.res=T
    if (save.res){
      if ((t+1) %in% save.res.T) {
        file.name <- paste("N",(t+1),".RData",sep="")
        save(N_M,N_F,file=paste(dir.res.name,file.name,sep="/"))
      }
    }
    
  }
  print(T_max)
  print("...done")
  if (save.res==F) return(list(N_M=N_M,N_F=N_F))
}
