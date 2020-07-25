sim.metapopgen.monoecious.multilocus <- function(init.par,
                                                 sigma,
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
  N1                        <- init.par$N1                         # Initial composition
  
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
  
  # Survival
  if (is.na(dim(sigma)[4])) {
      print("Augmenting sigma for time dimension")
    sigma <- array(rep(sigma,T_max),c(m,n,z,T_max))
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
    N <- N1
    dimnames(N) <- dimnames(N1)
    rm(N1)    
  } else {
    N       <- array(NA,dim=c(m,n,z,T_max))
    N[,,,1] <- N1
    dimnamesN1 <- dimnames(N1)
    rm(N1)
    dimnamesN1$time <- c(1:T_max)
    dimnames(N) <- dimnamesN1
    L       <- array(NA,dim=c(m,n,T_max))
    S       <- array(0,dim=c(m,n,T_max))
  }
  
  
  ##########################################################################
  # Define functions
  ##########################################################################

  # Recruitment
  # S        Number of settlers of all genotypes. Dimension: m
  # N        Number of adults of all genotyps and age-classes. Dimensions: m*z
  # m        Number of genotypes
  # kappa0   Carrying capacity. Scalar.
  # recr.dd  String
  # S[,i],Nprime[,i,],m,kappa0[i,t],recr.dd
  #S <- S[,i]
  #N <- Nprime[,i,]
  #kappa0 <- kappa0[i,t]
  
  recr <-
    function(S,N,m,kappa0,recr.dd) {
      switch(recr.dd,
             # Dependence on settler density
             settlers = {
               Ntot <- sum(S)
               sigma0 <- settler.survival(Ntot,kappa0)
             },
             # Dependence on adult density
             adults = {
               if (z==1) Ntot <- 0 else Ntot <- sum(N[,1:(z-1)]) # Does not count z because they will "shift out" with the aging function
               Stot <- sum(S)
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
      # Use recruitment probability to calculate the number of recruits
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
      save(N,file=paste(dir.res.name,file.name,sep="/"))
    }
  }
  
  
  print("Running simulation...")
  for (t in 1 : (T_max-1)) {
    if (t %% 10 == 0) print(t)
	
    # At each time-step, redefine variable Nprime
    # If save.res, redefine also larval and settlers numbers
    if (save.res) {
      Nprime        <- array(NA,dim=c(m,n,z))
      Nprimeprime   <- array(0,dim=c(m,n,z))
      L             <- array(NA,dim=c(m,n))
      S             <- array(0,dim=c(m,n))
    } else {
      Nprime       <- array(NA,dim=c(m,n,z))
      Nprimeprime  <- array(0,dim=c(m,n,z))
      L            <- array(NA,dim=c(m,n,T_max))
      S            <- array(0,dim=c(m,n,T_max))
    }
    
    
    ### Survival
    if (verbose) cat("t =",t,"Apply survival function \n")
    # If there is only one age-class, we must force the third dimension. What if only one year?
    if (length(dim(sigma))==2) dim(sigma)[3] <- 1
    for (i in 1 : n) {
      for (x in 1 : z) {
        for (k in 1 : m) {
          if (save.res){
            Nprime[k,i,x] = surv(sigma[k,i,x,t],N[k,i,x])
          } else {
            Nprime[k,i,x] = surv(sigma[k,i,x,t],N[k,i,x,t])
          }
		  if(is.na(Nprime[k,i,x])){
              Nprime[k,i,x] <- 0
              cat("NA adjustment in deme",i,"age",x,"genotype",k,"\n")
            }
        }
      }
    }
    
    # Adult dispersal
    if (verbose) cat("t =",t,"Apply adult dispersal function \n")
    for (i in 1 : n) {
        for (x in 1 : z) {
            for (k in 1 : m) {
                y = disp.ad(Nprime[k,i,x],delta.ad[,i,x,t])
                Nprimeprime[k,,x] <- Nprimeprime[k,,x] + y[1:n]       
                
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
        if (sum(Nprimeprime[,i,])==0) { # To save computing time
          L[,i] = 0
          next
        } else {
          L[,i] <- repr(Nprimeprime[,i,], NULL, phi_F[,i,,t], phi_M[,i,,t], l, m, z, meiosis_matrix, mat_geno_to_index_mapping)
        }
        
      } else {
        
        if (sum(Nprimeprime[,i,])==0) { # To save computing time
          L[,i,t] = 0
          next
        } else {
          L[,i,t] <- repr(Nprimeprime[,i,], NULL, phi_F[,i,,t], phi_M[,i,,t], l, m, z, meiosis_matrix, mat_geno_to_index_mapping)
        }
        
      }
    }
    
    # Propagule dispersal
    if (verbose) print("Apply propagule dispersal function")
    for (i in 1 : n) {
      for (k in 1 : m) {
        if (save.res) {
          y = disp(L[k,i],delta.prop[,i,t])
          S[k,] <- S[k,] + y[1:n]       
        } else {
          y = disp(L[k,i,t],delta.prop[,i,t])
          S[k,,t] <- S[k,,t] + y[1:n]
        }
      }
    }  
    
    # Recruitment
    if (verbose) print("Apply recruitment function")
    for (i in 1 : n) {
      if (save.res) {
        N[,i,1] <- recr(S[,i],Nprimeprime[,i,],m,kappa0[i,t],recr.dd)
      } else {
        N[,i,1,t+1] <- recr(S[,i,t],Nprimeprime[,i,],m,kappa0[i,t],recr.dd)
      }
    }
    
    
    if (verbose) print("Calculates N at t+1")
    for (i in 1 : n){
      for (x in 1 : z) {
        if (x == 1) next
        for (k in 1 : m) {
          if (save.res) {
            N[k,i,x] <- Nprimeprime[k,i,x-1]
          } else {
            N[k,i,x,t+1] <- Nprimeprime[k,i,x-1]
          }
        }
      }
    }
    
    # Save results if save.res=T
    if (save.res){
      if ((t+1) %in% save.res.T) {
        file.name <- paste("N",(t+1),".RData",sep="")
        save(N,file=paste(dir.res.name,file.name,sep="/"))
      }
    }
    
    
  }
  print(T_max)
  print("...done")
  if (save.res==F)
    return(N)
}
