sim.metapopgen.dioecious.multilocus <- function(init.par,
                                                sigma_F, sigma_M,
                                                phi_F, phi_M,
                                                fec.distr_F = "poisson", fec.distr_M = "poisson",
                                                migration = "forward", migr,
                                                delta.prop = NULL, delta.ad = NULL,
                                                recr.dd="settlers",
                                                T_max,
                                                save.res = F, save.res.T = seq(1:T_max),
                                                output.var = "N",
                                                verbose = F) {
  
  # Reading basic variables
  m                         <- init.par$m                          # Number of genotypes
  meiosis_matrix            <- init.par$meiosis_matrix             # Meiosis matrix
  mat_geno_to_index_mapping <- init.par$mat_geno_to_index_mapping  # Mapping of combinations of gametotypes to genotypes
  if (is.null(mat_geno_to_index_mapping)) mat_geno_to_index_mapping <- init.par$index_matr # This is for gamete-based representation
  l                         <- dim(meiosis_matrix)[1]              # Number of multilocus gametotypes
  n                         <- init.par$n                          # Number of demes
  z                         <- init.par$z                          # Number of age-classes
  kappa0                    <- init.par$kappa0                     # Carrying capacity
  N1_F                      <- init.par$N1_F                       # Initial composition
  N1_M                      <- init.par$N1_M                       # Initial composition
  
  # If only one age-class and recr.dd=="adults", gives an error
  if (z == 1 & recr.dd == "adults") {
    stop("Detected only one age class (z=1) and recruitment probability dependent on adult density (recr.dd == 'adults'). This combination is not supported. Use recr.dd == 'settlers' instead.")
  }
  
  # Check fec.distr_F and fec.distr_M
  if(!fec.distr_F %in% c("fixed","poisson")) stop(paste("Unknown parameter value for fec.distr_F:",fec.distr_F))
  if(!fec.distr_M %in% c("fixed","poisson")) stop(paste("Unknown parameter value for fec.distr_M:",fec.distr_M))

  # Check output.var
  for (i.var in 1 : length(output.var)) {
    if(!output.var[i.var] %in% c("N","Nprime","Nprimeprime","L","S")) stop(paste0("Unknown variable in output.var:",output.var[i.var]))
  }
  # Read output.var
  output.N <- FALSE
  output.Nprime <- FALSE
  output.Nprimeprime <- FALSE
  output.L <- FALSE
  output.S <- FALSE
  if("N" %in% output.var) output.N <- TRUE
  if("Nprime" %in% output.var) output.Nprime <- TRUE
  if("Nprimeprime" %in% output.var) output.Nprimeprime <- TRUE
  if("L" %in% output.var) output.L <- TRUE
  if("S" %in% output.var) output.S <- TRUE
  
  if (migration == "backward" & output.L == TRUE) stop("Variable \"L\" is not calculated with backward migration. Correct argument \"output.var\"")
  if (migration == "backward" & output.S == TRUE) stop("Variable \"S\" is not calculated with backward migration. Correct argument \"output.var\"")
  
  
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
  
  # Adult dispersal
  if (is.na(dim(delta.ad)[3])) {
    print("Augmenting delta.ad for age dimension")
    delta.ad <- array(rep(delta.ad,z),c(n,n,z))
  }
  if (is.na(dim(delta.ad)[4])) {
    print("Augmenting delta.ad for time dimension")
    delta.ad <- array(rep(delta.ad,T_max),c(n,n,z,T_max))
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
  
  # Propagule dispersal
  if (migration == "forward") {
    if (is.na(dim(delta.prop)[3])) {
      print("Augmenting delta.prop for time dimension")
      delta.prop <- array(rep(delta.prop,T_max),c(n,n,T_max))
    }
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
    
    dimnamesN1      <- dimnames(N1_F)
    dimnamesN1noage <- dimnames(N1_F)[c(1,2)]
    
    dimnames(N_F) <- dimnamesN1
    dimnames(N_M) <- dimnamesN1

    rm(N1_M, N1_F) 
    
  } else {
    N_F                <- array(NA, dim=c(m,n,z,T_max))
    N_F[,,,1] <- N1_F
    N_M                <- array(NA, dim=c(m,n,z,T_max))
    N_M[,,,1] <- N1_M
    Nprime_F           <- array(NA, dim=c(m,n,z,T_max))
    Nprime_M           <- array(NA, dim=c(m,n,z,T_max))
    Nprimeprime_F      <- array(0,  dim=c(m,n,z,T_max))
    Nprimeprime_M      <- array(0,  dim=c(m,n,z,T_max))
    L_F                <- array(NA, dim=c(m,n,T_max))
    L_M                <- array(NA, dim=c(m,n,T_max))
    S_F                <- array(0,  dim=c(m,n,T_max))
    S_M                <- array(0,  dim=c(m,n,T_max))
    
    dimnamesN1 <- dimnames(N1_M)
    dimnamesN1$time <- c(1:T_max)

    dimnames(N_F)           <- dimnamesN1
    dimnames(N_M)           <- dimnamesN1
    dimnames(Nprime_F)      <- dimnamesN1
    dimnames(Nprime_M)      <- dimnamesN1
    dimnames(Nprimeprime_F) <- dimnamesN1
    dimnames(Nprimeprime_M) <- dimnamesN1
    dimnamesN1$age <- NULL
    dimnames(L_F)           <- dimnamesN1
    dimnames(L_M)           <- dimnamesN1
    dimnames(S_F)           <- dimnamesN1
    dimnames(S_M)           <- dimnamesN1
    
    rm(N1_M, N1_F, dimnamesN1)

  }
  
  # Create output folder
  if (save.res){
    dir.res.name <- paste(getwd(),format(Sys.time(), "%Y-%b-%d-%H.%M.%S"),sep="/")
    dir.create(dir.res.name)
  }
  
  ##########################################################################
  # Simulate metapopulation genetics
  ##########################################################################

  
  print("Running simulation...")
  for (t in 1 : T_max) {
    if (t %% 10 == 0) print(t)
    
    # If save.res, redefine variables
    if (save.res) {
      Nprime_F  <- array(NA,dim=c(m,n,z))
      Nprime_M  <- array(NA,dim=c(m,n,z))
      Nprimeprime_F  <- array(0,dim=c(m,n,z))
      Nprimeprime_M  <- array(0,dim=c(m,n,z))
      L_F       <- array(NA,dim=c(m,n))
      L_M       <- array(NA,dim=c(m,n))
      S_F       <- array(0,dim=c(m,n))
      S_M       <- array(0,dim=c(m,n))
      
      dimnames(Nprime_F) <- dimnamesN1
      dimnames(Nprime_M) <- dimnamesN1
      dimnames(Nprimeprime_F) <- dimnamesN1
      dimnames(Nprimeprime_M) <- dimnamesN1
      dimnames(L_F) <- dimnamesN1noage
      dimnames(L_M) <- dimnamesN1noage
      dimnames(S_F) <- dimnamesN1noage
      dimnames(S_M) <- dimnamesN1noage
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
            Nprime_M[k,i,x] = surv(sigma_M[k,i,x,t], N_M[k,i,x])
            Nprime_F[k,i,x] = surv(sigma_F[k,i,x,t], N_F[k,i,x])
          } else {
            Nprime_M[k,i,x,t] = surv(sigma_M[k,i,x,t], N_M[k,i,x,t])
            Nprime_F[k,i,x,t] = surv(sigma_F[k,i,x,t], N_F[k,i,x,t])
          }
        }
      }
    }  
    
    # Adult dispersal
    if (verbose) cat("t =",t,"Apply adult dispersal function \n") 
    for (i in 1 : n) {
      for (x in 1 : z) {
        for (k in 1 : m) {
          if (save.res){
            y = disp.ad(Nprime_F[k,i,x], delta.ad[,i,x,t])
            Nprimeprime_F[k,,x] <- Nprimeprime_F[k,,x] + y[1:n]
            y = disp.ad(Nprime_M[k,i,x], delta.ad[,i,x,t])
            Nprimeprime_M[k,,x] <- Nprimeprime_M[k,,x] + y[1:n]
          } else {
            y = disp.ad(Nprime_F[k,i,x,t], delta.ad[,i,x,t])
            Nprimeprime_F[k,,x,t] <- Nprimeprime_F[k,,x,t] + y[1:n]
            y = disp.ad(Nprime_M[k,i,x,t], delta.ad[,i,x,t])
            Nprimeprime_M[k,,x,t] <- Nprimeprime_M[k,,x,t] + y[1:n]
          }
        }
      }
    } 
    
    # Reproduction
    if (migration == "forward") { 
      if (verbose) cat("t =",t,"Apply reproduction function \n")
      # If there is only one age-class, we must force the third dimension
      if (length(dim(phi_F))==2) dim(phi_F)[3] <- 1
      if (length(dim(phi_M))==2) dim(phi_M)[3] <- 1
      
      for (i in 1 : n) {
        if (save.res) {
          if (sum(Nprimeprime_F[,i,], Nprimeprime_M[,i,])==0) { # To save computing time. In the older version it was: if (sum(Nprime_M[,i,])==0 | sum(Nprime_F[,i,])==0)
            L_M[,i] = 0
            L_F[,i] = 0
            next
          } else {
            LL <- repr(Nprimeprime_F[,i,], Nprimeprime_M[,i,], phi_F[,i,,t], phi_M[,i,,t], l, m, z,
                       meiosis_matrix, mat_geno_to_index_mapping, fec.distr_F, fec.distr_M, migration)
            L_M[,i] <- LL[,1]
            L_F[,i] <- LL[,2]
            rm(LL)
          }
        } else {
          if (sum(Nprimeprime_F[,i,,t], Nprimeprime_M[,i,,t])==0) { # To save computing time. In the older version it was: if (sum(Nprime_M[,i,])==0 | sum(Nprime_F[,i,])==0)
            L_M[,i,t] = 0
            L_F[,i,t] = 0
            next
          } else {
            LL <- repr(Nprimeprime_F[,i,,t], Nprimeprime_M[,i,,t], phi_F[,i,,t], phi_M[,i,,t], l, m, z,
                       meiosis_matrix, mat_geno_to_index_mapping, fec.distr_F, fec.distr_M, migration)
            L_M[,i,t] <- LL[,1]
            L_F[,i,t] <- LL[,2]
          }
        }
      }
    }
    
    # Propagule dispersal
    if (migration == "forward") {
      if (verbose) cat("t =",t,"Apply propagule dispersal function \n")
      for (i in 1 : n) {
        for (k in 1 : m) {
          if (save.res) {
            y = disp(L_M[k,i],delta.prop[,i,t])
            S_M[k,] <- S_M[k,] + y[1:n]
            y = disp(L_F[k,i],delta.prop[,i,t])
            S_F[k,] <- S_F[k,] + y[1:n] 
          } else {
            y = disp(L_M[k,i,t], delta.prop[,i,t])
            S_M[k,,t] <- S_M[k,,t] + y[1:n]
            y = disp(L_F[k,i,t], delta.prop[,i,t])
            S_F[k,,t] <- S_F[k,,t] + y[1:n]
          }
        }
      }
    }
    
    # Save results if save.res=T
    if (save.res){
      if (output.N) {
        file.name <- paste0("N",t,".RData")
        save(N_F, N_M, file=paste(dir.res.name,file.name,sep="/"))
      }
      if (output.Nprime) {
        file.name <- paste0("Nprime",t,".RData")
        save(Nprime_F, Nprime_M, file=paste(dir.res.name,file.name,sep="/"))
      }
      if (output.Nprimeprime) {
        file.name <- paste0("Nprimeprime",t,".RData")
        save(Nprimeprime_F, Nprimeprime_M, file=paste(dir.res.name,file.name,sep="/"))
      }
      if (output.L) {
        file.name <- paste0("L",t,".RData")
        save(L_F, L_M, file=paste(dir.res.name,file.name,sep="/"))
      }
      if (output.S) {
        file.name <- paste0("S",t,".RData")
        save(S_F, S_M, file=paste(dir.res.name,file.name,sep="/"))
      }
    }
    
    
    # Recruitment
    if (t == T_max) break # otherwise attempts to write on T_max + 1
    
    # Recruitment with backward migration
    if (migration == "backward") {
      # Not implemented yet. See the corresponding section of sim.metapopgen.monoecious.multilocus
      next
    }
    
    if (verbose) cat("t =",t,"Apply recruitment function \n")
    for (i in 1 : n) {
      if (save.res) {
        Naged <- recr(N_F = array(Nprimeprime_F[,i,], dim=c(m,z)),
                      N_M = array(Nprimeprime_M[,i,], dim=c(m,z)),
                      S_F = array(S_F[,i], dim=c(m,1)),
                      S_M = array(S_M[,i], dim=c(m,1)),
                      m = m,
                      z = z,
                      kappa0 = kappa0[i,t+1],
                      recr.dd = recr.dd,
                      sexuality = "dioecious")
        N_F[,i,] <- Naged[,,1]
        N_M[,i,] <- Naged[,,2]
        
      } else {
        Naged <- recr(N_F = array(Nprimeprime_F[,i,,t], dim=c(m,z)),
                      N_M = array(Nprimeprime_M[,i,,t], dim=c(m,z)),
                      S_F = array(S_F[,i,t], dim=c(m,1)),
                      S_M = array(S_M[,i,t], dim=c(m,1)),
                      m = m,
                      z = z,
                      kappa0 = kappa0[i,t+1],
                      recr.dd = recr.dd,
                      sexuality = "dioecious")
        N_F[,i,,t+1] <- Naged[,,1] 
        N_M[,i,,t+1] <- Naged[,,2]
      }
    }
  }
  print("...done")
  
  if (save.res==F) {
    output.res <- list()
    if (output.N) {
      output.res$N_F <- N_F
      output.res$N_M <- N_M
    }
    if (output.Nprime) {
      output.res$Nprime_F <- Nprime_F
      output.res$Nprime_M <- Nprime_M
    }
    if (output.Nprimeprime) {
      output.res$Nprimeprime_F <- Nprimeprime_F
      output.res$Nprimeprime_M <- Nprimeprime_M
    }
    if (output.L) {
      output.res$L_F <- L_F
      output.res$L_M <- L_M
    }
    if (output.S) {
      output.res$S_F <- S_F
      output.res$S_M <- S_M
    }
    return(output.res)
  }
  
}
