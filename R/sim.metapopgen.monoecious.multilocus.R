sim.metapopgen.monoecious.multilocus <- function(init.par,
                                                 sigma,
                                                 phi_F, phi_M,
                                                 fec.distr_F, fec.distr_M,
                                                 migration="forward",
                                                 delta.prop = NULL, delta.ad = NULL, migr = NULL,
                                                 recr.dd="settlers",
                                                 T_max,
                                                 save.res=F, save.res.T=seq(1:T_max),
                                                 output.var="N",
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
  
  # Check fec.distr_F and fec.distr_M
  if(!fec.distr_F %in% c("fixed","poisson")) stop(paste("Unknown parameter value for fec.distr_F:",fec.distr_F))
  if(!fec.distr_M %in% c("fixed","poisson")) stop(paste("Unknown parameter value for fec.distr_M:",fec.distr_M))
  
  # Check output.var
  for (i.var in 1 : length(output.var)) {
    if(!output.var[i.var] %in% c("N","Nprime","Nprimeprime","G_F","G_M","L","S")) stop(paste0("Unknown variable in output.var:",output.var[i.var]))
  }
  # Read output.var
  output.N <- FALSE
  output.Nprime <- FALSE
  output.Nprimeprime <- FALSE
  output.G_F <- FALSE
  output.G_M <- FALSE
  output.L <- FALSE
  output.S <- FALSE
  if("N" %in% output.var) output.N <- TRUE
  if("Nprime" %in% output.var) output.Nprime <- TRUE
  if("Nprimeprime" %in% output.var) output.Nprimeprime <- TRUE
  if("G_F" %in% output.var) output.G_F <- TRUE # check that it's backward migration
  if("G_M" %in% output.var) output.G_M <- TRUE
  if("L" %in% output.var) output.L <- TRUE
  if("S" %in% output.var) output.S <- TRUE
  
  # .............................Check the existence of dispersal matrices
  if (migration == "forward") {
    if (is.null(delta.prop)) {
      print("Setting propagule dispersal probability (delta.prop)")
      delta.prop <- diag(1,n)
      delta.prop <- array(delta.prop,c(n,n,T_max))
    }
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
    N <- N1
    dimnames(N) <- dimnames(N1)
    rm(N1)    
  } else {
    N           <- array(NA, dim=c(m,n,z,T_max))
    N[,,,1] <- N1
    Nprime      <- array(NA, dim=c(m,n,z,T_max))
    Nprimeprime <- array(0,  dim=c(m,n,z,T_max))
    L           <- array(NA, dim=c(m,n,T_max))
    S           <- array(0, dim=c(m,n,T_max))
    if (output.G_F) {
      G_Ft <- array(NA,c(l,n,T_max))
    }
    if (output.G_M) {
      G_Mt <- array(NA,c(l,n,T_max))
    }
    
    dimnamesN1 <- dimnames(N1)
    dimnamesN1$time <- c(1:T_max)
    dimnames(N) <- dimnamesN1
    rm(N1)
    # Might want to set dimnames for the other objects too
    
  }
  
  
  ##########################################################################
  # Simulate metapopulation genetics
  ##########################################################################
  if (save.res){
    dir.res.name <- paste(getwd(),format(Sys.time(), "%Y-%b-%d-%H.%M.%S"),sep="/")
    dir.create(dir.res.name)

  }
  
  
  print("Running simulation...")
  for (t in 1 : T_max) {
    if (t %% 10 == 0) print(t)
	
    # If save.res, redefine variables
    if (save.res) {
      Nprime        <- array(NA,dim=c(m,n,z))
      Nprimeprime   <- array(0,dim=c(m,n,z))
      L             <- array(NA,dim=c(m,n))
      S             <- array(0,dim=c(m,n))
    }
    
    
    # Survival
    if (verbose) cat("t =",t,"Apply survival function \n")
    # If there is only one age-class, we must force the third dimension. What if only one year?
    if (length(dim(sigma))==2) dim(sigma)[3] <- 1
    for (i in 1 : n) {
      for (x in 1 : z) {
        for (k in 1 : m) {
          if (save.res){
            Nprime[k,i,x] = surv(sigma[k,i,x,t], N[k,i,x])
            if(is.na(Nprime[k,i,x])){ # Needed?
              Nprime[k,i,x] <- 0
              cat("NA adjustment in deme",i,"age",x,"genotype",k,"time",t,"\n")
            }
          } else {
            Nprime[k,i,x,t] = surv(sigma[k,i,x,t], N[k,i,x,t])
            if(is.na(Nprime[k,i,x,t])){ # Needed?
              Nprime[k,i,x,t] <- 0
              cat("NA adjustment in deme",i,"age",x,"genotype",k,"time",t,"\n")
            }
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
            y = disp.ad(Nprime[k,i,x], delta.ad[,i,x,t])
            Nprimeprime[k,,x] <- Nprimeprime[k,,x] + y[1:n] 
          } else {
            y = disp.ad(Nprime[k,i,x,t], delta.ad[,i,x,t])
            Nprimeprime[k,,x,t] <- Nprimeprime[k,,x,t] + y[1:n]     
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
      
      # if (migration == "backward") {
      #   G_F <- array(0,c(l,n))
      #   G_M <- array(0,c(l,n))
      # }
      
      for (i in 1 : n) {
        if (save.res) {
          if (sum(Nprimeprime[,i,])==0) { # To save computing time
            L[,i] = 0
            next
          } else {
            res.repr <- repr(Nprimeprime[,i,], NULL, phi_F[,i,,t], phi_M[,i,,t], l, m, z,
                             meiosis_matrix, mat_geno_to_index_mapping, fec.distr_F, fec.distr_M, migration)
            if (migration == "forward") {
              L[,i] <- res.repr
            } else {
              # G_F[,i] <- res.repr$G_F
              # G_M[,i] <- res.repr$G_M
            }
          }
          
        } else {
          
          if (sum(Nprimeprime[,i,,t])==0) { # To save computing time
            L[,i,t] = 0
            next
          } else {
            res.repr <- repr(Nprimeprime[,i,,t], NULL, phi_F[,i,,t], phi_M[,i,,t], l, m, z,
                             meiosis_matrix, mat_geno_to_index_mapping, fec.distr_F, fec.distr_M, migration)
            if (migration == "forward") {
              L[,i,t] <- res.repr
            } else {
              G_F[,i] <- res.repr$G_F
              G_M[,i] <- res.repr$G_M
              
              if (output.G_F) G_Ft[,,t] <- G_F
              if (output.G_M) G_Mt[,,t] <- G_M
            }
          }
          
        }
        
      }
    }
    
    # Propagule dispersal
    if (migration == "forward") {
      if (verbose) print("Apply propagule dispersal function")
      for (i in 1 : n) {
        for (k in 1 : m) {
          if (save.res) {
            y = disp(L[k,i], delta.prop[,i,t])
            S[k,] <- S[k,] + y[1:n]       
          } else {
            y = disp(L[k,i,t], delta.prop[,i,t])
            S[k,,t] <- S[k,,t] + y[1:n]
          }
        }
      }
    }
    
    
    # Save results if save.res=T
    if (save.res){
      if (output.N) {
        file.name <- paste0("N",t,".RData")
        save(N,file=paste(dir.res.name,file.name,sep="/"))
      }
      if (output.Nprime) {
        file.name <- paste0("Nprime",t,".RData")
        save(Nprime,file=paste(dir.res.name,file.name,sep="/"))
      }
      if (output.Nprimeprime) {
        file.name <- paste0("Nprimeprime",t,".RData")
        save(Nprimeprime,file=paste(dir.res.name,file.name,sep="/"))
      }
      if (output.G_F) {
        file.name <- paste0("G_F",t,".RData")
        save(G_F,file=paste(dir.res.name,file.name,sep="/"))
      }
      if (output.G_M) {
        file.name <- paste0("G_M",t,".RData")
        save(G_M,file=paste(dir.res.name,file.name,sep="/"))
      }
      if (output.L) {
        file.name <- paste0("L",t,".RData")
        save(L,file=paste(dir.res.name,file.name,sep="/"))
      }
      if (output.S) {
        file.name <- paste0("S",t,".RData")
        save(S,file=paste(dir.res.name,file.name,sep="/"))
      }
    }

    # Recruitment
    if (t == T_max) break # otherwise attempts to write on T_max + 1
    
    # Recruitment with backward migration
    if (migration == "backward") {
      # This gives the number of gametes of each type (column) produced per deme (row) 
      gamete_pool <- t(Nprimeprime[,,1,t]) %*% t(meiosis_matrix)# 1 bcs only one age class
      for (j in 1 : n) {
        res.recr <- recr.backward.migration(migr, l, m, n, z, j.local.deme=j, kappa0[j,t+1], # without G_F, G_M, 
                                            sexuality="monoecious", mat_geno_to_index_mapping, t=t+1,
                                            gamete_pool = gamete_pool)
        if (save.res) {
          N[,j,1] <- res.recr # 1 bcs only 1 age class is supported now
        } else {
          N[,j,1,t+1] <- res.recr # idem
        }

      }
      next
    }
    
    if (verbose) print("Apply recruitment function")
    for (i in 1 : n) {
      if (save.res) {
        # This overwrites N
        N[,i,] <- recr(N = array(Nprimeprime[,i,], dim=c(m,z)),
                       S = array(S[,i], dim=c(m,1)),
                       m = m,
                       z = z,
                       kappa0 = kappa0[i,t+1], # Pay attention here: using carrying capacity of next generation
                       recr.dd = recr.dd,
                       sexuality="monoecious")
      } else {
        N[,i,,t+1] <- recr(N = array(Nprimeprime[,i,,t], dim=c(m,z)),
                           S = array(S[,i,t], dim=c(m,1)),
                           m = m,
                           z = z,
                           kappa0 = kappa0[i,t+1],
                           recr.dd = recr.dd,
                           sexuality = "monoecious")
      }
    }

  }
  print("...done")
  
  if (save.res==F) {
    output.res <- list()
    if (output.N) {
      output.res$N <- N
    }
    if (output.Nprime) {
      output.res$Nprime <- Nprime
    }
    if (output.Nprimeprime) {
      output.res$Nprimeprime <- Nprimeprime
    }
    if (output.G_F) {
      output.res$G_F <- G_Ft
    }
    if (output.G_M) {
      output.res$G_M <- G_Mt
    }
    if (output.L) {
      output.res$L <- L
    }
    if (output.S) {
      output.res$S <- S
    }
    return(output.res)
  }
}
