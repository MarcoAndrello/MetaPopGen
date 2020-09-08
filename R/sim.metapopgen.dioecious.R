sim.metapopgen.dioecious <- function(input.type, demographic.data,
                                     N1_M, N1_F,
                                     sigma_M, sigma_F,
                                     phi_M, phi_F,
                                     fec.distr_F = "poisson", fec.distr_M = "poisson",
                                     mu,
                                     migration = "forward", migr=NULL,
                                     delta,
                                     recr.dd = "settlers",
                                     kappa0,
                                     T_max,
                                     save.res = F, save.res.T = seq(1,T_max),
                                     verbose = F) {
  
  
  ##########################################################################
  
  # Initial definitions
  
  ##########################################################################
  
  
  if (input.type=="data.frame") {
    
    print("Input type = data.frame")
    
    a <- metapopgen.input.convert.dioecious(demographic.data)
    N1_M    <- a[[1]]
    N1_F    <- a[[2]]
    sigma_M <- a[[3]]
    sigma_F <- a[[4]]
    phi_M <- a[[5]]
    phi_F <- a[[6]]
    rm(a)
    
    
  } else {
    if (input.type == "array") {
      
      print("Input type = array")
      
    } else {
      stop("Unknown value for argument input.type. It must be either data.frame, array or txt")
      
    }
  }
  

  # Define basic variables
  
  m <- dim(N1_M)[1]                 # Number of genotypes
  l <- (sqrt(1+8*m)-1)/2            # Nuber of alleles
  n <- dim(N1_M)[2]                 # Number of demes
  z <- dim(N1_M)[3]                 # Number of age-classes
  
  # Check fec.distr_F and fec.distr_M
  if(!fec.distr_F %in% c("fixed","poisson")) stop(paste("Unknown parameter value for fec.distr_F:",fec.distr_F))
  if(!fec.distr_M %in% c("fixed","poisson")) stop(paste("Unknown parameter value for fec.distr_M:",fec.distr_M))
  
  # If only one age-class and recr.dd=="adults", gives an error
  if (z == 1 & recr.dd == "adults") {
    stop("Detected only one age class (z=1) and recruitment probability dependent on adult density (recr.dd == 'adults'). This combination is not supported. Use recr.dd == 'settlers' instead.")
  }
  
  
  ##########################################################################
  
  # Check if input data are time-dependent or not
  
  ##########################################################################
  
  # Male Survival
  if (is.na(dim(sigma_M)[4])) {
    sigma_M <- array(rep(sigma_M,T_max),c(m,n,z,T_max))
  }
  
  # Female Survival
  if (is.na(dim(sigma_F)[4])) {
    sigma_F <- array(rep(sigma_F,T_max),c(m,n,z,T_max))
  }
  
  # Female fecundity
  if (is.na(dim(phi_F)[4])) {
    phi_F <- array(rep(phi_F,T_max),c(m,n,z,T_max))
  }
  
  # Male fecundity
  if (is.na(dim(phi_M)[4])) {
    phi_M <- array(rep(phi_M,T_max),c(m,n,z,T_max))
  }
  
  # Dispersal
  if (is.na(dim(delta)[3])) {
    delta <- array(rep(delta,T_max),c(n,n,T_max))
  }
  
  # Carrying capacity
  if (is.vector(kappa0)) {
    kappa0 <- array(rep(kappa0,T_max),c(n,T_max))
  }
  
  
  ##########################################################################
  
  # Initialize state variables
  
  ##########################################################################
  
  print("Initializing variables...")
  if (save.res){
    N_M <- N1_M
    N_F <- N1_F
    rm(N1_M,N1_F)    
  } else {
    N_M       <- array(NA,dim=c(m,n,z,T_max))
    N_M[,,,1] <- N1_M
    N_F       <- array(NA,dim=c(m,n,z,T_max))
    N_F[,,,1] <- N1_F
    rm(N1_M,N1_F)
    L_M       <- array(NA,dim=c(m,n,T_max))
    L_F       <- array(NA,dim=c(m,n,T_max))
    S_M       <- array(0,dim=c(m,n,T_max))
    S_F       <- array(0,dim=c(m,n,T_max))
  }
  
  
  ##########################################################################
  
  # Define functions
  
  ##########################################################################
  
  

  
  
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

    # At each time-step, redefine variable Nprime
    # If save.res, redefine also larval and settlers numbers
    if (save.res) {
      Nprime_M  <- array(NA,dim=c(m,n,z))
      Nprime_F  <- array(NA,dim=c(m,n,z))
      L_M       <- array(NA,dim=c(m,n))
      L_F       <- array(NA,dim=c(m,n))
      S_M       <- array(0,dim=c(m,n))
      S_F       <- array(0,dim=c(m,n))
    } else {
      Nprime_M  <- array(NA,dim=c(m,n,z))
      Nprime_F  <- array(NA,dim=c(m,n,z))
    }
    
    
    ### Survival
    
    # If there is only one age-class, we must force the third dimension. What if only one year?
    if (length(dim(sigma_M))==2) dim(sigma_M)[3] <- 1
    
    # If there is only one age-class, we must force the third dimension. What if only one year?
    if (length(dim(sigma_F))==2) dim(sigma_F)[3] <- 1
    
    if (verbose) print("Apply survival function")
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
        }
      }
    }  
    
    if (verbose) print("Apply reproduction function")

    
    # If there is only one age-class, we must force the third dimension
    if (length(dim(phi_F))==2) dim(phi_F)[3] <- 1
    if (length(dim(phi_M))==2) dim(phi_M)[3] <- 1
    
    for (i in 1 : n) {
      
      if (save.res) {
        if (sum(Nprime_M[,i,])==0 | sum(Nprime_F[,i,])==0) { # To save computing time
          L_M[,i] = 0
          L_F[,i] = 0
          next
        } else {
          LL <- repr.dioecious.onelocus(Nprime_M, Nprime_F, phi_F[,,,t], phi_M[,,,t],
                                        mu, i, l, m, n, z,
                                        fec.distr_F, fec.distr_M, migration, verbose)
          L_M[,i] <- LL[,1]
          L_F[,i] <- LL[,2]
          rm(LL)
        }
        
      } else {
        if (sum(Nprime_M[,i,])==0 | sum(Nprime_F[,i,])==0) { # To save computing time
          L_M[,i,t] = 0
          L_F[,i,t] = 0
          next
        } else {
          LL <- repr.dioecious.onelocus(Nprime_M, Nprime_F, phi_F[,,,t], phi_M[,,,t],
                                        mu, i, l, m, n, z,
                                        fec.distr_F, fec.distr_M, migration, verbose)
          L_M[,i,t] <- LL[,1]
          L_F[,i,t] <- LL[,2]
        }
        
      }
    }
    
    
    if (verbose) print("Apply dispersal function")
    for (i in 1 : n) {
      for (k in 1 : m) {
        if (save.res) {
          y = disp(L_M[k,i],delta[,i,t])
          S_M[k,] <- S_M[k,] + y[1:n]
          y = disp(L_F[k,i],delta[,i,t])
          S_F[k,] <- S_F[k,] + y[1:n] 
        } else {
          y = disp(L_M[k,i,t],delta[,i,t])
          S_M[k,,t] <- S_M[k,,t] + y[1:n]
          y = disp(L_F[k,i,t],delta[,i,t])
          S_F[k,,t] <- S_F[k,,t] + y[1:n]
        }
      }
    }  
    
    
    if (verbose) print("Apply recruitment function")
    for (i in 1 : n) {
      if (save.res) {
        Naged <- recr(N_F = array(Nprime_F[,i,], dim=c(m,z)),
                      N_M = array(Nprime_M[,i,], dim=c(m,z)),
                      S_F = array(S_F[,i], dim=c(m,1)),
                      S_M = array(S_M[,i], dim=c(m,1)),
                      m = m,
                      z = z,
                      kappa0 = kappa0[i,t],
                      recr.dd = recr.dd,
                      sexuality = "dioecious")
        N_F[,i,] <- Naged[,,1]
        N_M[,i,] <- Naged[,,2]
      } else {
        Naged <- recr(N_F = array(Nprime_F[,i,], dim=c(m,z)),
                      N_M = array(Nprime_M[,i,], dim=c(m,z)),
                      S_F = array(S_F[,i,t], dim=c(m,1)),
                      S_M = array(S_M[,i,t], dim=c(m,1)),
                      m = m,
                      z = z,
                      kappa0 = kappa0[i,t],
                      recr.dd = recr.dd,
                      sexuality = "dioecious")
        N_F[,i,,t+1] <- Naged[,,1] 
        N_M[,i,,t+1] <- Naged[,,2]
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
