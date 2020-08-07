sim.metapopgen.monoecious <- function(input.type,demographic.data,
                                      N1,
                                      sigma,
                                      phi_F, phi_M,
                                      fec.distr_F = "poisson", fec.distr_M = "poisson",
                                      mu,
                                      migration="forward", 
                                      delta, migr=NULL,
                                      recr.dd="settlers",
                                      kappa0,
                                      T_max,save.res=F,save.res.T=seq(1,T_max),
                                      verbose=F) {

  ##########################################################################
  
  # Initial definitions
  
  ##########################################################################
  
  
  if (input.type=="data.frame") {
    
    print("Input type = data.frame")
    
    a <- metapopgen.input.convert.monoecious(demographic.data)
    N1    <- a[[1]]
    sigma <- a[[2]]
    phi_M <- a[[3]]
    phi_F <- a[[4]]
    rm(a)
    
    
  } else {
    if (input.type == "array") {
      
      print("Input type = array")
      
    } else {
      stop("Unknown value for argument input.type. It must be either data.frame, array or txt")
      
    }
  }
  

  
  # Define basic variables
  
  m <- dim(N1)[1]                 # Number of genotypes
  l <- (sqrt(1+8*m)-1)/2          # Nuber of alleles
  n <- dim(N1)[2]                 # Number of demes
  z <- dim(N1)[3]                 # Number of age-classes
  
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
  
  # Survival
  if (is.na(dim(sigma)[4])) {
      sigma <- array(rep(sigma,T_max),c(m,n,z,T_max))
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
    N <- N1
    rm(N1)    
  } else {
    N       <- array(NA,dim=c(m,n,z,T_max))
    N[,,,1] <- N1
    rm(N1)
    L       <- array(NA,dim=c(m,n,T_max))
    S       <- array(0,dim=c(m,n,T_max))
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
  for (t in 1 : (T_max)) {
    
    if (t %% 10 == 0) print(t)

    # At each time-step, redefine variable Nprime
    # If save.res, redefine also larval and settlers numbers
    if (save.res) {
      Nprime  <- array(NA,dim=c(m,n,z))
      L       <- array(NA,dim=c(m,n))
      S       <- array(0,dim=c(m,n))
    } else {
      Nprime  <- array(NA,dim=c(m,n,z))
    }
    
    
    ### Survival
    
    # If there is only one age-class, we must force the third dimension. What if only one year?
    if (length(dim(sigma))==2) dim(sigma)[3] <- 1
    
    if (verbose) print("Apply survival function")
    for (i in 1 : n) {
      for (x in 1 : z) {
        for (k in 1 : m) {
          
          if (save.res){
            Nprime[k,i,x] = surv(sigma[k,i,x,t],N[k,i,x])
          } else {
            Nprime[k,i,x] = surv(sigma[k,i,x,t],N[k,i,x,t])
          }
        }
      }
    }
    
    
    ## Reproduction
    
    if (verbose) print("Apply reproduction function")

    # If there is only one age-class, we must force the third dimension
    if (length(dim(phi_F))==2) dim(phi_F)[3] <- 1
    if (length(dim(phi_M))==2) dim(phi_M)[3] <- 1
    
    if (migration == "backward") {
      G_F <- array(0,c(l,n))
      G_M <- array(0,c(l,n))
    }
    
    for (i in 1 : n) {
      
      if (save.res) {
        
        if (sum(Nprime[,i,])==0) { # To save computing time
          L[,i] = 0
          next
        } else {
          res.repr <- repr.monoecious.onelocus(Nprime, phi_F[,,,t], phi_M[,,,t], mu, i, l, m, n, z,
                                 fec.distr_F = fec.distr_F, fec.distr_M = fec.distr_M, migration=migration, verbose=verbose)
          if (migration == "forward") {
            L[,i] <- res.repr
          } else {
            G_F[,i] <- res.repr$G_F
            G_M[,i] <- res.repr$G_M
          }
        }
        
      } else {
        
        if (sum(Nprime[,i,])==0) { # To save computing time
          L[,i,t] = 0
          next
        } else {
          res.repr <- repr.monoecious.onelocus(Nprime,phi_F[,,,t], phi_M [,,,t],mu, i, l, m, n, z,
                                   fec.distr_F=fec.distr_F, fec.distr_M=fec.distr_M, migration=migration, verbose=verbose)
          if (migration == "forward") {
            L[,i,t] <- res.repr
          } else {
            G_F[,i] <- res.repr$G_F
            G_M[,i] <- res.repr$G_M

          }
        }
        
      }
    }
    
    
    ## Propagule dispersal
    if (migration == "forward") {
      if (verbose) print("Apply dispersal function")
      for (i in 1 : n) {
        for (k in 1 : m) {
          if (save.res) {
            y = disp(L[k,i],delta[,i,t])
            S[k,] <- S[k,] + y[1:n]       
          } else {
            y = disp(L[k,i,t],delta[,i,t])
            S[k,,t] <- S[k,,t] + y[1:n]
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
    
    # Recruitment
    if (t == T_max) break # otherwise attempts to write on T_max + 1
    
    # Recruitment with backward migration
    if (migration == "backward") {
      for (j in 1 : n) {
        res.recr <- recr.backward.migration.onelocus(G_F, G_M, migr, l, m, n, z, j.local.deme=j, kappa0[j,t+1],
                                            sexuality="monoecious", t=t+1)
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
        N[,i,] <- recr(N = array(Nprime[,i,], dim=c(m,z)),
                       S = array(S[,i], dim=c(m,1)),
                       m = m,
                       z = z,
                       kappa0 = kappa0[i,t+1], # Pay attention here: using carrying capacity of next generation
                       recr.dd = recr.dd,
                       sexuality="monoecious")
      } else {
        N[,i,,t+1] <- recr(N = array(Nprime[,i,], dim=c(m,z)),
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
  if (save.res==F) return(N)
}
