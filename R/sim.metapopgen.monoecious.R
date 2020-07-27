sim.metapopgen.monoecious <- function(input.type,demographic.data,
                                      N1,sigma,phi_F,phi_M,
                                      mu,delta,recr.dd="settlers",
                                      kappa0,T_max,save.res=F,save.res.T=seq(1,T_max),verbose=F) {

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
  
  # Define functions
  
  ##########################################################################
  
  # Reproduction
  repr.onelocus <-
    function(Nprime,phi_F,phi_M,mu,i,l,m) {
	
		# Adjust dimensions if there is only one age-class and/or one deme
		#if (length(dim(phi_F))==2) dim(phi_F)[3]<-1
		#if (length(dim(phi_M))==2) dim(phi_M)[3]<-1
      dim(phi_F) <- c(m,n,z)
      dim(phi_M) <- c(m,n,z)
	    # n and z are not passed to the function? But it works...
      
      if (verbose) print("Calculates total number of female gametes")
      
      fecx<- array(0,dim=c(m,z))	# Number of female gametes produced by all the individuals of each genotype in each age class
	  fec<- array(0,dim=m)			# Number of female gametes produced by all the individuals of each genotype
      
	  for (k in 1 : m) {
		for (x in 1 : z) {
			fecx[k,x] <- sum(as.numeric(rpois(Nprime[k,i,x],phi_F[k,i,x])))	# This is the contribution of variation in reproductive success among individuals to genetic drift
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
			fecx[k,x] <- sum(as.numeric(rpois(Nprime[k,i,x],phi_M[k,i,x])))	# This is the contribution of variation in reproductive success among individuals to genetic drift
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
    
    
    
    if (verbose) print("Apply reproduction function")

    
    # If there is only one age-class, we must force the third dimension
    if (length(dim(phi_F))==2) dim(phi_F)[3] <- 1
    if (length(dim(phi_M))==2) dim(phi_M)[3] <- 1
    
    for (i in 1 : n) {
      
      if (save.res) {
        
        if (sum(Nprime[,i,])==0) { # To save computing time
          L[,i] = 0
          next
        } else {
          L[,i] <- repr.onelocus(Nprime,phi_F[,,,t],phi_M[,,,t],mu,i,l,m)
        }
        
      } else {
        
        if (sum(Nprime[,i,])==0) { # To save computing time
          L[,i,t] = 0
          next
        } else {
          L[,i,t] <- repr.onelocus(Nprime,phi_F[,,,t],phi_M [,,,t],mu,i,l,m)
        }
        
      }
    }
    
    
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
    
    
    if (verbose) print("Apply recruitment function")
    for (i in 1 : n) {
      if (save.res) {
        N[,i,] <- recr(N = array(Nprime[,i,], dim=c(m,z)),
                        S = array(S[,i], dim=c(m,1)),
                        m = m,
                        z = z,
                        kappa0 = kappa0[i,t],
                        recr.dd = recr.dd,
                        sexuality="monoecious")
      } else {
        N[,i,,t+1] <- recr(N = array(Nprime[,i,], dim=c(m,z)),
                            S = array(S[,i,t], dim=c(m,1)),
                            m = m,
                            z = z,
                            kappa0 = kappa0[i,t],
                            recr.dd = recr.dd,
                            sexuality = "monoecious")
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
  if (save.res==F) return(N)
}
