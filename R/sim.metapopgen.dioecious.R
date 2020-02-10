sim.metapopgen.dioecious <- function(input.type,demographic.data,
                                     N1_M,N1_F,
                                     sigma_M,sigma_F,
                                     phi_M,phi_F,mu,
                                     delta,
                                     recr.dd="settlers",kappa0,
                                     T_max,save.res=F,save.res.T=seq(1,T_max),
                                     verbose=F) {
  
  
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
  
  
  # Reproduction
  repr.onelocus <-
    function(Nprime_M,Nprime_F,phi_F,phi_M,mu,i,l,m) {
      	
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
			fecx[k,x] <- sum(as.numeric(rpois(Nprime_F[k,i,x],phi_F[k,i,x])))	# This is the contribution of variation in reproductive success among individuals to genetic drift
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
			  fecx[k,x] <- sum(as.numeric(rpois(Nprime_M[k,i,x],phi_M[k,i,x])))	# This is the contribution of variation in reproductive success among individuals to genetic drift
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
               sigma <- settler.survival(Ntot,kappa0)
             },
             
             # Dependence on adult density
             adults = {
               Ntot <- sum(N_M[,1:(z-1)] + N_F[,1:(z-1)])
               Stot <- sum(S+S_othersex)
               Recr <- kappa0 - Ntot
               
               if (Recr <= 0){
                 sigma <- 0
               } else {
                 sigma <- Recr / Stot
               }
               
               if (sigma > 1) sigma <- 1
             },
             
             # No-match: error
             stop("Unknown value for argument recr.dd. Valid values: 'settlers', 'adults'")
      )
      
      # Use recruitment probability to calculate the number of recruits
      R <- array(0,dim=m)
      for (k in 1 : m) {
        R[k] <- rbinom(1,S[k],sigma)
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
    
    print(t)
    
    
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
          LL <- repr.onelocus(Nprime_M,Nprime_F,phi_F[,,,t],phi_M[,,,t],mu,i,l,m)
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
          LL <- repr.onelocus(Nprime_M,Nprime_F,phi_F[,,,t],phi_M[,,,t],mu,i,l,m)
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
        N_M[,i,1] <- recr(S_M[,i],S_F[,i],Nprime_M[,i,],Nprime_F[,i,],m,kappa0[i,t],recr.dd)  # Pass the abundance of other classes too
        N_F[,i,1] <- recr(S_F[,i],S_M[,i],Nprime_M[,i,],Nprime_F[,i,],m,kappa0[i,t],recr.dd)
      } else {
        N_M[,i,1,t+1] <- recr(S_M[,i,t],S_F[,i,t],Nprime_M[,i,],Nprime_F[,i,],m,kappa0[i,t],recr.dd)
        N_F[,i,1,t+1] <- recr(S_F[,i,t],S_M[,i,t],Nprime_M[,i,],Nprime_F[,i,],m,kappa0[i,t],recr.dd)
      }
    }
    
    if (verbose) print("Calculates N at t+1")
    for (i in 1 : n){
      for (x in 1 : z) {
        if (x == 1) next
        for (k in 1 : m) {
          if (save.res) {
            N_M[k,i,x] <- Nprime_M[k,i,x-1]
            N_F[k,i,x] <- Nprime_F[k,i,x-1]
          } else {
            N_M[k,i,x,t+1] <- Nprime_M[k,i,x-1]
            N_F[k,i,x,t+1] <- Nprime_F[k,i,x-1]
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
