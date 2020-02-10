# Pairwise FST 
theta.pairwise.biallelic = function(N,i,j,xi,xj,ti,tj,si,sj) {
  
  # Define 1-D array of genotype counts for the first and second group of individuals
  # Calculate number of individuals in each group 
  
  if (is.character(N)) {
    dir.res.name <- paste(getwd(),N,sep="/")
    
    # Load file for first group of individuals
    name.file <- paste(dir.res.name,"/N",ti,".RData",sep="")
    load(name.file)
    if (si == "M") {
      Ni <- N_M[,i,xi]
      n_i <- sum(N_M[,i,xi])
    } else {
      Ni <- N_F[,i,xi]
      n_i <- sum(N_F[,i,xi])
    }
    
    # Load file for second group of individuals
    name.file <- paste(dir.res.name,"/N",tj,".RData",sep="")
    load(name.file)
    if (sj == "M") {
      Nj <- N_M[,j,xj]
      n_j  <- sum(N_M[,j,xj])
    } else {
      Nj <- N_F[,j,xj]
      n_j  <- sum(N_F[,j,xj])
    }
    
    
  } else {
    if (si == "M") {
      Ni <- N[[1]][,i,xi,ti]
      n_i <- sum(N[[1]][,i,xi,ti])
    } else {
      Ni <- N[[2]][,i,xi,ti]
      n_i <- sum(N[[2]][,i,xi,ti])
    }
    
    if (sj == "M") {
      Nj <- N[[1]][,j,xj,tj]
      n_j  <- sum(N[[1]][,j,xj,tj])
    } else {
      Nj <- N[[2]][,j,xj,tj]
      n_j  <- sum(N[[2]][,j,xj,tj])
    }
  }
  
  
  # Calculate number of individuals in the total "population"
  # Population means the two groups together
  n_T	<- n_i + n_j
  
  # Calculate allele frequencies in each group and in the total "population"
  p_i	<- freq.all(Ni)
  p_j	<- freq.all(Nj)
  
  #Define basic variables
  
  m<-length(Nj)#Number of genotype
  r<-2#Number of demes
  l<-(sqrt(1+8*m)-1)/2#Number of alleles
  
  #Average sample size
  n_t=c(n_i,n_j)
  n_barre=mean(n_t)
  
  #Variation of sample size
  n_c<-((r*n_barre)-((sum(n_t^2/(r*n_barre)))))/(r-1)
  
  #Average sample frequency of allele A
  p_t<-rbind(p_i,p_j)
  p_barre=mean(c(p_i[1],p_j[1]))
  
  #Sample variance of allele A frequencies over populations
  s_2<-(sum(n_t*(p_t[,1]-p_barre)^2))/((r-1)*n_barre)
  
  # Calculate expected heterozygosities in each group and in the total "population"
  H_S_i	<- 1 - sum(p_i^2)
  H_S_j	<- 1 - sum(p_j^2)
  
  H_S_t=c(H_S_i,H_S_j)
  h_barre=sum(n_t*H_S_t)/(r*n_barre)
  
  
  a=(n_barre/n_c)*(s_2-1/(n_barre-1)*(p_barre*(1-p_barre)-((r-1)/r)*s_2-(h_barre/4)))
  b=(n_barre/(n_barre-1))*((p_barre*(1-p_barre))-(((r-1)/r)*s_2)-(((2*n_barre-1)/(4*n_barre))*h_barre))
  c=h_barre/2
  
  theta=a/(a+b+c)
  return(theta)
}