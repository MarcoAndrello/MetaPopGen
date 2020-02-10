# theta.global.biallelic<-function(N,s,t){
theta.global.biallelic<-function(N){
  
  # J'ai simplifié: la fonction prend en input la variable N (genotype X deme X classe d'age)
  
  
#   if (is.character(N)) {
#     dir.res.name <- paste(getwd(),N,sep="/")
#     name.file <- paste(N,"/N",t,".RData",sep="")
#     load(name.file)
#     
#     if (s == "M") {
#       N <- N_M
#     } else {
#       N <- N_F
#     }		
#   } else {
#     
#     if (s == "M") {
#       N <- N[[1]][,,,t]
#     } else {
#       N <- N[[2]][,,,t]
#     }
#   }
  
  #Define basic variables
  
  m<-dim(N)[1]#Number of genotype
  r<-dim(N)[2]#Number of demes
  l<-(sqrt(1+8*m)-1)/2#Number of alleles
  
  # Genotype arrays and number of individuals for each deme
  N_i <- apply(N,c(1,2),sum)
  n_i <- apply(N,2,sum)##Sample size in each deme
  
  # If any deme is empty (N_i==0), leave it out
  id.0 <- which(n_i==0)
  if (length(id.0)>=(r-1)) cat(length(id.0),"demes are empty. Stopping \n") # and gives error
  n_i <- n_i[-id.0]
  N_i <- N_i[,-id.0]
  r<-dim(N_i)[2]
  
  #Average sample size
  n_barre<-mean(n_i)
  
  #Variation of sample size
  n_c<-((r*n_barre)-((sum(n_i^2/(r*n_barre)))))/(r-1)
  
  
  # Calculate allele frequencies in each group
  p_i <- array(NA,dim=c(r,l))
  for (i in 1 : r) {
    p_i[i,] <- freq.all(N_i[,i])    
  }
  #???p_T	<- freq.all(N_T)
  
  #Average sample frequency of allele A
  p_barre<-sum(n_i*p_i[,1])/(r*n_barre)
  
  #Sample variance of allele A frequencies over populations
  s_2<-(sum(n_i*(p_i[,1]-p_barre)^2))/((r-1)*n_barre)
  
  # Calculate expected heterozygosities in each group and in the total population
#   H_S_i <- array(NA,dim=r)
#   for (i in 1 : r) {
#     H_S_i[i] <- 1 - sum(p_i[i,]^2)
#   } 
  # Dans l'article de W&C, il parlent de "observed proportion of individuals heterozygous"
  # Il faut donc l'hétérozygotie observée et pas l'hétérozygotie attendue:
  H_O_i <- N_i[2,] / n_i
  
  #Average heterozygote frequency for allele A
  #h_barre=sum(n_i*H_S_i)/(r*n_barre)
  h_barre=sum(n_i*H_O_i)/(r*n_barre)
  
  # Weir and Cockerham's (1984) equations 2, 3 and 4:
  a <- (n_barre/n_c)*(s_2-1/(n_barre-1)*(p_barre*(1-p_barre)-((r-1)/r)*s_2-(h_barre/4)))
  b <- (n_barre/(n_barre-1))*((p_barre*(1-p_barre))-(((r-1)/r)*s_2)-(((2*n_barre-1)/(4*n_barre))*h_barre))
  c <- h_barre/2
  
  theta=a/(a+b+c)
  return(theta)
}