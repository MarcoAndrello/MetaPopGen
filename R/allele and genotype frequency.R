# 
# #Calcul allele and genotypes frenquencies in a population  
# 
# Freq.genotype<- function(N){
#   m<-length(N)
#   P<-array(NA,m)
#   ntot<-sum(N)
#   for(i in 1:m){
#     P[i]<-N[i]/ntot
#   }
#   return(P)
# }
# 
# Freq.allele<-function(N,allele_vec,index){
#   P <- array(NA,dim=sum(allele_vec))
#   GENO <- create.genotype.matrix(allele_vec, index)
#   ntot <- sum(N)
#   k <- 1
#   for (i.locus in 1 : length(allele_vec)){
#     for(i.allele in 1 : allele_vec[i.locus]){
#       pos <- Position(GENO, i.locus, i.allele)
#       P[k] <- 0.5*sum(N[pos])/ntot
#       k <- k + 1
#     }
#   }
#   return(P)
# }
# 
# 
# #gives the position in the Genotype matrix GENO, with the number of the gene and the number
# #of the allele 
# 
# Position<-function(GENO,Gene,allele){
#   i<-Gene
#   pos<-c()
#   for(k in 1:length(GENO[,1])){
#     if(GENO[k,i]==allele){
#       pos<-c(pos,k)
#     }
#     if (GENO[k,i+2]==allele){
#       pos<-c(pos,k)
#     }
#   }
#   return(pos)
# }
# 

