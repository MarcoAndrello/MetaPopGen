# Convert survival, fecundity and N1 data from dataset to array form
metapopgen.input.convert.dioecious <- function(a) {

  # a must be a data.frame
  
m <- max(a$Genotype)
n <- max(a$Deme)
z <- max(a$Age)
tmax <- max(a$Time)

sigma_M <- array(NA,dim=c(m,n,z,tmax))
for (i in 1 : dim(a)[1]) sigma_M[a$Genotype[i],a$Deme[i],a$Age[i],a$Time[i]] <- a$Male_survival[i]

sigma_F <- array(NA,dim=c(m,n,z,tmax))
for (i in 1 : dim(a)[1]) sigma_F[a$Genotype[i],a$Deme[i],a$Age[i],a$Time[i]] <- a$Female_survival[i]

phi_M <- array(NA,dim=c(m,n,z,tmax))
for (i in 1 : dim(a)[1]) phi_M[a$Genotype[i],a$Deme[i],a$Age[i],a$Time[i]] <- a$Male_fecundity[i]

phi_F <- array(NA,dim=c(m,n,z,tmax))
for (i in 1 : dim(a)[1]) phi_F[a$Genotype[i],a$Deme[i],a$Age[i],a$Time[i]] <- a$Female_fecundity[i]

N1_M <- array(NA,dim=c(m,n,z))
for (i in 1 : dim(a)[1]) if (a$Time[i]==1) N1_M[a$Genotype[i],a$Deme[i],a$Age[i]] <- a$N1_M[i]

N1_F <- array(NA,dim=c(m,n,z))
for (i in 1 : dim(a)[1]) if (a$Time[i]==1) N1_F[a$Genotype[i],a$Deme[i],a$Age[i]] <- a$N1_F[i]

return(list(N1_M,N1_F,sigma_M,sigma_F,phi_M,phi_F))

}

