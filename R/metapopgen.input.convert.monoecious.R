# Convert survival, fecundity and N1 data from dataset to array form
metapopgen.input.convert.monoecious <- function(a) {

  # a must be a data.frame
  
m <- max(a$Genotype)
n <- max(a$Deme)
z <- max(a$Age)
tmax <- max(a$Time)

sigma <- array(NA,dim=c(m,n,z,tmax))
for (i in 1 : dim(a)[1]) sigma[a$Genotype[i],a$Deme[i],a$Age[i],a$Time[i]] <- a$Survival[i]

phi_M <- array(NA,dim=c(m,n,z,tmax))
for (i in 1 : dim(a)[1]) phi_M[a$Genotype[i],a$Deme[i],a$Age[i],a$Time[i]] <- a$Male_fecundity[i]

phi_F <- array(NA,dim=c(m,n,z,tmax))
for (i in 1 : dim(a)[1]) phi_F[a$Genotype[i],a$Deme[i],a$Age[i],a$Time[i]] <- a$Female_fecundity[i]

N1 <- array(NA,dim=c(m,n,z))
for (i in 1 : dim(a)[1]) if (a$Time[i]==1) N1[a$Genotype[i],a$Deme[i],a$Age[i]] <- a$N1[i]

return(list(N1,sigma,phi_M,phi_F))

}
