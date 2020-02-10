plot.pop.size <- function(name.dir,n,m,z,T_max){
  
       
  N.size <- array(NA, dim=c(n,T_max),dimnames=list(deme=c(1:n),time=c(1:T_max)))
  for (t in 1 : T_max){
    load(paste0(name.dir,"/N",t,".RData"))
    for (i in 1 : n){
      N.size[i,t] <- sum(N[,i,])
    }
  }

  ylim <- c(0,max(N.size))
  plot(c(1:T_max),N.size[1,],ylim=ylim,type="l",ylab="Population size",xlab="Time")
  for (i in 2 : n){
    matplot(c(1:T_max),N.size[i,],add=T,type="l")
  }
  
}

