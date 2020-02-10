
separation_loci<-function(N_F,N_M,list_name,List_gene){
    
    n<-dim(N_F)[2]
    m<-dim(N_F)[1]
    N<-N_F+N_M
    N_z<-apply(N,c(1,2),sum)
    
    mA<-List_gene[1]*(List_gene[1]+1)/2
    mB<-List_gene[2]*(List_gene[2]+1)/2
    N_A<-array(NA,dim=c(mA,n))
    N_B<-array(NA,dim=c(mB,n))
    
    name_A<-c()
    name_B<-c()
    for (k in 1:m){
        name_A<-c(name_A,paste0(substr(list_name[k],1,2),substr(list_name[k],6,7)))
        name_B<-c(name_B,paste0(substr(list_name[k],3,4),substr(list_name[k],8,9)))
    }
    k<-1
    while(length(name_A)!=mA){
        kk<-1
        while(kk<=length(name_A)){
            if (name_A[k]==name_A[kk]){
                name_A<-name_A[-kk]
            }
            kk<-kk+1
        }
        k<-k+1
    }
    
    k<-1
    while(length(name_B)!=mB){
        kk<-1
        while(kk<=length(name_B)){
            if (name_B[k]==name_B[kk]){
                name_B<-name_B[-kk]
            }
            kk<-kk+1
        }
        k<-k+1
    }
    
    dimnames(N_A)<-list(name_A,c(1:n))
    dimnames(N_B)<-list(name_B,c(1:n))
    
    
    
    for (k in 1:mA){
        A<-which(rownames(N_A)[k]== paste0(substr(list_name,1,2),substr(list_name,6,7)))
        N_A[k,]<-apply(N_z[A,],2,sum)
    }
    
    for (k in 1:mB){
        B<-which(rownames(N_B)[k]== paste0(substr(list_name,3,4),substr(list_name,8,9)))
        N_B[k,]<-apply(N_z[B,],2,sum)
    }
    return(c(N_A,N_B))
}

    
