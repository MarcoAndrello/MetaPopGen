# Multilocus functions

# expand.grid.unique
# def_genotype.name.locus
# genotype.index.multilocus
# create.meiosis.matrix
# create.mat_geno_to_index_mapping
# initialize.multilocus
# Type_gamete
# create.multilocus.survival


# Function expand.grid.unique found on stack.overflow from user Ferdinand.kraft
# https://stackoverflow.com/questions/17171148/non-redundant-version-of-expand-grid
# It is used to couple the alleles of the same locus to form unique genotypes
expand.grid.unique <- function(x, y, include.equals=TRUE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i) {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

# Define genotype names for each locus
def_genotype.name.locus <- function(allele_vec){
  genotype.name.locus <- list() # For each locus (dimension), define the names of the single-locus genotypes
  for (loc in 1 : nloc) {
    alleles <- vector("character",allele_vec[loc])
    for (i in 1 : allele_vec[loc]){
      alleles[i] <- paste0(LETTERS[loc],i)
    }
    genotype.name.locus[[loc]] <- apply(
      expand.grid.unique(alleles, alleles),
      1,paste0,collapse="")
  }
  genotype.name.locus
}

genotype.index.multilocus <- function(allele_vec,r,method){
    
    # .................... Initial checkings
    # Number of loci
    if (!is.numeric(allele_vec)) stop("allele_vec must be a numeric vector")
    if (length(allele_vec) > 26) stop("The length of allele_vec must be between 2 and 26")
    if (length(allele_vec)<2) stop("The length of allele_vec must be between 2 and 26")
    
    # Recombination rate
    if (length(r) > 1) stop("The recombination rate must be a single number")
    if (r < 0 | r > 0.5) stop("The recombination rate must be 0 < r <= 0.5")
    if ( r < 0.5 & length(allele_vec) > 2 ) stop("If number of loci is >2, then recombination rate must be 0.5")
    
    if (r < 0.5 & method != "gamete") stop("Method must be 'gamete' for two loci with r < 0.5")
    # ..................... End initial checkings
    
    # Given a vector of number of alleles per locus, the number of (multilocus) gametotypes is prod(allele_vec)
    l <- prod(allele_vec)
    
    # Set up the matrix that will contain the index of genotypes for each pair of gametotypes
    y <- array(NA,dim=c(l,l))
    
    # Two methods
    # 1) gamete-based: gives a a 2*2 matrix
    # 2) locus-based: gives a nloc-dimensional array
    
    
    # 1) gamete-based
    if ( method == "gamete") {
        # Set up the diagonal of the matrix: it corresponds to homozygote multilocus genotypes 
        diagy <- array(NA,dim=l)
        
        # The first element of the diagonal is always 1
        # The subsequent elements are found by summing the index of the diagonal of the previous (i-1) line
        # plus the number of new genotypes that the previous (i-1) line contains
        # Ex: 
        #       A1B1 A2B1 A1B2 A2B2
        # A1B1    1    2    3    4
        # A2B1    2    5    6    7
        # A1B2    3    6    8    9
        # A2B2    4    7    9   10
        # The first line contains 4 new genotypes, the second line contains 3 new genotypes, the third contains 2 new genotypes...
        # That is, line i contains (l-i) new genotypes. If we do not count the diagonal, then it is (l-(i-1)) new genotypes.
        # Thus when working on line i, the previous line contains (l-((i-1)-1)) genotypes = (l-(i-2)) genotypes
        diagy[1] <- 1
        for (i in 2 : l) {
            diagy[i] <- diagy[i-1] + (l-(i-2))
        }
        diag(y)<-diagy
        
        #The the subdiagonal elements. The number of new genotypes of line i is (l-i)
        # The matrix is symmetric
        for (i in 1 : (l-1)) {
            y[i,i:l] <- seq(y[i,i],(y[i,i]+l-i))
            y[i:l,i] <- seq(y[i,i],(y[i,i]+l-i))
        }
        
        # Set dimension names, for <=26 loci
        nloc <- length(allele_vec)
        alleles_list <- list()
        for (loc in 1 : nloc) { 
            alleles <- vector("character",allele_vec[loc])
            for (i in 1 : allele_vec[loc]){
                alleles[i] <- paste0(LETTERS[loc],i)
            }
            alleles_list[[loc]] <- alleles
        }
        # ....These two lines lead to former loci varying more slowly
        alleles_list <- expand.grid( rev(alleles_list) )
        alleles_list <- alleles_list[,c(nloc:1)]
        # ....
        alleles_list <- apply(alleles_list,2,as.character)
        
        dn <- apply(alleles_list,1,paste0,collapse="")
        dimnames(y) <- list(dn,dn)
        
    }
    
    # 2) Locus-based
    if ( method == "locus") {    

        nloc <- length(allele_vec)
        # Number of single-locus genotypes per locus
        m.locus <- vector()  
        for (i in 1 : nloc) m.locus[i] <- allele_vec[i] * (allele_vec[i] + 1) / 2
        
        # Define matrix
        y <- array(c(1:prod(m.locus)), rev(m.locus)) # rev is to ensure that former loci vary more slowly
        y <- aperm(y,c(nloc:1))                      # To ensure that former loci vary more slowly
        
        # Define dimnames
        dimnames(y) <- def_genotype.name.locus(allele_vec)
        
    } 
    return(y)
}


create.meiosis.matrix <- function(index_matr, allele_vec, r=0.5, mu, method, return_submatrices = F){
    # Creation of the meiosis matrix
    # This matrix gives the probabilities for each gamete type as a function of the parental genotypes and mutation,
    
    # index_matr    n-dimensional array of indices of multilocus genotypes
    #               (obtained from genotype.index.multilocus in the nloc-dimensional form for the locus-based method or in the 2D form for the gamete-based method)
    # allele_vec    vector of number of alleles at each locus
    # r = 0.5       recombination rate
    # mu            Vector giving the mutation rate at each locus
    # method can be "gamete" or "genotype" 
    
    if (method == "locus" ) {
        if(any(allele_vec>9)) stop("More than 10 alleles per locus is not implemented yet")
        
        nloc <- length(allele_vec) # Number of loci
        num.genotypes <- max(index_matr)       # Number of multi-locus genotypes
        num.gametes <- prod(allele_vec)      # Number of gametotypes
        
        # Set the names of the multilocus gametotypes (for <=26 loci) in the good order
        # These will be the rownames of the meiosis_matrix
        alleles_list <- list()
        for (i.locus in 1 : nloc) { 
            alleles <- vector("character",allele_vec[i.locus])
            for (i in 1 : allele_vec[i.locus]){
                alleles[i] <- paste0(LETTERS[i.locus],i)
            }
            alleles_list[[i.locus]] <- alleles
        }
        # ....These two lines lead to former loci varying more slowly
        alleles_table <- expand.grid( rev(alleles_list) )
        alleles_table <- alleles_table[,c(nloc:1)]
        # ....
        alleles_table  <- apply(alleles_table ,2,as.character)
        gametotype.names <- apply(alleles_table ,1,paste0,collapse="")
        
        # Set the names of the multilocus genotypes (for <=26 loci) in the good order
        # These will be the colnames of the meiosis_matrix
        genotype.names <- rep(NA,num.genotypes)
        allele_locus_1 <- list() # list containing, for each genotype, the first allele of each locus
        allele_locus_2 <- list() # list containing, for each genotype, the second allele of each locus
        for (i.geno in 1 : num.genotypes) {
            pos.geno <- which(index_matr==i.geno,arr.ind=T)
            allele_locus_1[[i.geno]] <- allele_locus_2[[i.geno]] <- vector()
            for (i.locus in 1 : length(pos.geno)) {
                allele_locus_1[[i.geno]][i.locus] <- substr( dimnames(index_matr)[[i.locus]][pos.geno[i.locus]],1,2) # Only for < 10 alleles
                allele_locus_2[[i.geno]][i.locus] <- substr( dimnames(index_matr)[[i.locus]][pos.geno[i.locus]],3,4) # Only for < 10 alleles
            }
            allele_multi_1 <- paste0(allele_locus_1[[i.geno]],collapse="")
            allele_multi_2 <- paste0(allele_locus_2[[i.geno]],collapse="")
            genotype.names[i.geno] <- paste0(allele_multi_1,"/",allele_multi_2)
        }
        
        # Set up the list containing, for each locus, the contribution of each allele to each gametotype (as a matrix)
        contr_matr <- list() #  Each element of the list corresponds to a locus, and is a matrix with allele_vec[i.locus] columns
        for (i.locus in 1 : nloc) { 
            contr_matr[[i.locus]] <- array(0,c(num.gametes,allele_vec[i.locus]))
            dimnames(contr_matr[[i.locus]]) <- list(gametotype = gametotype.names, allele = alleles_list[[i.locus]])
            for (i.allele in 1 : allele_vec[i.locus]) {
                id <- grep(colnames(contr_matr[[i.locus]])[i.allele],rownames(contr_matr[[i.locus]])) # Indices of the gametotypes containing the focal allele
                # Without mutation:
                (contr_matr[[i.locus]][id, i.allele] <- 1)
                # With mutation (has not undergone *extensive* checking):
                contr_matr[[i.locus]][id, i.allele] <- (1-mu[i.locus])
                contr_matr[[i.locus]][-id, i.allele] <- mu[i.locus] / (allele_vec[i.locus]-1) # Mutation to any other allele with the same probability
            }
        }
        names(contr_matr) <- paste0("Locus",LETTERS[1:nloc])
        
        # Set up meiosis_matrix
        meiosis_matrix <- matrix(0, nrow = num.gametes, ncol = num.genotypes)
        dimnames(meiosis_matrix) <- list(gamete = gametotype.names, genotype = genotype.names)
        
        for (i.geno in 1 : num.genotypes) {
            contr_geno <- array(NA,c(num.gametes,nloc)) # This contains the contribution of each locus
            dimnames(contr_geno) <- list(gametotype = gametotype.names, Locus = paste0("Locus",LETTERS[1:nloc])) # Commented to speed up computation; decomment for debugging
            for (i.locus in 1 : nloc) {
                contr_geno[,i.locus] <-                    
                    0.5 *                                                            # Add the contribution of the first allele at this locus; 0.5 because of diploidy 
                    contr_matr[[i.locus]][ , allele_locus_1[[i.geno]][i.locus] ] +   # This vector gives the unitary contribution of the first allele at this locus of this genotype to all gametotypes
                    0.5 *                                                            # Add the contribution of the second allele at this locus; 0.5 because of diploidy 
                    contr_matr[[i.locus]][ , allele_locus_2[[i.geno]][i.locus] ]     # This vector gives the unitary contribution of the second allele at this locus of this genotype to all gametotypes
            }
            meiosis_matrix[,i.geno] <- apply(contr_geno,1,prod)                      # prod to filter only the possible combinations
        }
        
        return(meiosis_matrix)
    }
    
    if (method == "gamete") {
        num.loci <- length(allele_vec) # Number of loci
        num.genotypes <- max(index_matr)    # Number of multi-locus genotypes
        num.gametes <- prod(allele_vec)   # Number of multi-locus alleles
        
        # Create GENO matrix to know which allele comes from which genotype
        # GENO <- array(NA,dim=c(num.genotypes,2*num.loci))    # Matrix: for each genotype (lines), gives which allele at each locus
        # #(locus 1: 1st and 2nd column; locus 2: rd and 4th column; etc)
        # names <- c()
        # for(i in 1 : num.gametes){
        #   for(ii in 1 : num.gametes){
        #     if(is.na(GENO[index_matr[i,ii],1])){
        #       w <- 0
        #       for (h in 1 : num.loci){
        #         a <- c(substr(rownames(index_matr)[i],2*h,2*h),substr(colnames(index_matr)[ii],2*h,2*h))
        #         GENO[index_matr[i,ii],c(h+w,2*h)] <- as.numeric(a)
        #         w <- w+1
        #       }
        #       names <- c(names,paste0(rownames(index_matr)[i],"/",colnames(index_matr)[ii]))
        #     }
        #   }
        # }
        #dimnames(GENO) <- list(names,rep(1:2,num.loci))
        GENO <- array(NA,c(num.genotypes,4)) # Works only for 2 loci. 4 is (num.loci*2)
        rownames(GENO) <- rep("",num.genotypes)
        i.geno <- 1
        for (i.geno in 1 : num.genotypes){
            indices <- which(index_matr == i.geno, arr.ind = T)
            if (dim(indices)[1] > 1) indices <- indices[dim(indices)[1],]
            id.row <- indices[1]
            row <- rownames(index_matr)[id.row]
            GENO[i.geno,1] <- as.numeric(substr(row,2,2))
            GENO[i.geno,3] <- as.numeric(substr(row,4,4))
            id.col <- indices[2]
            col <- colnames(index_matr)[id.col]
            GENO[i.geno,2] <- as.numeric(substr(col,2,2))
            GENO[i.geno,4] <- as.numeric(substr(col,4,4))
            rownames(GENO)[i.geno] <- paste0(row,"/",col)
        }
        # colnames(GENO) <- c("gam1_A","gam1_B","gam2_A","gam2_B")
        
        
        # Create RECOMB matrix
        RECOMB <- array(0, dim=c(num.gametes,num.genotypes)) # Matrix giving the probability of production of each gameto-type (lines) for each genotype (columns)
        dimnames(RECOMB) <- list(colnames(index_matr), rownames(GENO))
        # if (length(allele_vec)==2){
        # This gives the gameto-type: which allele for each locus for each line of RECOMB
        # Ensures that first locus varies more slowly, as in GENO
        TYPE <- expand.grid(1:allele_vec[2], 1:allele_vec[1])
        TYPE <- TYPE[,c(2:1)]
        colnames(TYPE) <- LETTERS[1:num.loci] ## I added this
        rownames(TYPE) <- rownames(RECOMB) ### I added this 07/07/2020
        for (k in 1 : num.genotypes){
            genotype<-GENO[k,]
            for(j in 1 : num.gametes){
                if ((genotype[1]== TYPE[j,1]) & (genotype[3]==TYPE[j,2])){
                    RECOMB[j,k] <- RECOMB[j,k] + ((1/2)*(1-r))
                }
                if ((genotype[1]== TYPE[j,1]) & (genotype[4]==TYPE[j,2])){
                    RECOMB[j,k] <- RECOMB[j,k] + (r/2)
                }
                if ((genotype[2]== TYPE[j,1]) & (genotype[3]==TYPE[j,2])){
                    RECOMB[j,k] <- RECOMB[j,k] + (r/2)
                }
                if ((genotype[2]== TYPE[j,1]) & (genotype[4]==TYPE[j,2])){
                    RECOMB[j,k] <- RECOMB[j,k] + ((1/2)*(1-r))
                }
            }
        }
        # } else {
        # p<-2^num.loci
        # 
        # for(k in 1:num.genotypes){
        #   if(substr(rownames(GENO)[k],1,2*num.loci)==substr(rownames(GENO)[k],2*num.loci+2,4*num.loci+1)){ # If it is a homozygote
        #     P<-which(colnames(index_matr)==substr(rownames(GENO)[k],1,2*num.loci)) # .. it produces only one gametotype, equal to either gamete forming the genotype
        #     RECOMB[P,k]<-1
        #   }
        #   else {
        #     alleles_list <- list()
        #     w<-0
        #     for (h in 1 : num.loci) {
        #       alleles <- vector("character",2)
        #       alleles[1] <- paste0(LETTERS[h],GENO[k,h+w])
        #       alleles[2]<-paste0(LETTERS[h],GENO[k,h+1+w])
        #       w<-w+1
        #       alleles_list[[h]] <- alleles
        #     }
        #     alleles_list <- expand.grid(alleles_list)
        #     alleles_list<-apply(alleles_list,1,paste0,collapse="")
        #     for (h in 1: length(alleles_list)){
        #       P<-which(rownames(index_matr)==alleles_list[h])
        #       RECOMB[P,k]<-RECOMB[P,k]+1/p
        #     }
        #   }
        # }
        # }
        
        MU <- array(NA,dim=c(num.gametes,num.gametes)) # Mutation probabilities between gametes
        gametes <- rownames(index_matr)
        dimnames(MU) <- list(gamete_source = gametes, gamete_mutated = gametes)
        for(i in 1 : num.gametes){
            for(j in 1 : num.gametes){
                p <- rep(NA,num.loci)
                for (i.locus in 1 : num.loci){
                    if( substr(gametes[i],2*i.locus,2*i.locus) == substr(gametes[j],2*i.locus,2*i.locus) ){ # Check if the locus is homozygous for this (i,j) combination of gametes
                        p[i.locus] <- 1 - mu[i.locus]
                    }
                    else {
                        p[i.locus] <- mu[i.locus] / (allele_vec[i.locus]-1) # Mutation to a different allele is mu/(number of alleles minus 1)
                    }
                }
                MU[i,j] <- prod(p)
            }
        }
        
          PROBA<-array(0,dim=c(num.gametes,num.genotypes))
          dimnames(PROBA)<-list(type=rownames(index_matr),genotype_parent=rownames(GENO))
          
          for (j in 1 : num.genotypes){
            type <- c()
            p <- c()
            for(i in 1 : num.gametes){
              if(RECOMB[i,j] == 0) {
                next
              } else {
                type <- c(type,i)
                p <- c(p,RECOMB[i,j])
              }
            }
            p_i <- c()
            for(i in 1 : length(type)){
              p_i <- p[i]*MU[type[i],]
              PROBA[,j] <- PROBA[,j] + c(p_i)
            }
          }

          # Slower but clearer (07/07/2020) -- try transforming the inner for loop with matrix multiplications
          # PROBA<-array(0,dim=c(num.gametes,num.genotypes))
          # dimnames(PROBA)<-list(type=rownames(index_matr),genotype_parent=rownames(GENO))
          # for (j in 1 : ncol(PROBA)) {
          #   for (k1 in 1 : nrow(PROBA)) {
          #     for (k2 in 1 : nrow(MU)) {
          #       PROBA[k1,j] <- PROBA[k1,j] +  MU[k2,k1] * RECOMB[k2,j]
          #     }
          #   }
          # }

        
        if (return_submatrices) return(list(RECOMB = RECOMB, MU = MU, PROBA = PROBA)) else return(PROBA)
        
    }
    
}


create.mat_geno_to_index_mapping <- function(allele_vec, Proba, index_matr) {
    
    # mat_geno_to_index_mapping : una tabella di corrispondenza tra la matrice mat_geno (2*2, femmine*maschi) e il numero del genotipo multilocus dato da index_matr, che ? quello in cui i genotipi sono salvati in metapopgen
    # viene creata offline all'inizio delle simulazioni
    # Sarebbe valida anchebcon il method "old" ?
    
    nloc <- length(allele_vec) # nloc is the number of loci
    l <- dim(Proba)[1] # l is the number of gametotypes
    mat_geno_to_index_mapping <- array(NA,c(l,l)) 
    dimnames(mat_geno_to_index_mapping) <- list(female=rownames(Proba), male=rownames(Proba))
    pos.locus <- vector()
    
    # la logica di quello che segue ?:
    # Per ogni combinazione di gameti devo trovare il genotipo multilocus che esce, quindi
    # per ogni locus devo trovare il genotipo di quel locus: guardo come si chiamano i gameti a quel locus (substr...), e ne deduco il nome del genotipo (paste0...)
    # leggo i nomi dei genotipi di quel locus (dimnames(index_matr))
    # e poi faccio un match
    # il comando index_matr[matrix(pos.locus,ncol=nloc)] ? un'indicizzazione che equivale a index_matr[pos.locus[1],pos.locus[2],..], ma indipendente dal numero di loci !
    for (i in 1 : l) {               # l is dim(mat_geno)[1]
        for (j in 1 : l) {           # l is dim(mat_geno)[2]
            # Note that mat_geno is always a 2*2 matrix (female*male) and always square, so dim(mat_geno)[1] = dim(mat_geno)[2]
            for (i.locus in 1 : nloc) {
                all.female <- as.numeric( substr(rownames(Proba)[i],(i.locus*2),(i.locus*2)) ) # Valid only for <10 alleles, otherwise we have three characters per locus
                all.male   <- as.numeric( substr(rownames(Proba)[j],(i.locus*2),(i.locus*2)) ) # Valid only for <10 alleles, otherwise we have three characters per locus
                name.genotype  <- paste0(LETTERS[i.locus], min(all.female,all.male), LETTERS[i.locus], max(all.female,all.male))
                name.genotypes <- dimnames(index_matr)[[i.locus]]
                pos.locus[i.locus] <- match(name.genotype, name.genotypes) # pos_locus[i.locus] indica a quale genotipo (monolocus) corrisponde la combinazione i,j di mat_geno 
            }
            mat_geno_to_index_mapping[i,j] <- index_matr[matrix(pos.locus,ncol=nloc)]
        }
    }
    return(mat_geno_to_index_mapping)
}

# initialize.multilocus
initialize.multilocus <- function(allele_vec, r, mu=rep(0,length(allele_vec)), n, z, kappa0, init.state="fixed", sexuality=NULL, return_submatrices = F) {
  
  # allele_vec    Vector of number of alleles at each locus
  # r = 0.5       Recombination rate
  # mu            Vector giving the mutation rate at each locus
  # n             Number of demes
  # z             Number of age.classes
  # kappa0        Number of individual per deme (can be a vector of length n)
  # init.state    Can be "fixed" or "random"
  # allele_vec=c(2,3)
  # r=0.5
  # mu=rep(0.01,length(allele_vec))
  # n=3
  # z=1
  # kappa0=100
  # init.state="fixed"
  # return_submatrices = F
  
  # Checking
  
  # Number of loci
  if (!is.numeric(allele_vec)) stop("allele_vec must be a numeric vector")
  if (length(allele_vec) > 26) stop("The length of allele_vec must be between 2 and 26")
  if (length(allele_vec)<2) stop("The length of allele_vec must be between 2 and 26")
  # Length of mu and allele_vec
  if (length(allele_vec) != length(mu)) stop("allele_vec and mu must have the same length")
  # Recombination rate
  if (length(r) > 1) stop("The recombination rate must be a single number")
  if (r < 0 | r > 0.5) stop("The recombination rate must be 0 < r <= 0.5")
  if ( r < 0.5 & length(allele_vec) > 2 ) stop("If number of loci is >2, then recombination rate must be 0.5")
  # Method
  if (r < 0.5) method <- "gamete" else method <- "locus"
  # Sexuality
  if (is.null(sexuality)) {
    cat("Assuming monoecious, as sexuality is not specified \n")
    sexuality <- "monoecious"
  }
  if (all(sexuality != "monoecious" && sexuality != "dioecious")) {
    stop("Sexuality must be either 'monoecious' or 'dioecious' \n")
  }
  
  cat("Generating genotype index \n")
  index_matr <- genotype.index.multilocus(allele_vec,r,method)
  cat("Generating meiosis matrix \n")
  ret_meiosis_matrix <- create.meiosis.matrix(index_matr, allele_vec, r, mu, method, return_submatrices)
  if (return_submatrices && method == "gamete") {
    meiosis_matrix <- ret_meiosis_matrix$PROBA
    RECOMB <- ret_meiosis_matrix$RECOMB
    MU <- ret_meiosis_matrix$MU
  } else {
    meiosis_matrix <- ret_meiosis_matrix
  }
  if (method == "locus") {
    cat("Generating genotype mapping \n")
    mat_geno_to_index_mapping <- create.mat_geno_to_index_mapping(allele_vec, meiosis_matrix, index_matr)
  } else {
    mat_geno_to_index_mapping <- NULL
  }
  
  
  cat("Generating N1 \n")
  if (sexuality == "monoecious") n.sexes <- 1 else n.sexes <- 2
  m <- dim(meiosis_matrix)[2] # Number of multilocus genotypes
  if (init.state == "fixed") {
    N1 <- array(round(kappa0/(m*z*n.sexes)),dim=c(m,n,z))
  }
  if (init.state == "random") {
    # Random in all demes
    N1 <- array(0,dim=c(m,n,z))
    for (deme in 1 : n){
      for (age in 1 : z) {
        for (genotype in 1 : m) {
          N1[genotype,deme,age] <- round(runif(1)*((kappa0*2)/(m*z*n.sexes))) # Expected value of runif is 0.5, so I multiply kappa0 by 2. What is k??
        }
      }
    }
  }
  
  # Adding dimanames to N1
  dimnames(N1) <- list(genotype = colnames(meiosis_matrix),
                       deme = c(1:n),
                       age = c(1:z))
  
  cat("Done")
  out <- list(allele_vec = allele_vec,
              r = r,
              mu = mu,
              m = dim(meiosis_matrix)[2],
              n = n,
              z = z,
              kappa0 = kappa0,
              index_matr = index_matr,
              meiosis_matrix = meiosis_matrix,
              mat_geno_to_index_mapping = mat_geno_to_index_mapping,
              N1 = N1,
              method = method)
  if (sexuality == "dioecious") {
    names(out)[which(names(out) == "N1")] <- "N1_F"
    out$N1_M <- out$N1_F
  }
  if (return_submatrices && method == "gamete") {
    out$RECOMB <- RECOMB
    out$MU <- MU
  }
  return(out)
}


# Type gamete define the number of gametes of each gametotype given by a genotype
Type_gamete <- function(fec,Proba){
    
    l <- dim(Proba)[1] # Number of gametotypes
    m <- dim(Proba)[2] # Number of genotypes
    
    Gamete <- array(0,dim = c(l,m))
    for (k in 1 : m){
        err <- try (as.vector(rmultinom(1, fec[k], Proba[,k])),silent=T) # meiosis
        if (class(err)=="try-error") {
            prob.mvr <- fec[k] * Proba[,k]				# Vector of means of the multivariate normal distribution
            var.mvr <- fec[k]* Proba[,k] * (1-Proba[,k])	# Vector of variances of the multivariate normal distribution
            sigma.mvr <- diag(var.mvr, l)			# Variance-covariance matrix of the multivariate normal distribution
            
            for (i.mvr in 1 : l){
                for (j.mvr in 1 : l) {
                    if (i.mvr == j.mvr) next
                    sigma.mvr[i.mvr,j.mvr] <- -fec[k] * Proba[i.mvr,k] * Proba[j.mvr,k]
                }
            }
            err<-as.vector(round(loinorm(1,prob.mvr,sigma.mvr)))
        }
        Gamete[,k]<-err
    }
    type<-array(0,l)
    for (j in 1:l){
        type[j]<-sum(Gamete[j,])
    }
    return(type)
}

# Create multilocus survival
create.multilocus.rate <-function(allele_vec, init.par, n.sel, x, xi, rate.max, omega=1, T_max){
    # This function was used in Andrello et al, submitted to Mol Ecol Res
    # Checking
    # Number of loci
    if (!is.numeric(allele_vec)) stop("allele_vec must be a numeric vector")
    if (length(allele_vec) > 26) stop("The length of allele_vec must be between 2 and 26")
    if (length(allele_vec)<2) stop("The length of allele_vec must be between 2 and 26")
    # nsel
    if (n.sel > length(allele_vec)) stop(paste0("n.sel = ",n.sel," cannot be larger than the number of loci, length(allele_vec) = ",length(allele_vec)))
    if (any(allele_vec[1:n.sel]>2)) stop(paste0("Some of the ",n.sel," loci under selection have more than two alleles. The loci under selection must all be biallelic"))
    # xi
    if (length(xi)!=(2*n.sel + 1)) stop(paste0("xi is of length ",length(xi),", but it must be of length ",(2*n.sel + 1),
                                               " (number of phenotypes with ",n.sel," biallelic loci under selection)"))
    # omega
    if (omega<=0) stop(paste0("omega = ",omega," must be strictly positive"))
    
    # Parameters
    nloci <- length(allele_vec) 
    l <- dim(init.par$meiosis_matrix)[1] # Number of multilocus gametes
    m <- init.par$m # Number of genotypes
    n <- init.par$n # Number of demes
    z <- init.par$z
    K = 2*n.sel + 1 # 2 is because all loci under selection are biallelic. Next version of the function could consider cases of more than 2 alleles at selected loci
    
    # Initialize survival
    rate_OUT <- array(0,dim=c(m,n,z,T_max))
    # Defining the rate as a function of genotype class and environmental variable x using an exponential function
    # See equation in the manuscript
    Rate <- array(NA,dim=c(n,K))
    for(i in 1 : K) {
        Rate[,i] <- rate.max * exp( -(xi[i] - x)^2 / (2*omega^2) )
    }
    
    # Computing the number of "+" alleles  per multilocus genotype
    # This code checks the even positions of the name of the genotypes to identify if the allele of each locus under selection is of type "1"
    # The first n.sel loci are the ones under selection, so the for loops stop at n.sel
    # The first for loop is for the first gamete (to the left hand side of the slash) and the second for loop is for the second gamete (to the right hand side of the slash)
    al <- array(0,dim=m)
    for (k in 1:m){
        # First gamete of the genotype
        for (i.pos in 1 : n.sel) {
            pos <- (2*i.pos)
            if(substr(colnames(init.par$meiosis_matrix)[k],pos,pos)=="1"){
                al[k]<-al[k]+1
            }
        }
        # Second gamete of the genotype
        for (i.pos in 1 : n.sel) {
            pos <- (nloci*2) + 1 + (2*i.pos)
            if(substr(colnames(init.par$meiosis_matrix)[k],pos,pos)=="1"){
                al[k]<-al[k]+1
            }
        }
    }
    
    # Assigning vital rates to multilocus genotypes as a function of the number of "+" alleles
    for (x in 1 : z) {
        for (t in 1 : T_max) {
            for(k in 1:m) {
                rate_OUT[k,,x,t] <- Rate[,al[k]+1]
            }
        }
    }
    dimnames(rate_OUT) <- list(genotype = dimnames(init.par$meiosis_matrix)$genotype,
                               deme = c(1:init.par$n),
                               age = c(1:init.par$z),
                               time = c(1:T_max))
    return(rate_OUT)
}

