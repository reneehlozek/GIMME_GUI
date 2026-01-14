#230222: AR=.6; Contemp = .35; lagged = .45, noise=.01
rm(list=ls())
library('SyncRNG')
seed=123456
set.seed(seed)
s <- SyncRNG(seed=seed)
universal_rng <- SyncRNG(seed=seed)


#######################################################
#### 1. load functions needed to simulate the data ####
#######################################################

syncrng.box.muller <- function(mu, sigma, n, seed=0, rng=NULL) {
  "Generate random samples from a normal distribution using the Box-Muller transform.
  
  Args:
     mu: The mean of the normal distribution.
     sigma: The standard deviation of the normal distribution.
     n: The number of random samples to generate.
     seed: An optional seed value for the random number generator (default is 0).
     rng: An optional instance of SyncRNG. If NULL, a new instance will be created using the seed.
  
   Returns:
     A vector of n random samples from the specified normal distribution."
  
  rng <- if (is.null(rng)) universal_rng else rng
  
  two.pi <- 2 * pi
  ngen <- ceiling(n / 2)
  out <- replicate(2 * ngen, 0.0)
  
  for (i in 1:ngen) {
    u1 <- 0.0
    u2 <- 0.0
    
    while (u1 == 0) { u1 <- rng$rand(); }
    while (u2 == 0) { u2 <- rng$rand(); }
    
    mag <- sigma * sqrt(-2.0 * log(u1))
    z0 <- mag * cos(two.pi * u2) + mu
    z1 <- mag * sin(two.pi * u2) + mu
    
    out[2*i - 1] = z0;
    out[2*i] = z1;
  }
  return(out[1:n]);
}

check_reverse_offdiagonal <- function(mat){
  regen <- FALSE
  for (i in 1:nrow(mat)) {
    for (j in 1:i) {
      if ((j!=i) && (mat[i, j] != 0.0) && (mat[j, i] != 0.0)) {
        #cat('i=', i, ', j=', j, ', matij=', mat[i, j], ', matji=', mat[j, i], ', we must regenerate \n')
        regen <- TRUE
        break
      }
    }
  }
  return(regen)
  
}

# #check_symmetric_indices <- function(mat) {
#   # Get the non-zero indices of the matrix using arr.ind = TRUE
#   non_zero_indices <- which(mat != 0, arr.ind = TRUE)
#   new_nonz_indices <- non_zero_indices
#   # Initialize a vector to store results
#   results <- character()
#   # Check for each non-zero index if its corresponding off-diagonal symmetric index is in the list
#   countnonzero=0
#   for (i in seq(1, nrow(non_zero_indices), by = 2)) {
#     row_index <- non_zero_indices[i, "row"]
#     col_index <- non_zero_indices[i, "col"]
#     
#     # Check if the symmetric counterpart is off-diagonal and in the list
#     symmetric_index <- c(col_index, row_index)
#     if (row_index != col_index && any(apply(non_zero_indices, 1, function(x) all(x == symmetric_index)))) {
#       countnonzero <-countnonzero+1
#       indexarr <- convert_to_arr_ind_format(col_index,row_index,c(6,12))
#       mat[indexarr]<-0
#     }
#   }
#   
#   return(list(mat=mat,nonzeroels=countnonzero))
# }

#Group Relations for Set 1: AR: V2, V4, V6; Contemporary Connections: V5 -> V3, V6 -> V4; Lagged Connections: V1lagged -> V6, V4lagged -> V5

mat.generate.asw <- function(p.con,nvar,ar.val,dens,cnt.group,lag.group,con.b,lag.b){
  repeat{
    A <- matrix(0, ncol = nvar, nrow = nvar, )        ###create null contemporaneous matrix
    Phi <- matrix(0, ncol = nvar, nrow = nvar)        ###create null lagged matrix (remember -> ARs in diagonal)
    cntpos.all <- (nrow(Phi)*ncol(Phi)) - nrow(Phi)   ###list all available (non-ar) spots in contem matrix
    lagpos.all <- (nrow(A)*ncol(A)) - nrow(A)         ### list all available spots in lagged matrix (excluding AR)
    cnt.all <-cntpos.all+lagpos.all
    cntpos.used <-dens*cntpos.all
    lagpos.used <-dens*lagpos.all
    pos.used <- cntpos.used + lagpos.used         ## number of non-zero paths you want (to match GUI)
    
    #cnt.group <- 0.1
    #lag.group<-0.3
    cntg.all <- cnt.group*cntpos.used                   ###number of contemporaneous group paths given density and proportion of group paths
    lagg.all <- lag.group*lagpos.used
    indices <- which(Phi == 0, arr.ind = TRUE)  ###create data frame that lists all combinations of columns/rows in matrices. 
                                                ## Run "indices <- indices[which(indices[,1] != indices[,2]), ]" to kick out diagonal
    
    
    #  ## LAGGED MATRIX SPECIFICATION
    
    lagdiag.set<-c(8, 22, 36)                            ###select AR paths for the lagged matrix
    laggrp.diag <- lagdiag.set                           ###name AR paths
    
    laggroup.random <- FALSE ## specify if you want to assign your own group level or not
    
    ## CHANGE THIS IN THE FUNCTION CALL
    if(laggroup.random){
      lagrow.col = s$shuffle(c(1:(nvar*nvar+1)))[1:round(lagg.all)]
    } else {
  
      lagrow.col<-c(17, 24, 31)
    }
    
    
    grp.lag <- lagrow.col   ###name lagged group paths
    print('lagged group')
    print(grp.lag)
    
    
    ## CONTEMP MATRIX SPECIFICATION
    
    cntgroup.random <- FALSE                              ## specify if you want to assign your own group level or not
    if(cntgroup.random){
      cntrow.col = s$shuffle(c(1:(nvar*nvar+1)))[1:round(cntg.all)]
      #s$shuffle(c((nvar*nvar+1):(2*nvar*nvar)))[1:round(cntg.all)]
    } else {
      
      cntrow.col<-c(28)
    }
    
    grp.con <- cntrow.col                             ###name contemporaneous group paths
    print('contemp group')
    print(grp.con)
  
    Phi[indices[grp.lag,]] <- lag.b                   ###set lagged paths (Identified below)
    A[indices[grp.con,]] <- con.b                     ###set contemporaneous paths (Identified below)
    Phi[indices[laggrp.diag,]] <- ar.val                     ###set AR paths (Identified below)
    break      
    
  }
  
  all <- cbind(Phi, A)                                          ###combine contemporaneous and lagged matrix 
  ind.pres <- which(all != 0, arr.ind = T)                      ###identify the group paths 
  level <- "grp"
  all.lvl <- matrix(NA, ncol = ncol(all), nrow = nrow(all))
  all.lvl[ind.pres] <- level
  
  all_sub1 <- all
  all_lvl1 <- all.lvl
  
  res <- list(sub1 = all_sub1, 
              lvl1 = all_lvl1)
  return(res)
}

#This function generates group-level relations, adds individual-level relations to the matrix, 
#adds noise, and simulates the time series (unnecessary for our purposes but helpful for checks!)
ts.generate.asw <- function (mat, lvl, t,dens,cnt.group,lag.group,con.b,lag.b,p.con) {
  repeat {
    repeat{
      v <- ncol(mat)/2                                                                              ###calculate number of variables in matrix
      Phi <- mat[, 1:v]                                                                             ###pull out Phi
      A <- mat[, (v+1):(v*2)]                                                                       ###pull out A
      A_ind <- matrix(0, ncol=v, nrow=v)                                                            ###set A indices to zero
      group_A_reverse <-check_reverse_offdiagonal(A)
      if(group_A_reverse) {print('something is wrong with your group level contem paths')}
      indices.A <- which(A==0,arr.ind=T)                                                            ###finds indices of zero elements in matrix
      indices.Phi <- which(Phi==0,arr.ind=T)
      group_phi_reverse <-check_reverse_offdiagonal(Phi)
      if(group_phi_reverse) {print('something is wrong with your group level lagged paths')}
      
      indices.A <- indices.A[which(indices.A[,1]!=indices.A[,2]),]                                  ###kick out diagonal (from A matrix only)
      pos.all<-(v*v-v)*2      ###determine the number of individual paths to add
      
      cntpos.all <- (nrow(Phi)*ncol(Phi)) - nrow(Phi)  
      lagpos.all <- (nrow(A)*ncol(A)) - nrow(A) 
      
      cntpos.used <-dens*cntpos.all
      lagpos.used <-dens*lagpos.all
      
      cntg.all <- cnt.group*cntpos.used                   ###number of contemporaneous group paths given density and proportion of group paths
      lagg.all <- lag.group*lagpos.used
      
      
      cnt.all <- cntg.all+lagg.all 
      

      vtest <-c(round(cnt.all/2),(round(cnt.all)-round(cnt.all/2)))
      randtest <- s$shuffle(vtest)[1:2]
      # print(randtest)
      # print('randtest')
      # print(rand)
      # print('rand')
      
      # row.tol.A      <- sample(1:nrow(indices.A), rand[1], replace = F)
      # row.col.Phi      <- sample(1:nrow(indices.Phi), rand[2], replace = F)
      # print(row.col.A)
      # print(row.col.Phi)
      # print('--')
      atest <-1:nrow(indices.A)
      row.col.A <-s$shuffle(atest)[1:randtest[1]]
      atest <-1:nrow(indices.Phi)
      row.col.Phi   <- s$shuffle(atest)[1:randtest[2]]
      
      # print(row.col.A)
      # print(row.col.Phi)
      # print('hi')
      # print(huh)
      Phi[indices.Phi[row.col.Phi,]] <- lag.b                            
      ###set betas for lagged and contemporaneous paths
      A[indices.A[row.col.A,]]     <- con.b
      
      
      regen <-check_reverse_offdiagonal(A)
      
     
      while (regen==TRUE){
        
        A <- mat[, (v+1):(v*2)] # reset A to the group level
        atest <-1:nrow(indices.A)
        randtest <- s$shuffle(vtest)[1:2]
        row.col.A <-s$shuffle(atest)[1:randtest[1]]
        con.b <- s$shuffle(contest)[1]
        A[indices.A[row.col.A,]]     <- con.b
        regen <-check_reverse_offdiagonal(A)
        #print('we need to regen contemp')
        #print(con.b)
        }
      
      
      regen <-check_reverse_offdiagonal(Phi)
      while (regen==TRUE){
        
        Phi <- mat[, 1:v]
        atest <-1:nrow(indices.Phi)
        randtest <- s$shuffle(vtest)[1:2]
        row.col.Phi   <- s$shuffle(atest)[1:randtest[2]]
        ## RH added these again 
        
        lag.b <- s$shuffle(lagtest)[1]
        Phi[indices.Phi[row.col.Phi,]] <- lag.b 
        regen <-check_reverse_offdiagonal(Phi)
        #print('we need to regen lag')
        #print(lag.b)
        }
      
      nonoisepaths  <- cbind(Phi, A) ### bind Phi and A before adding noise to them to check
      noise.inds      <- which(A != 0, arr.ind = TRUE)                                            
      ####add noise to A betas, SD =.1
      A[noise.inds]   <- A[noise.inds] + syncrng.box.muller(0, 0.01, n=nrow(noise.inds))
      noise.inds      <- which(Phi != 0, arr.ind = TRUE)                                          
      ###add noise to Phi betas, SD =.1
      Phi[noise.inds]   <- Phi[noise.inds] + syncrng.box.muller(0, 0.01, n=nrow(noise.inds))
      break
    }
    
    st <- (t+50)  
    ###This chunk is from Alex's code and is used to simulate the time series 
    noise <- matrix(rnorm(v*st,0,1),v) #
    I     <- diag(v) # identity matrix
    time  <- matrix(0,nrow=v, ncol=(st+1))
    time1 <- matrix(0,nrow=v, ncol=st)
    
    for (i in 1:st){                                                                                ###simulate data points for each time step
      time1[,i]  <- solve(I-A)%*%(Phi%*%time[,i] + noise[,i])
      time[,i+1] <- time1[,i]
    }               
    time1  <- time1[,(51:(50+t))]                                                                   
    series <- t(time1)
    paths  <- cbind(Phi, A)
    if (abs(max(series, na.rm = TRUE)) < 20 & abs(min(series, na.rm = TRUE)) > .01 
        & abs(min(series, na.rm = TRUE)) < 20) break
  }
  
  lvl[is.na(lvl) & paths != 0] <- "ind"
  
  list   <- list("series"  = series,
                 "paths"   = paths,
                 "nonoisepaths" = nonoisepaths,
                 "levels"  = lvl)
  return(list)
}


#######################################################
#### 2. Simulate Data# ################################
#######################################################

# enter simulation parameters
v             <- c(6) # Number of variables
n             <- c(1) # number of individuals
t             <- c(200) # Number of time points
rep           <- seq(1) # replications per condition 
ar            <-c(.6) # ar paths to try
conditions    <- expand.grid(t, n, v, ar,rep)
all           <- rbind(conditions)
colnames(all) <- c("t", "n", "v", "ar", "rep")
cntbeta <- 0.35
lagbeta <-0.45
negcon <- c(cntbeta, -1*cntbeta)                                     #to create the negative numbers
neglag <- c(lagbeta, -1*lagbeta)                                     #to create the negative numbers

# This loop generates a folder name for each condition and replication
for (i in 1:nrow(all)){
  all$folder[i] <- paste("t",all$t[i],
                         "ar",all$ar[i],"rep",all$rep[i],sep="_")
}

rownames(all) <- NULL
colnames(all) <- c("t","n","v","ar","rep","folder")
all$t         <- as.numeric(as.character(all$t))
all$n         <- as.numeric(as.character(all$n))
all$ar         <- as.numeric(as.character(all$ar))


# name directories to place simulated data in (Change)
dir.create('/Users/reneehlozek/Dropbox/CIFAR_Sims/14012026_Missing')
data.path <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/14012026_Missing/data'
true.path <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/14012026_Missing/true'
true_noisefree.path <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/14012026_Missing/true_noisefree'
level.path <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/14012026_Missing/levels'

dir.create(data.path)
dir.create(true.path)
dir.create(true_noisefree.path)
dir.create(level.path)

# Creates actual folders for each iteration
folders <- all$folder
for (i in 1:nrow(all)){
  data <- file.path(data.path, folders[i])
  dir.create(data)
  true <- file.path(true.path, folders[i])
  dir.create(true)
  true_noisefree <- file.path(true_noisefree.path, folders[i])
  dir.create(true_noisefree)
  level <- file.path(level.path, folders[i])
  dir.create(level)
}

# Does the simulations
for (i in 1:nrow(all)){
  if (!length(list.files(file.path(data.path,all$folder[i]))) %in% c(50)) {
    
    contest<-negcon
    lagtest<-neglag
    con.b <- s$shuffle(contest)[1]
    lag.b <- s$shuffle(lagtest)[1]   # check that this is working
    #print('inside simulation', i)
    #print(con.b)
    #print(lag.b)
    
    # generate group matrix for each simulated data set
    
    # 4a = p.con
    # 4 = dens
    # 5bi1 = p.group
  
    res <- mat.generate.asw(p.con = .50, 
                            nvar = all$v[i], ar.val=all$ar[i],
                            cnt.group=0.1,lag.group=0.3,
                            dens = .20,
                            con.b = con.b, lag.b = lag.b)
    
    #for each individual,generate matrix and ts
    for (a in 1:all$n[i]){
      
      out <- ts.generate.asw(mat = res$sub1,
                             lvl = res$lvl1,
                             t   = all$t[i],
                             cnt.group=0.1,lag.group=0.3,
                             dens = .20,
                             con.b = con.b, lag.b = lag.b,
                             p.con = .50)
      
      out$series <- round(out$series,digits=5)
      
      ### Code for missing data
      
      
      drop_pct_row <-  0.30                #Natasha added missingness (row-level)
      target_median_block <- 2            #days in a row
      
      n_rows <- nrow(out$series)
      n_total_rows_to_drop <- round(drop_pct_row * n_rows)
      
      rows_to_drop <- c()
      counter <- 0
      
      while (counter < n_total_rows_to_drop) {
        block_size <- max(1, round(rnorm(1, mean = target_median_block, sd = 1)))
        available_rows <- setdiff(1:n_rows, rows_to_drop)
        if (length(available_rows) == 0) break
        
        start_row <- sample(available_rows, 1)
        block_rows <- start_row:(start_row + block_size - 1)
        block_rows <- block_rows[block_rows <= n_rows]
        
        rows_to_drop <- c(rows_to_drop, block_rows)
        counter <- length(unique(rows_to_drop))
      }
      
      rows_to_drop <- unique(rows_to_drop)
      out$series[rows_to_drop, ] <- NA
      
      drop_pct_random <-  0.05   #Add random non-overlapping missingness
      non_missing_mask <- !is.na(out$series)
      available_indices <- which(non_missing_mask, arr.ind = TRUE)
      n_total_available <- nrow(available_indices)
      n_drop_random <- round(drop_pct_random * prod(dim(out$series)))
      
      if (n_total_available > 0) {
        sampled_indices <- available_indices[
          sample(1:n_total_available, min(n_drop_random, n_total_available)),
          , drop = FALSE
        ]
        for (idx in seq_len(nrow(sampled_indices))) {
          out$series[sampled_indices[idx, 1], sampled_indices[idx, 2]] <- NA
        }
      }
      
      
      
      write.csv(out$series,
                file.path(data.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$paths,
                file.path(true.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$nonoisepaths,
                file.path(true_noisefree.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$levels,
                file.path(level.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
    }
    
  }
}

