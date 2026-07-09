#NC to rerun, making AR=.6 and make contemp/lagged=.45, and noise=.05
#230328 Update: AR=.5; Contemp = .3; lagged=.4, and noise=.01

rm(list=ls())

#######################################################
#### 1. load functions needed to simulate the data ####
#######################################################

#Group Relations for Set 1: AR: V2, V4, V6; Contemporary Connections: V5 -> V3, V6 -> V4; Lagged Connections: V1lagged -> V6, V4lagged -> V5

mat.generate.asw <- function(p.con,nvar,AR,dens,p.group,con.b,lag.b){
  repeat{
    A <- matrix(0, ncol = nvar, nrow = nvar, )        ###create null contemporaneous matrix
    Phi <- matrix(0, ncol = nvar, nrow = nvar)        ###create null lagged matrix (remember -> ARs in diagonal)
    pos.all <- 2*(nrow(A)*ncol(A))                    ###list all available (non-ar) spots (You could add - 2*nrow(A) to remove AR spots (see AW code))
    cnt.all <- dens*p.group*pos.all                   ###number of group paths given density and proportion of group paths
    indices <- which(Phi == 0, arr.ind = TRUE)        ###create data frame that lists all combinations of columns/rows in matrices. Run "indices <- indices[which(indices[,1] != indices[,2]), ]" to kick out diagonal
    row.col<-c(17, 24, 31, 28)                        ###select group-level associations
    diag.set<-c(8, 22, 36)                            ###select AR paths 
    n.p.1<-row.col[1:round(p.con*length(row.col))]    ###identify contemporaneous group paths 
    n.p.2<-row.col[(length(n.p.1)+1):length(row.col)] ###identify lagged group paths 
    grp.con <- n.p.1                                  ###name contemporaneous group paths
    grp.lag <- n.p.2                                  ###name lagged group paths
    grp.diag <- diag.set                              ###name AR paths
    Phi[indices[grp.lag,]] <- lag.b                   ###set lagged paths (Identified below)
    A[indices[grp.con,]] <- con.b                     ###set contemporaneous paths (Identified below)
    Phi[indices[grp.diag,]] <- AR                     ###set AR paths (Identified below)
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
ts.generate.asw <- function (mat, lvl, t,dens,p.group,con.b,lag.b,p.con) {
  repeat {
    
    repeat{
      v <- ncol(mat)/2                                                                              ###calculate number of variables in matrix
      Phi <- mat[, 1:v]                                                                             ###pull out Phi
      A <- mat[, (v+1):(v*2)]                                                                       ###pull out A
      A_ind <- matrix(0, ncol=v, nrow=v)                                                            ###set A indices to zero
      indices.A <- which(A==0,arr.ind=T)                                                            ###finds indices of zero elements in matrix
      indices.Phi <- which(Phi==0,arr.ind=T)
      indices.A <- indices.A[which(indices.A[,1]!=indices.A[,2]),]                                  ###kick out diagonal (from A matrix only)
      pos.all<-(v*v-v)*2                                                                            ###determine the number of individual paths to add
      cnt.all <- dens*p.group*pos.all  
      rand<-sample(c(round(cnt.all/2),(round(cnt.all)-round(cnt.all/2))),size=2,replace=FALSE)      ###randomly select row/col combinations for individual-level paths
      row.col.A      <- sample(1:nrow(indices.A), rand[1], replace = F) 
      row.col.Phi      <- sample(1:nrow(indices.Phi), rand[2], replace = F) 
      Phi[indices.Phi[row.col.Phi,]] <- lag.b                                                       ###set betas for lagged and contemporaneous paths
      A[indices.A[row.col.A,]]     <- con.b
      noise.inds      <- which(A != 0, arr.ind = TRUE)                                              ####add noise to A betas, SD =.01
      A[noise.inds]   <- A[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .01)
      noise.inds      <- which(Phi != 0, arr.ind = TRUE)                                            ###add noise to Phi betas, SD =.01
      Phi[noise.inds] <- Phi[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .01)
      break
    }
    
    st <- (t+50)                                                                                    ###This chunk is from Alex's code and is used to simulate the time series 
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
                 "levels"  = lvl)
  return(list)
}


#######################################################
#### 2. Simulate Data# ################################
#######################################################

# enter simulation parameters
v             <- c(6) # Number of variables
n             <- c(25) # number of individuals
t             <- c(200) # Number of time points
rep           <- seq(1) # replications per condition 
ar            <-c(.5) # ar paths to try
conditions    <- expand.grid(t, n, v, ar,rep)
all           <- rbind(conditions)
colnames(all) <- c("t", "n", "v", "ar", "rep")
negcon <- c(.3, -.3)                                     #to create the negative numbers
neglag <- c(.4, -.4)                                     #to create the negative numbers

# This loop generates a folder name for each condition and replication
for (i in 1:nrow(all)) {
  all$folder[i] <- paste("t", all$t[i], "ar", all$ar[i], "rep", all$rep[i], sep = "_")
}

rownames(all) <- NULL
colnames(all) <- c("t", "n", "v", "ar", "rep", "folder")
all$t  <- as.numeric(as.character(all$t))
all$n  <- as.numeric(as.character(all$n))
all$ar <- as.numeric(as.character(all$ar))

# Create parent directories (change if needed)
dir.create('~/Dropbox (Personal)/Research Projects/CIFAR/simset1missing25')
data.path  <- '~/Dropbox (Personal)/Research Projects/CIFAR/simset1missing25/data'
true.path  <- '~/Dropbox (Personal)/Research Projects/CIFAR/simset1missing25/true'
level.path <- '~/Dropbox (Personal)/Research Projects/CIFAR/simset1missing25/levels'

dir.create(data.path)
dir.create(true.path)
dir.create(level.path)

# Create subfolders for each condition
folders <- all$folder
for (i in 1:nrow(all)) {
  dir.create(file.path(data.path,  folders[i]))
  dir.create(file.path(true.path,  folders[i]))
  dir.create(file.path(level.path, folders[i]))
}

# Run simulations
for (i in 1:nrow(all)) {
  if (!length(list.files(file.path(data.path, all$folder[i]))) %in% c(50)) {
    
    # Generate matrix
    res <- mat.generate.asw(
      p.con = 0.50,
      nvar = all$v[i],
      AR = all$ar[i],
      p.group = 0.50,
      dens = 0.20,
      con.b = sample(negcon, 1),
      lag.b = sample(neglag, 1)
    )
    
    for (a in 1:all$n[i]) {
      
      out <- ts.generate.asw(
        mat  = res$sub1,
        lvl  = res$lvl1,
        t    = all$t[i],
        p.group = 0.50,
        dens = 0.20,
        con.b = sample(negcon, 1),
        lag.b = sample(neglag, 1),
        p.con = 0.50
      )
      
      out$series <- round(out$series, digits = 5)
      
      drop_pct_row <- 0.70                #Natasha added missingness (row-level)
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
      
      drop_pct_random <- 0.05   #Add random non-overlapping missingness
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

      # Save simulated files
      write.csv(out$series,
                file.path(data.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$paths,
                file.path(true.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$levels,
                file.path(level.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
    }
  }
}

