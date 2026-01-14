library('SyncRNG')
seed = 123456
set.seed(seed)
s <- SyncRNG(seed = seed)
universal_rng <- SyncRNG(seed = seed)

syncrng.box.muller <- function(mu, sigma, n, seed = 0, rng = NULL) {
  rng <- if (is.null(rng)) universal_rng else rng
  two.pi <- 2 * pi
  ngen <- ceiling(n / 2)
  out <- replicate(2 * ngen, 0.0)
  for (i in 1:ngen) {
    u1 <- 0.0; u2 <- 0.0
    while (u1 == 0) { u1 <- rng$rand() }
    while (u2 == 0) { u2 <- rng$rand() }
    mag <- sigma * sqrt(-2.0 * log(u1))
    z0 <- mag * cos(two.pi * u2) + mu
    z1 <- mag * sin(two.pi * u2) + mu
    out[2*i - 1] = z0
    out[2*i] = z1
  }
  return(out[1:n])
}

check_reverse_offdiagonal <- function(mat) {
  regen <- FALSE
  for (i in 1:nrow(mat)) {
    for (j in 1:i) {
      if ((j != i) && (mat[i, j] != 0.0) && (mat[j, i] != 0.0)) {
        regen <- TRUE
        break
      }
    }
  }
  return(regen)
}

mat.generate.asw <- function(nvar, ar.value, dens, p.group,
                             lag.b, con.b,
                             ar.sd = 0.01,
                             lag.sd = 0.01,
                             con.sd = 0.01,
                             user_ar_indices = NULL,
                             laggroup.random = TRUE,
                             cntgroup.random = TRUE,
                             user_lag_indices = NULL,
                             user_con_indices = NULL) {
  repeat {
    Phi <- matrix(0, nrow = nvar, ncol = nvar)
    A <- matrix(0, nrow = nvar, ncol = nvar)
    
    n_lag_paths <- nvar * nvar - nvar
    n_con_paths <- nvar * (nvar * (nvar - 1))
    total_paths <- n_lag_paths + n_con_paths
    
    total_active_paths <- round(dens * total_paths)
    n_active_lag <- round(n_lag_paths / total_paths * total_active_paths)
    n_active_con <- total_active_paths - n_active_lag
    
    n_group_lag <- round(p.group * n_active_lag)
    n_group_con <- round(p.group * n_active_con)
    
    if (!is.null(user_ar_indices) && nrow(user_ar_indices) > 0) {
      for (i in seq_len(nrow(user_ar_indices))) {
        r <- user_ar_indices[i, 1]
        c <- user_ar_indices[i, 2]
        if (r == c && r <= nvar) {
          Phi[r, c] <- ar.value + syncrng.box.muller(0, ar.sd, 1, rng = s)
        }
      }
    }
    
    lag_indices <- which(Phi == 0, arr.ind = TRUE)
    lag_indices <- lag_indices[lag_indices[, 1] != lag_indices[, 2], , drop = FALSE]
    if (laggroup.random) {
      idx <- s$shuffle(1:nrow(lag_indices))[1:n_group_lag]
      selected_lag <- lag_indices[idx, , drop = FALSE]
    } else {
      selected_lag <- user_lag_indices
    }
    lag_betas <- syncrng.box.muller(lag.b, lag.sd, nrow(selected_lag), rng = s)
    for (i in seq_len(nrow(selected_lag))) {
      Phi[selected_lag[i, 1], selected_lag[i, 2]] <- lag_betas[i]
    }
    
    con_indices <- which(A == 0, arr.ind = TRUE)
    con_indices <- con_indices[con_indices[, 1] != con_indices[, 2], , drop = FALSE]
    if (cntgroup.random) {
      idx <- s$shuffle(1:nrow(con_indices))[1:n_group_con]
      selected_con <- con_indices[idx, , drop = FALSE]
    } else {
      selected_con <- user_con_indices
    }
    con_betas <- syncrng.box.muller(con.b, con.sd, nrow(selected_con), rng = s)
    for (i in seq_len(nrow(selected_con))) {
      A[selected_con[i, 1], selected_con[i, 2]] <- con_betas[i]
    }
    
    all_mat <- cbind(Phi, A)
    ind.pres <- which(all_mat != 0, arr.ind = TRUE)
    all_lvl <- matrix(NA, nrow = nvar, ncol = 2 * nvar)
    all_lvl[ind.pres] <- "grp"
    res <- list(sub1 = all_mat, lvl1 = all_lvl)
    break
  }
  return(res)
}

ts.generate.asw <- function (mat, lvl, t) {
  repeat {
    v <- ncol(mat)/2
    Phi <- mat[, 1:v]
    A <- mat[, (v+1):(2*v)]
    
    st <- t + 50
    noise <- matrix(rnorm(v * st, 0, 1), v)
    I <- diag(v)
    time <- matrix(0, nrow = v, ncol = st + 1)
    time1 <- matrix(0, nrow = v, ncol = st)
    
    for (i in 1:st) {
      time1[, i] <- solve(I - A) %*% (Phi %*% time[, i] + noise[, i])
      time[, i + 1] <- time1[, i]
    }
    
    time1 <- time1[, 51:(50 + t)]
    series <- t(time1)
    paths <- cbind(Phi, A)
    if (abs(max(series, na.rm = TRUE)) < 20 & abs(min(series, na.rm = TRUE)) > .01 & abs(min(series, na.rm = TRUE)) < 20) break
  }
  
  lvl[is.na(lvl) & paths != 0] <- "ind"
  
  list(
    series = series,
    paths = paths,
    nonoisepaths = paths,
    levels = lvl
  )
}
