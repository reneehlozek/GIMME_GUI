# GUI for Generating Time Series Matrices with Missing Data
# Combined with ASW simulation functions for correct levels and true output

# Required libraries
library(shiny)
library(ggplot2)
library(reshape2)
library(shinyBS)
library(shinyFiles)
library(SyncRNG)

# Initialize random number generators
seed <- 163422312
set.seed(seed)
s <- SyncRNG(seed = seed)
universal_rng <- SyncRNG(seed = seed)

#######################################################
#### 1. Simulation Functions (from working code) ######
#######################################################

syncrng.box.muller <- function(mu, sigma, n, seed = 0, rng = NULL) {
  # Code to ensure that we can have simultaneous random generation from python and/or R
  rng <- if (is.null(rng)) universal_rng else rng
  
  two.pi <- 2 * pi
  ngen <- ceiling(n / 2)
  out <- replicate(2 * ngen, 0.0)
  
  for (i in 1:ngen) {
    u1 <- 0.0
    u2 <- 0.0
    
    while (u1 == 0) { u1 <- rng$rand() }
    while (u2 == 0) { u2 <- rng$rand() }
    
    mag <- sigma * sqrt(-2.0 * log(u1))
    z0 <- mag * cos(two.pi * u2) + mu
    z1 <- mag * sin(two.pi * u2) + mu
    
    out[2 * i - 1] <- z0
    out[2 * i] <- z1
  }
  return(out[1:n])
}

# Code to check if we have dual paths (forward and back) and if so to regenerate
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
mat.generate.asw <- function(p.con, nvar, ar.val, dens, cnt.group, lag.group, con.b, lag.b,
                             ar_paths = NULL, lag_paths = NULL, con_paths = NULL,
                             laggroup.random = TRUE, cntgroup.random = TRUE) {
  repeat {
    A <- matrix(0, ncol = nvar, nrow = nvar)
    Phi <- matrix(0, ncol = nvar, nrow = nvar)
    
    cntpos.all <- (nrow(Phi) * ncol(Phi)) - nrow(Phi)
    lagpos.all <- (nrow(A) * ncol(A)) - nrow(A)
    
    cntpos.used <- dens * cntpos.all
    lagpos.used <- dens * lagpos.all
    
    cntg.all <- cnt.group * cntpos.used
    lagg.all <- lag.group * lagpos.used
    
    indices <- which(Phi == 0, arr.ind = TRUE)
    
    # Set the AR paths (diagonal of lagged matrix)
    if (!is.null(ar_paths) && nrow(ar_paths) > 0) {
      lagdiag.set <- c()
      for (i in 1:nrow(ar_paths)) {
        idx <- which(indices[, 1] == ar_paths[i, 1] & indices[, 2] == ar_paths[i, 2])
        if (length(idx) > 0) lagdiag.set <- c(lagdiag.set, idx)
      }
    } else {
      lagdiag.set <- which(indices[, 1] == indices[, 2])
    }
    laggrp.diag <- lagdiag.set
    
    # If you do not specify the lagged group paths specifically, then this will generate 
    # random group lagged paths    
    if (laggroup.random) {
      lagrow.col <- s$shuffle(c(1:(nvar * nvar)))[1:max(1, round(lagg.all))]
    } else {
      if (!is.null(lag_paths) && nrow(lag_paths) > 0) {
        lagrow.col <- c()
        for (i in 1:nrow(lag_paths)) {
          idx <- which(indices[, 1] == lag_paths[i, 1] & indices[, 2] == lag_paths[i, 2])
          if (length(idx) > 0) lagrow.col <- c(lagrow.col, idx)
        }
      } else {
        lagrow.col <- s$shuffle(c(1:(nvar * nvar)))[1:max(1, round(lagg.all))]
      }
    }
    grp.lag <- lagrow.col
   
# If you do not specify the contemporary paths specifically, then this will generate 
# random contemporary paths
    if (cntgroup.random) {
      cntrow.col <- s$shuffle(c(1:(nvar * nvar)))[1:max(1, round(cntg.all))]
    } else {
      if (!is.null(con_paths) && nrow(con_paths) > 0) {
        cntrow.col <- c()
        for (i in 1:nrow(con_paths)) {
          idx <- which(indices[, 1] == con_paths[i, 1] & indices[, 2] == con_paths[i, 2])
          if (length(idx) > 0) cntrow.col <- c(cntrow.col, idx)
        }
      } else {
        cntrow.col <- s$shuffle(c(1:(nvar * nvar)))[1:max(1, round(cntg.all))]
      }
    }
    grp.con <- cntrow.col
    
    if (length(grp.lag) > 0) Phi[indices[grp.lag, , drop = FALSE]] <- lag.b
    if (length(grp.con) > 0) A[indices[grp.con, , drop = FALSE]] <- con.b
    if (length(laggrp.diag) > 0) Phi[indices[laggrp.diag, , drop = FALSE]] <- ar.val
    
    break
  }
  
  all <- cbind(Phi, A)
  ind.pres <- which(all != 0, arr.ind = TRUE)
  level <- "grp"
  all.lvl <- matrix(NA, ncol = ncol(all), nrow = nrow(all))
  all.lvl[ind.pres] <- level
  
  res <- list(
    sub1 = all,
    lvl1 = all.lvl,
    Phi = Phi,
    A = A
  )
  return(res)
}
# Iterate over all the people and options 
ts.generate.asw <- function(mat, lvl, t, dens, cnt.group, lag.group, con.b, lag.b, p.con, noise_sd = 0.01) {
  repeat {
    repeat {
      v <- ncol(mat) / 2
      Phi <- mat[, 1:v]
      A <- mat[, (v + 1):(v * 2)]
      
      group_A_reverse <- check_reverse_offdiagonal(A)
      if (group_A_reverse) print('Warning: issue with group level contemp paths')
      
      indices.A <- which(A == 0, arr.ind = TRUE)
      indices.Phi <- which(Phi == 0, arr.ind = TRUE)
      
      group_phi_reverse <- check_reverse_offdiagonal(Phi)
      if (group_phi_reverse) print('Warning: issue with group level lagged paths')
      
      indices.A <- indices.A[which(indices.A[, 1] != indices.A[, 2]), , drop = FALSE]
      
      cntpos.all <- (nrow(Phi) * ncol(Phi)) - nrow(Phi)
      lagpos.all <- (nrow(A) * ncol(A)) - nrow(A)
      
      cntpos.used <- dens * cntpos.all
      lagpos.used <- dens * lagpos.all
      
      cntg.all <- cnt.group * cntpos.used
      lagg.all <- lag.group * lagpos.used
      
      cnt.all <- cntg.all + lagg.all
      
      vtest <- c(round(cnt.all / 2), (round(cnt.all) - round(cnt.all / 2)))
      randtest <- s$shuffle(vtest)[1:2]
      
      if (nrow(indices.A) > 0 && randtest[1] > 0) {
        atest <- 1:nrow(indices.A)
        n_to_select <- min(randtest[1], length(atest))
        if (n_to_select > 0) {
          row.col.A <- s$shuffle(atest)[1:n_to_select]
          A[indices.A[row.col.A, , drop = FALSE]] <- con.b
        }
      }
      
      if (nrow(indices.Phi) > 0 && randtest[2] > 0) {
        atest <- 1:nrow(indices.Phi)
        n_to_select <- min(randtest[2], length(atest))
        if (n_to_select > 0) {
          row.col.Phi <- s$shuffle(atest)[1:n_to_select]
          Phi[indices.Phi[row.col.Phi, , drop = FALSE]] <- lag.b
        }
      }
      
      regen <- check_reverse_offdiagonal(A)
      attempt <- 100 # Changed from 0 and attempt from 100
      while (regen && attempt < 100) {
        A <- mat[, (v + 1):(v * 2)]
        if (nrow(indices.A) > 0 && randtest[1] > 0) {
          atest <- 1:nrow(indices.A)
          randtest <- s$shuffle(vtest)[1:2]
          n_to_select <- min(randtest[1], length(atest))
          if (n_to_select > 0) {
            row.col.A <- s$shuffle(atest)[1:n_to_select]
            A[indices.A[row.col.A, , drop = FALSE]] <- con.b
          }
        }
        regen <- check_reverse_offdiagonal(A)
        attempt <- attempt + 1
      }
      
      regen <- check_reverse_offdiagonal(Phi)
      attempt <- 100 # Changed from 0 and attempt from 100
      while (regen && attempt < 100) {
        Phi <- mat[, 1:v]
        if (nrow(indices.Phi) > 0 && randtest[2] > 0) {
          atest <- 1:nrow(indices.Phi)
          randtest <- s$shuffle(vtest)[1:2]
          n_to_select <- min(randtest[2], length(atest))
          if (n_to_select > 0) {
            row.col.Phi <- s$shuffle(atest)[1:n_to_select]
            Phi[indices.Phi[row.col.Phi, , drop = FALSE]] <- lag.b
          }
        }
        regen <- check_reverse_offdiagonal(Phi)
        attempt <- attempt + 1
      }
      
      nonoisepaths <- cbind(Phi, A)
      
      noise.inds <- which(A != 0, arr.ind = TRUE)
      if (nrow(noise.inds) > 0) {
        A[noise.inds] <- A[noise.inds] + syncrng.box.muller(0, noise_sd, n = nrow(noise.inds))
      }
      
      noise.inds <- which(Phi != 0, arr.ind = TRUE)
      if (nrow(noise.inds) > 0) {
        Phi[noise.inds] <- Phi[noise.inds] + syncrng.box.muller(0, noise_sd, n = nrow(noise.inds))
      }
      
      break
    }
    
    st <- (t + 50)
    noise <- matrix(rnorm(v * st, 0, 1), v)
    I <- diag(v)
    time <- matrix(0, nrow = v, ncol = (st + 1))
    time1 <- matrix(0, nrow = v, ncol = st)
    
    for (i in 1:st) {
      time1[, i] <- solve(I - A) %*% (Phi %*% time[, i] + noise[, i])
      time[, i + 1] <- time1[, i]
    }
    
    time1 <- time1[, (51:(50 + t))]
    series <- t(time1)
    paths <- cbind(Phi, A)
    
    if (abs(max(series, na.rm = TRUE)) < 20 & abs(min(series, na.rm = TRUE)) > .01 &
        abs(min(series, na.rm = TRUE)) < 20) break
  }
  
  lvl[is.na(lvl) & paths != 0] <- "ind"
  
  list(
    series = series,
    paths = paths,
    nonoisepaths = nonoisepaths,
    levels = lvl,
    Phi = Phi,
    A = A
  )
}

# Code to drop certain variables (or whole days) to simulate observational conditions
add_missing_data <- function(series, drop_pct_row = 0.20, target_median_block = 2, drop_pct_random = 0.05) {
  n_rows <- nrow(series)
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
  series[rows_to_drop, ] <- NA
  
  non_missing_mask <- !is.na(series)
  available_indices <- which(non_missing_mask, arr.ind = TRUE)
  n_total_available <- nrow(available_indices)
  n_drop_random <- round(drop_pct_random * prod(dim(series)))
  
  if (n_total_available > 0 && n_drop_random > 0) {
    sampled_indices <- available_indices[
      sample(1:n_total_available, min(n_drop_random, n_total_available)),
      , drop = FALSE
    ]
    for (idx in seq_len(nrow(sampled_indices))) {
      series[sampled_indices[idx, 1], sampled_indices[idx, 2]] <- NA
    }
  }
  
  return(series)
}

#######################################################
#### 2. Shiny UI ######################################
#######################################################

ui <- fluidPage(
  titlePanel("Time Series Matrix Simulator"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n_people", "Number of People in Study:", value = 10, min = 1),
      numericInput("time_steps", "Number of Time Steps:", value = 200, min = 1),
      numericInput("n_vars", "Number of Variables (N):", value = 6, min = 2),
      bsTooltip("n_vars", "Contemporaneous total = N*(N-1), Lagged total = N*N.", "right", options = list(container = "body")),
      
      numericInput("network_density", "Total Network Density (%)", value = 20, min = 0, max = 100),
      bsTooltip("network_density", "Percent of total possible paths (contemp + lagged) that are active.", "right", options = list(container = "body")),
      
      numericInput("cnt_group", "Proportion of Contemp Paths that are Group-level (%)", value = 10, min = 0, max = 100),
      
      numericInput("lag_group", "Proportion of Lagged Paths that are Group-level (%)", value = 30, min = 0, max = 100),
      # lag_group not used!
      hr(),
      h3("Lagged Matrix (Φ)"),
      
      textAreaInput("ar_paths", "User-defined AR Paths (row,col pairs)", placeholder = "1,1; 2,2; 3,3"),
      numericInput("ar_coeff", "AR Coefficient:", value = 0.6, step = 0.1),
      
      radioButtons("group_model_lag", "Group-level Path Model (Lagged):",
                   choices = c("Random" = "random", "User-specified" = "confirm")),
      
      conditionalPanel(
        condition = "input.group_model_lag == 'confirm'",
        textAreaInput("group_indices_lag", "Group-level Lagged Paths (row,col pairs)", placeholder = "1,2; 3,4")
      ),
      
      numericInput("lag_beta", "Lagged Path Coefficient:", value = 0.45, step = 0.05),
      checkboxInput("lag_negative", "Include negative lagged coefficients", value = TRUE),
      
      hr(),
      h3("Contemporaneous Matrix (A)"),
      
      radioButtons("group_model_con", "Group-level Path Model (Contemporaneous):",
                   choices = c("Random" = "random", "User-specified" = "confirm")),
      
      conditionalPanel(
        condition = "input.group_model_con == 'confirm'",
        textAreaInput("group_indices_con", "Group-level Contemporaneous Paths (row,col pairs)", placeholder = "1,2; 2,3")
      ),
      
      numericInput("con_beta", "Contemporaneous Path Coefficient:", value = 0.35, step = 0.05),
      checkboxInput("con_negative", "Include negative contemp coefficients", value = TRUE),
      
      hr(),
      h3("Noise & Missing Data Settings"),
      numericInput("coeff_noise_sd", "Coefficient Noise SD:", value = 0.01, min = 0, step = 0.01),
      numericInput("drop_pct_row", "Percentage of days dropped (%):", value = 20, min = 0, max = 100),
      numericInput("target_median_block", "Median number of days dropped in a row:", value = 2, min = 1),
      numericInput("drop_pct_random", "Random variable dropout (%):", value = 5, min = 0, max = 100),
      
      hr(),
      h3("Output Settings"),
      # Number of repetitions (replaces single rep number) 
      numericInput("n_reps", "Number of Repetitions:", value = 1, min = 1),
      numericInput("random_seed", "Random Seed (for first repetition):", value = 163422312),

      
      actionButton("run_sim", "Run Simulation", class = "btn-primary"),
      br(), br(),
      downloadButton("download_all", "Download All Data (ZIP)"),
      br(), br(),
      shinyDirButton("folder", "Choose Output Folder", "Select a folder to save files"),
      actionButton("save_to_folder", "Save to Selected Folder", class = "btn-success"),
      verbatimTextOutput("folder_path")
    ),
    
    mainPanel(
      h4("Simulation Status"),
      verbatimTextOutput("sim_status"),
      
      hr(),
      uiOutput("rep_selector"),
      h4("Select Person to View"),
      uiOutput("person_selector"),
      
      tabsetPanel(
        tabPanel("Individual Matrices",
                 h4("Individual Lagged Matrix (Φ) - with noise"),
                 plotOutput("phi_indiv_plot"),
                 h4("Individual Contemporaneous Matrix (A) - with noise"),
                 plotOutput("a_indiv_plot"),
                 h4("Path Levels (grp/ind)"),
                 plotOutput("levels_plot")
        ),
        tabPanel("Group-level Matrices",
                 h4("Group-level Lagged Matrix (Φ)"),
                 plotOutput("phi_group_plot"),
                 h4("Group-level Contemporaneous Matrix (A)"),
                 plotOutput("a_group_plot")
        ),
        tabPanel("Time Series",
                 h4("Simulated Time Series (Selected Person)"),
                 plotOutput("ts_plot"),
                 h4("Time Series Overview (First 5 People)"),
                 plotOutput("ts_overview_plot", height = "600px")
        ),
        tabPanel("Output Preview",
                 h4("Data File Preview"),
                 tableOutput("data_preview"),
                 h4("True Matrix Preview"),
                 tableOutput("true_preview"),
                 h4("Levels Matrix Preview"),
                 tableOutput("levels_preview")
        )
      )
    )
  )
)

#######################################################
#### 3. Shiny Server ##################################
#######################################################

server <- function(input, output, session) {
  # Choose your home directory
  volumes <- c(Home = fs::path_home(), getVolumes()())
  shinyDirChoose(input, "folder", roots = volumes, session = session)
  
  selected_folder <- reactive({
    if (is.integer(input$folder)) return(NULL)
    parseDirPath(volumes, input$folder)
  })
  
  output$folder_path <- renderText({
    folder <- selected_folder()
    if (is.null(folder) || length(folder) == 0) {
      "No folder selected"
    } else {
      paste("Selected folder:", folder)
    }
  })
  
  parse_indices <- function(txt) {
    if (nzchar(txt)) {
      pairs <- unlist(strsplit(txt, ";"))
      result <- do.call(rbind, lapply(pairs, function(pair) {
        as.integer(unlist(strsplit(trimws(pair), ",")))
      }))
      if (is.null(dim(result))) result <- matrix(result, ncol = 2, byrow = TRUE)
      return(result)
    } else {
      return(matrix(nrow = 0, ncol = 2))
    }
  }
  
  # ── Store results for ALL repetitions ─────────────────────────────────────
  # sim_results is a list indexed by repetition number.
  # Each element has the same structure as before.
  sim_results <- reactiveValues(
    reps = list(),   # sim_results$reps[[rep_idx]] = list(group_mat, ..., individual_results, n_people)
    n_reps_done = 0
  )
  
  
  # ── Helper: run one full repetition ───────────────────────────────────────
  run_one_rep <- function(rep_idx, rep_seed, input) {
    
    set.seed(rep_seed)
    s <<- SyncRNG(seed = rep_seed)
    universal_rng <<- SyncRNG(seed = rep_seed)
    
    n_people <- input$n_people
    v        <- input$n_vars
    t        <- input$time_steps
    
    con_beta <- input$con_beta
    lag_beta <- input$lag_beta
    
    negcon <- if (input$con_negative) c(con_beta, -con_beta) else c(con_beta)
    neglag <- if (input$lag_negative) c(lag_beta, -lag_beta) else c(lag_beta)
    
    con.b <- s$shuffle(negcon)[1]
    lag.b <- s$shuffle(neglag)[1]
    
    ar_paths  <- parse_indices(input$ar_paths)
    lag_paths <- parse_indices(input$group_indices_lag)
    con_paths <- parse_indices(input$group_indices_con)
    
    res <- mat.generate.asw(
      p.con         = 0.50,
      nvar          = v,
      ar.val        = input$ar_coeff,
      cnt.group     = input$cnt_group / 100,
      lag.group     = input$lag_group / 100,
      dens          = input$network_density / 100,
      con.b         = con.b,
      lag.b         = lag.b,
      ar_paths      = if (nrow(ar_paths)  > 0) ar_paths  else NULL,
      lag_paths     = if (nrow(lag_paths) > 0) lag_paths else NULL,
      con_paths     = if (nrow(con_paths) > 0) con_paths else NULL,
      laggroup.random = (input$group_model_lag == "random"),
      cntgroup.random = (input$group_model_con == "random")
    )
    
    individual_results <- vector("list", n_people)
    
    for (a in 1:n_people) {
      out <- ts.generate.asw(
        mat       = res$sub1,
        lvl       = res$lvl1,
        t         = t,
        cnt.group = input$cnt_group / 100,
        lag.group = input$lag_group / 100,
        dens      = input$network_density / 100,
        con.b     = con.b,
        lag.b     = lag.b,
        p.con     = 0.50,
        noise_sd  = input$coeff_noise_sd
      )
      
      out$series <- round(out$series, digits = 5)
      out$series <- add_missing_data(
        out$series,
        drop_pct_row        = input$drop_pct_row / 100,
        target_median_block = input$target_median_block,
        drop_pct_random     = input$drop_pct_random / 100
      )
      colnames(out$series) <- paste0("V", 1:v)
      individual_results[[a]] <- out
    }
    
    list(
      group_mat          = res$sub1,
      group_lvl          = res$lvl1,
      group_Phi          = res$Phi,
      group_A            = res$A,
      individual_results = individual_results,
      n_people           = n_people,
      rep_seed           = rep_seed
    )
  }
  # ──────────────────────────────────────────────────────────────────────────
  
  # ── Main simulation: loop over all repetitions ────────────────────────────
  observeEvent(input$run_sim, {
    n_reps     <- input$n_reps
    base_seed  <- input$random_seed
    
    # Clear previous results
    sim_results$reps       <- list()
    sim_results$n_reps_done <- 0
    
    withProgress(message = 'Running simulations...', value = 0, {
      for (rep_idx in 1:n_reps) {
        
        # Each repetition gets a deterministic but distinct seed
        rep_seed <- base_seed + (rep_idx - 1) * 1000
        
        incProgress(
          amount  = 1 / n_reps,
          detail  = paste0("Repetition ", rep_idx, " of ", n_reps, "...")
        )
        
        sim_results$reps[[rep_idx]] <- run_one_rep(rep_idx, rep_seed, input)
      }
      sim_results$n_reps_done <- n_reps
    })
  })
  # ──────────────────────────────────────────────────────────────────────────
  
  # ── Subfolder name now includes repetition number ─────────────────────────
  get_subfolder_name <- function(rep_idx) {
    paste0("t_", input$time_steps, "_ar_", input$ar_coeff, "_rep_", rep_idx)
  }
  # ──────────────────────────────────────────────────────────────────────────
  
  # ── Save files for one repetition ─────────────────────────────────────────
  save_rep_files <- function(base_dir, rep_idx) {
    rep_data       <- sim_results$reps[[rep_idx]]
    subfolder_name <- get_subfolder_name(rep_idx)
    
    data_dir            <- file.path(base_dir, "data",            subfolder_name)
    true_dir            <- file.path(base_dir, "true",            subfolder_name)
    true_noisefree_dir  <- file.path(base_dir, "true_noisefree",  subfolder_name)
    levels_dir          <- file.path(base_dir, "levels",          subfolder_name)
    
    dir.create(data_dir,           showWarnings = FALSE, recursive = TRUE)
    dir.create(true_dir,           showWarnings = FALSE, recursive = TRUE)
    dir.create(true_noisefree_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(levels_dir,         showWarnings = FALSE, recursive = TRUE)
    
    for (a in 1:rep_data$n_people) {
      out <- rep_data$individual_results[[a]]
      write.csv(out$series,       file.path(data_dir,           paste0("ind_", a, ".csv")), row.names = FALSE)
      write.csv(out$paths,        file.path(true_dir,           paste0("ind_", a, ".csv")), row.names = FALSE)
      write.csv(out$nonoisepaths, file.path(true_noisefree_dir, paste0("ind_", a, ".csv")), row.names = FALSE)
      write.csv(out$levels,       file.path(levels_dir,         paste0("ind_", a, ".csv")), row.names = FALSE)
    }
  }
  # ──────────────────────────────────────────────────────────────────────────
  
  #  Save simulation files (all repetitions + parameters) 
  save_simulation_files <- function(base_dir) {
    for (rep_idx in seq_len(sim_results$n_reps_done)) {
      save_rep_files(base_dir, rep_idx)
    }
    
    params <- data.frame(
      Parameter = c("n_people", "time_steps", "n_vars", "network_density",
                    "cnt_group", "lag_group", "ar_coeff", "con_beta", "lag_beta",
                    "coeff_noise_sd", "drop_pct_row", "target_median_block",
                    "drop_pct_random", "n_reps", "base_random_seed","ar_paths", "group_indices_lag","group_indices_con"),
      Value = c(input$n_people, input$time_steps, input$n_vars, input$network_density,
                input$cnt_group, input$lag_group, input$ar_coeff, input$con_beta, input$lag_beta,
                input$coeff_noise_sd, input$drop_pct_row, input$target_median_block,
                input$drop_pct_random, input$n_reps, input$random_seed,input$ar_paths, input$group_indices_lag,input$group_indices_con)
    )
    write.csv(params, file.path(base_dir, "simulation_parameters.csv"), row.names = FALSE)
    
    return(base_dir)
  }
  
  # Display output
  output$sim_status <- renderText({
    if (sim_results$n_reps_done == 0) {
      "No simulation run yet. Click 'Run Simulation' to start."
    } else {
      paste0(
        "Simulation complete ", sim_results$n_reps_done, " repetition(s) run.\n",
        sim_results$reps[[1]]$n_people, " people | ",
        input$time_steps, " time steps | ",
        input$n_vars, " variables.\n",
        "Folders: data/t_<t>_ar_<ar>_rep_<r>/ind_i.csv"
      )
    }
  })
  
  # Select which iteration you want to display
  output$rep_selector <- renderUI({
    req(sim_results$n_reps_done > 0)
    selectInput("selected_rep", "Select Repetition:",
                choices  = 1:sim_results$n_reps_done,
                selected = 1)
  })
  # ####
  
  output$person_selector <- renderUI({
    req(sim_results$n_reps_done > 0)
    n_people <- sim_results$reps[[1]]$n_people
    selectInput("selected_person", "Select Person:",
                choices = 1:n_people, selected = 1)
  })
  
  # Display currently selected rep data 
  current_rep <- reactive({
    req(sim_results$n_reps_done > 0)
    req(input$selected_rep)
    sim_results$reps[[as.integer(input$selected_rep)]]
  })
  
  # write_gimme_file <- function(){
  #   
  #   
  # }
  
  draw_matrix_plot <- function(mat, title, lag_labels = FALSE, show_labels = TRUE) {
    n <- nrow(mat)
    
    if (ncol(mat) == 2 * n) {
      rownames(mat) <- paste0("V", 1:n)
      colnames(mat) <- c(paste0("V", 1:n, "_lag"), paste0("V", 1:n))
    } else {
      rownames(mat) <- paste0("V", 1:n)
      colnames(mat) <- if (lag_labels) paste0("V", 1:n, "_lag") else paste0("V", 1:n)
    }
    
    melted <- melt(mat)
    
    if (show_labels) {
      melted$label <- ifelse(!is.na(melted$value) & melted$value != 0,
                             round(as.numeric(as.character(melted$value)), 4), "")
    } else {
      melted$label <- ifelse(!is.na(melted$value), as.character(melted$value), "")
    }
    
    melted$Var1 <- factor(melted$Var1, levels = rev(paste0("V", 1:n)))
    
    if (is.numeric(melted$value)) {
      ggplot(melted, aes(Var2, Var1, fill = value)) +
        geom_tile(color = "black") +
        geom_text(aes(label = label), size = 3) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey90") +
        labs(title = title) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 0),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
        ) +
        scale_x_discrete(position = "top")
    } else {
      melted$value_factor <- factor(melted$value, levels = c("grp", "ind", NA))
      ggplot(melted, aes(Var2, Var1, fill = value_factor)) +
        geom_tile(color = "black") +
        geom_text(aes(label = label), size = 3) +
        scale_fill_manual(values = c("grp" = "steelblue", "ind" = "coral", "NA" = "white"),
                          na.value = "white", name = "Level") +
        labs(title = title) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 0),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
        ) +
        scale_x_discrete(position = "top")
    }
  }
  
  # Make a plot of current_rep() 
  output$phi_group_plot <- renderPlot({
    draw_matrix_plot(current_rep()$group_Phi, "Group-level Lagged Matrix (Φ)", lag_labels = TRUE)
  })
  
  output$a_group_plot <- renderPlot({
    draw_matrix_plot(current_rep()$group_A, "Group-level Contemporaneous Matrix (A)")
  })
  
  output$phi_indiv_plot <- renderPlot({
    req(input$selected_person)
    p   <- as.integer(input$selected_person)
    out <- current_rep()$individual_results[[p]]
    v   <- input$n_vars
    draw_matrix_plot(out$paths[, 1:v], paste("Person", p, "- Lagged Matrix (Φ)"), lag_labels = TRUE)
  })
  
  output$a_indiv_plot <- renderPlot({
    req(input$selected_person)
    p   <- as.integer(input$selected_person)
    out <- current_rep()$individual_results[[p]]
    v   <- input$n_vars
    draw_matrix_plot(out$paths[, (v + 1):(2 * v)], paste("Person", p, "- Contemporaneous Matrix (A)"))
  })
  
  output$levels_plot <- renderPlot({
    req(input$selected_person)
    p   <- as.integer(input$selected_person)
    out <- current_rep()$individual_results[[p]]
    draw_matrix_plot(out$levels, paste("Person", p, "- Path Levels"), show_labels = FALSE)
  })
  
  output$ts_plot <- renderPlot({
    req(input$selected_person)
    p      <- as.integer(input$selected_person)
    series <- current_rep()$individual_results[[p]]$series
    df     <- as.data.frame(series)
    df$Time <- 1:nrow(df)
    df_melt <- melt(df, id.vars = "Time")
    
    ggplot(df_melt, aes(x = Time, y = value, color = variable)) +
      geom_line(na.rm = TRUE) +
      labs(title = paste("Simulated Time Series - Person", p), y = "Value", x = "Time") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$ts_overview_plot <- renderPlot({
    n_show   <- min(5, current_rep()$n_people)
    all_data <- do.call(rbind, lapply(1:n_show, function(p) {
      df         <- as.data.frame(current_rep()$individual_results[[p]]$series)
      df$Time    <- 1:nrow(df)
      df$Person  <- paste("Person", p)
      df
    }))
    df_melt <- melt(all_data, id.vars = c("Time", "Person"))
    
    ggplot(df_melt, aes(x = Time, y = value, color = variable)) +
      geom_line(na.rm = TRUE) +
      facet_wrap(~Person, ncol = 1, scales = "free_y") +
      labs(title = "Time Series Overview", y = "Value", x = "Time") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  
  output$data_preview <- renderTable({
    req(input$selected_person)
    p <- as.integer(input$selected_person)
    head(current_rep()$individual_results[[p]]$series, 10)
  })
  
  output$true_preview <- renderTable({
    req(input$selected_person)
    p <- as.integer(input$selected_person)
    current_rep()$individual_results[[p]]$paths
  })
  
  output$levels_preview <- renderTable({
    req(input$selected_person)
    p <- as.integer(input$selected_person)
    current_rep()$individual_results[[p]]$levels
  })
  
  # Download all the iterations in a zip folder
  output$download_all <- downloadHandler(
    filename = function() {
      paste0("simulation_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      req(sim_results$n_reps_done > 0)
      
      temp_dir  <- tempdir()
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      base_dir  <- file.path(temp_dir, paste0("simulation_", timestamp))
      dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
      
      withProgress(message = 'Creating ZIP file...', value = 0, {
        save_simulation_files(base_dir)
        incProgress(0.8)
        old_wd <- getwd()
        setwd(temp_dir)
        zip(file, paste0("simulation_", timestamp))
        setwd(old_wd)
        incProgress(0.2)
      })
    }
  )
  
  observeEvent(input$save_to_folder, {
    req(sim_results$n_reps_done > 0)
    folder <- selected_folder()
    
    if (is.null(folder) || length(folder) == 0) {
      showNotification("Please select a folder first!", type = "error")
      return()
    }
    
    withProgress(message = 'Saving files...', value = 0, {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      base_dir  <- file.path(folder, paste0("simulation_", timestamp))
      save_simulation_files(base_dir)
      incProgress(1)
    })
    
    showNotification(paste("Files saved to:", base_dir), type = "message", duration = 10)
  })
}

shinyApp(ui = ui, server = server)
