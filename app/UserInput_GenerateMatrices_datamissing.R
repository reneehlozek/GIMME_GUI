# GUI for Generating Time Series Matrices with Missing Data, based on Chaku R simulations where observations are dropped
# Required libraries
library(shiny)
library(ggplot2)
library(reshape2)
library(shinyBS)

droprows<- function(series, drop_pct_row=25, target_median_block=1,drop_pct_random=5){
  
  # checking for fractions
  if (drop_pct_row>1){
    drop_pct_row<-drop_pct_row/100
  }
  if (drop_pct_random>1){
    drop_pct_random<-drop_pct_random/100
  }
  
  
    
  #drop_pct_row <- 0.30                #Natasha added missingness (row-level)
  #target_median_block <- 2            #days in a row
  
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
  out$series[rows_to_drop, ] <- NA
#  drop_pct_random <- 0.05   #Add random non-overlapping missingness
  non_missing_mask <- !is.na(series)
  available_indices <- which(non_missing_mask, arr.ind = TRUE)
  n_total_available <- nrow(available_indices)
  n_drop_random <- round(drop_pct_random * prod(dim(out$series)))
  
  if (n_total_available > 0) {
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

#syncrng.box.muller <- function(mu, sigma, n, seed=0, rng=NULL)
# --- UI ---
ui <- fluidPage(
  titlePanel("Time Series Matrix Simulator"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n_people", "Number of People in Study:", value = 10, min = 1),
      numericInput("time_steps", "Number of Time Steps:", value = 100, min = 1),
      numericInput("n_vars", "Number of Variables (N):", value = 6, min = 2),
      bsTooltip("n_vars", "Contemporaneous total = N*(N-1), Lagged total = N*N.", "right", options = list(container = "body")),
      
      numericInput("network_density", "Total Network Density (%)", value = 20, min = 0, max = 100),
      bsTooltip("network_density", "Percent of total possible paths (contemp + lagged) that are active.", "right", options = list(container = "body")),
      
      numericInput("percent_contemp", "Percentage of Active Paths that are Contemporaneous (%)", value = 40, min = 0, max = 100),
      bsTooltip("percent_contemp", "Split of active paths: % contemporaneous vs % lagged.", "right", options = list(container = "body")),
      
      hr(),
      h3("Lagged Matrix (Φ)"),
      
      textAreaInput("ar_paths", "User-defined AR Paths (row,col pairs)", placeholder = "1,1; 2,2"),
      numericInput("ar_coeff", "AR Coefficient:", value = 0.5),
      
      radioButtons("group_model_lag", "Group-level Path Model (Lagged):",
                   choices = c("Random" = "random", "User-specified" = "confirm")),
      
      conditionalPanel(
        condition = "input.group_model_lag == 'confirm'",
        textAreaInput("group_indices_lag", "Group-level Lagged Paths (row,col pairs)", placeholder = "1,2; 3,4")
      ),
      
      numericInput("group_lag_beta", "Group-level Lagged Beta:", value = 0.3),
      numericInput("group_lag_sd", "Group-level Lagged Beta SD:", value = 0.01),
      
      numericInput("indiv_lag_beta", "Individual-level Lagged Beta:", value = 0.3),
      numericInput("indiv_lag_sd", "Individual-level Lagged Beta SD:", value = 0.01),
      
      hr(),
      h3("Contemporaneous Matrix (A)"),
      
      radioButtons("group_model_con", "Group-level Path Model (Contemporaneous):",
                   choices = c("Random" = "random", "User-specified" = "confirm")),
      
      conditionalPanel(
        condition = "input.group_model_con == 'confirm'",
        textAreaInput("group_indices_con", "Group-level Contemporaneous Paths (row,col pairs)", placeholder = "1,2; 2,3")
      ),
      
      numericInput("group_con_beta", "Group-level Contemporaneous Beta:", value = 0.3),
      numericInput("group_con_sd", "Group-level Contemporaneous Beta SD:", value = 0.01),
      
      numericInput("indiv_con_beta", "Individual-level Contemporaneous Beta:", value = 0.3),
      numericInput("indiv_con_sd", "Individual-level Contemporaneous Beta SD:", value = 0.01),
      
      hr(),
      h3("Timeseries Noise Settings"),
      numericInput("overall_noise", "Overall Timeseries Noise SD:", value = 0.1, min = 0),
      numericInput("drop_pct_row", "Percentage of days dropped (%):", value = 25, min = 0),
      numericInput("target_median_block", "Median number of days dropped in a row:", value = 2, min = 0), 
      numericInput("drop_pct_random", "Percentage of observations where random variables are dropped (%):", value = 5, min = 0), 
      
      actionButton("run_sim", "Run Simulation")
    ),
    
    mainPanel(
      h4("Lagged Matrix Output"),
      plotOutput("phi_plot"),
      h4("Contemporaneous Matrix Output"),
      plotOutput("a_plot"),
      h4("Simulated Time Series Output"),
      plotOutput("ts_plot")
    )
  )
)

# --- SERVER ---
server <- function(input, output) {
  
  parse_indices <- function(txt) {
    if (nzchar(txt)) {
      pairs <- unlist(strsplit(txt, ";"))
      do.call(rbind, lapply(pairs, function(pair) as.integer(unlist(strsplit(trimws(pair), ",")))))
    } else {
      matrix(nrow = 0, ncol = 2)
    }
  }
  
  generate_matrices <- eventReactive(input$run_sim, {
    v <- input$n_vars
    A <- matrix(0, v, v)
    Phi <- matrix(0, v, v)
    
    # Total possible paths
    total_contemp_paths <- v * (v - 1)
    total_lagged_paths <- v * v
    total_possible_paths <- total_contemp_paths + total_lagged_paths
    
    # Number of active paths based on network density
    n_active_paths <- round(input$network_density / 100 * total_possible_paths)
    
    # Split active paths into contemporaneous and lagged
    n_contemp_active <- round(input$percent_contemp / 100 * n_active_paths)
    n_lagged_active <- n_active_paths - n_contemp_active
    
    # -- AR paths
    ar_coords <- parse_indices(input$ar_paths)
    for (i in seq_len(nrow(ar_coords))) {
      Phi[ar_coords[i, 1], ar_coords[i, 2]] <- input$ar_coeff
    }
    
    ### --- LAGGED PATHS ---
    pos_lag <- which(matrix(1:(v*v), v) > 0, arr.ind = TRUE)
    
    if (input$group_model_lag == "confirm") {
      group_lag_coords <- parse_indices(input$group_indices_lag)
      n_group_lag <- nrow(group_lag_coords)
      for (i in seq_len(n_group_lag)) {
        Phi[group_lag_coords[i, 1], group_lag_coords[i, 2]] <- rnorm(1, mean = input$group_lag_beta, sd = input$group_lag_sd)
      }
    } else {
      n_group_lag <- round(0.5 * n_lagged_active)
      selected_lag <- pos_lag[sample(nrow(pos_lag), n_group_lag), , drop = FALSE]
      for (i in seq_len(nrow(selected_lag))) {
        Phi[selected_lag[i, 1], selected_lag[i, 2]] <- rnorm(1, mean = input$group_lag_beta, sd = input$group_lag_sd)
      }
    }
    
    remaining_lag <- setdiff(1:nrow(pos_lag), selected_lag)
    n_indiv_lag <- n_lagged_active - n_group_lag
    selected_indiv_lag <- pos_lag[sample(remaining_lag, n_indiv_lag), , drop = FALSE]
    for (i in seq_len(nrow(selected_indiv_lag))) {
      Phi[selected_indiv_lag[i, 1], selected_indiv_lag[i, 2]] <- rnorm(1, mean = input$indiv_lag_beta, sd = input$indiv_lag_sd)
    }
    
    ### --- CONTEMPORANEOUS PATHS ---
    pos_con <- which(row(A) != col(A), arr.ind = TRUE)
    
    if (input$group_model_con == "confirm") {
      group_con_coords <- parse_indices(input$group_indices_con)
      n_group_con <- nrow(group_con_coords)
      for (i in seq_len(n_group_con)) {
        A[group_con_coords[i, 1], group_con_coords[i, 2]] <- rnorm(1, mean = input$group_con_beta, sd = input$group_con_sd)
      }
    } else {
      n_group_con <- round(0.5 * n_contemp_active)
      selected_con <- pos_con[sample(nrow(pos_con), n_group_con), , drop = FALSE]
      for (i in seq_len(nrow(selected_con))) {
        A[selected_con[i, 1], selected_con[i, 2]] <- rnorm(1, mean = input$group_con_beta, sd = input$group_con_sd)
      }
    }
    
    remaining_con <- setdiff(1:nrow(pos_con), selected_con)
    n_indiv_con <- n_contemp_active - n_group_con
    selected_indiv_con <- pos_con[sample(remaining_con, n_indiv_con), , drop = FALSE]
    for (i in seq_len(nrow(selected_indiv_con))) {
      A[selected_indiv_con[i, 1], selected_indiv_con[i, 2]] <- rnorm(1, mean = input$indiv_con_beta, sd = input$indiv_con_sd)
    }
    
    list(A = A, Phi = Phi)
  })
  
  simulate_timeseries <- eventReactive(input$run_sim, {
    mat <- generate_matrices()
    A <- mat$A
    Phi <- mat$Phi
    v <- nrow(A)
    t <- input$time_steps
    st <- t + 50
    
    noise <- matrix(rnorm(v * st, 0, input$overall_noise), v)
    I <- diag(v)
    time <- matrix(0, nrow = v, ncol = (st + 1))
    time1 <- matrix(0, nrow = v, ncol = st)
    
    for (i in 1:st) {
      time1[, i] <- solve(I - A) %*% (Phi %*% time[, i] + noise[, i])
      time[, i + 1] <- time1[, i]
    }
    
    time1 <- time1[, 51:(50 + t)]
    series <- t(time1)
    series <- round(series,digits=5)
    
    series <- droprows(series, input$drop_pct_row, input$target_median_block,input$drop_pct_random)
    
    colnames(series) <- paste0("V", 1:v)
    series_df <- as.data.frame(series)
    series_df$Time <- 1:t
    series_df
  })
  
  draw_matrix_plot <- function(mat, title, lag_labels = FALSE) {
    n <- nrow(mat)
    rownames(mat) <- paste0("V", 1:n)
    colnames(mat) <- if (lag_labels) paste0("V", 1:n, "_lag") else paste0("V", 1:n)
    melted <- melt(mat)
    melted$label <- ifelse(melted$value != 0, round(melted$value, 2), "")
    melted$Var1 <- factor(melted$Var1, levels = rev(levels(melted$Var1)))
    
    ggplot(melted, aes(Var2, Var1, fill = value)) +
      geom_tile(color = "black") +
      geom_text(aes(label = label), size = 4) +
      scale_fill_gradient(low = "white", high = if (lag_labels) "darkred" else "steelblue") +
      labs(title = title) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(hjust = 0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x.top = element_text(),
        axis.ticks.length = unit(0, "pt")
      ) +
      scale_x_discrete(position = "top")
  }
  
  output$phi_plot <- renderPlot({
    req(generate_matrices())
    draw_matrix_plot(generate_matrices()$Phi, "Lagged Matrix", lag_labels = TRUE)
  })
  
  output$a_plot <- renderPlot({
    req(generate_matrices())
    draw_matrix_plot(generate_matrices()$A, "Contemporaneous Matrix")
  })
  
  output$ts_plot <- renderPlot({
    req(simulate_timeseries())
    df <- simulate_timeseries()
    df[is.na(df)] <- 0
    df_melt <- melt(df, id.vars = "Time")
    ggplot(df_melt, aes(x = Time, y = value, color = variable)) +
      geom_line() +
      labs(title = "Simulated Time Series", y = "Value", x = "Time") +
      theme_minimal()
  })
}

# Run the app
shinyApp(ui = ui, server = server)
