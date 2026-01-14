source("HelperFunctions.R")
library(shiny)
library(ggplot2)
library(reshape2)
library(shinyBS)

# Parse user input index strings like "1,2; 2,3"
parse_indices <- function(txt) {
  if (nzchar(txt)) {
    pairs <- unlist(strsplit(txt, ";"))
    do.call(rbind, lapply(pairs, function(pair) as.integer(unlist(strsplit(trimws(pair), ",")))))
  } else {
    matrix(nrow = 0, ncol = 2)
  }
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Time Series Matrix Simulator"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n_people", "Number of People in Study:", value = 10, min = 1),
      numericInput("time_steps", "Number of Time Steps:", value = 100, min = 1),
      numericInput("n_vars", "Number of Variables (N):", value = 6, min = 2),
      bsTooltip("n_vars", "Contemporaneous total = N*(N-1), Lagged total = N*N.", "right", options = list(container = "body")),
      
      numericInput("network_density", "Total Network Density (Decimal):", value = 0.2, min = 0, max = 1, step = 0.01),
      bsTooltip("network_density", "What percentage of total possible paths do you want to be non-zero?", "right", options = list(container = "body")),
      
      numericInput("p_group", "Proportion of Group-Level Paths (Decimal):", value = 0.4, min = 0, max = 1, step = 0.01),
      bsTooltip("p_group", "What percentage of active paths do you want at the group level?", "right", options = list(container = "body")),
      
      hr(),
      h3("Lagged Matrix (Φ)"),
      textAreaInput("ar_paths", "User-defined AR Paths (row,col pairs)", placeholder = "1,1; 2,2"),
      numericInput("ar_coeff", "AR Coefficient:", value = 0.5),
      numericInput("ar_sd", "AR Coefficient SD:", value = 0.01),
      bsTooltip("ar_sd", "Gaussian noise is added to beta coefficients with this standard deviation", "right", options = list(container = "body")),
      
      radioButtons("group_model_lag", "Group-level Path Model (Lagged):",
                   choices = c("Random" = "random", "User-specified" = "confirm")),
      
      conditionalPanel(
        condition = "input.group_model_lag == 'confirm'",
        textAreaInput("group_indices_lag", "Group-level Lagged Paths (row,col pairs)", placeholder = "1,2; 3,4")
      ),
      
      numericInput("group_lag_beta", "Group-level Lagged Beta:", value = 0.3),
      numericInput("group_lag_sd", "Group-level Lagged Beta SD:", value = 0.01),
      bsTooltip("group_lag_sd", "Gaussian noise is added to beta coefficients with this standard deviation", "right", options = list(container = "body")),
      
      numericInput("indiv_lag_beta", "Individual-level Lagged Beta:", value = 0.3),
      numericInput("indiv_lag_sd", "Individual-level Lagged Beta SD:", value = 0.01),
      bsTooltip("indiv_lag_sd", "Gaussian noise is added to beta coefficients with this standard deviation", "right", options = list(container = "body")),
      
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
      bsTooltip("group_con_sd", "Gaussian noise is added to beta coefficients with this standard deviation", "right", options = list(container = "body")),
      
      numericInput("indiv_con_beta", "Individual-level Contemporaneous Beta:", value = 0.3),
      numericInput("indiv_con_sd", "Individual-level Contemporaneous Beta SD:", value = 0.01),
      bsTooltip("indiv_con_sd", "Gaussian noise is added to beta coefficients with this standard deviation", "right", options = list(container = "body")),
      
      hr(),
      h3("Timeseries Noise Settings"),
      numericInput("overall_noise", "Overall Timeseries Noise SD:", value = 0.1, min = 0),
      
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

# ---- SERVER ----
server <- function(input, output) {
  
  mat_obj <- eventReactive(input$run_sim, {
    mat.generate.asw(
      nvar = input$n_vars,
      ar.value = input$ar_coeff,
      dens = input$network_density,
      p.group = input$p_group,
      lag.b = input$group_lag_beta,
      con.b = input$group_con_beta,
      ar.sd = input$ar_sd,
      lag.sd = input$group_lag_sd,
      con.sd = input$group_con_sd,
      user_ar_indices = parse_indices(input$ar_paths),
      laggroup.random = (input$group_model_lag == "random"),
      cntgroup.random = (input$group_model_con == "random"),
      user_lag_indices = parse_indices(input$group_indices_lag),
      user_con_indices = parse_indices(input$group_indices_con)
    )
  })
  
  ts_obj <- eventReactive(input$run_sim, {
    ts.generate.asw(
      mat = mat_obj()$sub1,
      lvl = mat_obj()$lvl1,
      t = input$time_steps
    )
  })
  
  draw_matrix_plot <- function(mat, title, lag_labels = FALSE) {
    n <- nrow(mat)
    rownames(mat) <- paste0("V", 1:n)
    colnames(mat) <- if (lag_labels) paste0("V", 1:n, "_lag") else paste0("V", 1:n)
    melted <- melt(mat)
    melted$label <- ifelse(melted$value != 0, round(melted$value, 2), "")
    melted$Var1 <- factor(melted$Var1, levels = rev(levels(factor(melted$Var1))))
    
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
    req(ts_obj())
    draw_matrix_plot(ts_obj()$paths[, 1:input$n_vars], "Lagged Matrix", lag_labels = TRUE)
  })
  
  output$a_plot <- renderPlot({
    req(ts_obj())
    draw_matrix_plot(ts_obj()$paths[, (input$n_vars + 1):(2 * input$n_vars)], "Contemporaneous Matrix")
  })
  
  output$ts_plot <- renderPlot({
    req(ts_obj())
    df <- ts_obj()$series
    df <- as.data.frame(df)
    df$Time <- 1:nrow(df)
    df_melt <- melt(df, id.vars = "Time")
    ggplot(df_melt, aes(x = Time, y = value, color = variable)) +
      geom_line() +
      labs(title = "Simulated Time Series", y = "Value", x = "Time") +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)
