library(shiny)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Define UI
ui <- fluidPage(
  titlePanel("SCR-B-analysis"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Configuration"),
      textInput("work_dir", "Working Directory (containing CH folders)", value = getwd()),
      actionButton("load_btn", "Load / Refresh Data", class = "btn-primary"),
      hr(),
      h4("Parameter Selection"),
      selectInput("param_select", "Choose Parameter:", choices = NULL),
      hr(),
      helpText("Ensure your working directory contains subfolders like '...CH1', '...CH2' with 'MLD' csv files.")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Trace Plots", 
                 plotOutput("trace_plot", height = "500px"),
                 h4("Detailed Chain Inspection"),
                 plotOutput("density_plot", height = "300px")
        ),
        tabPanel("Estimates & Stats",
                 h3("Parameter Estimates"),
                 tableOutput("est_table"),
                 hr(),
                 h3("Convergence Diagnostics (Gelman-Rubin)"),
                 verbatimTextOutput("gelman_print"),
                 plotOutput("gelman_plot", height = "300px")
        )
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Reactive values to store data
  values <- reactiveValues(
    chains = NULL,
    mcmc_list = NULL,
    params = NULL,
    summary_df = NULL
  )
  
  # Load Data Function
  observeEvent(input$load_btn, {
    req(input$work_dir)
    
    withProgress(message = 'Loading Chains...', value = 0, {
      
      main_path <- input$work_dir
      if(!dir.exists(main_path)) {
        showNotification("Directory does not exist!", type = "error")
        return()
      }
      
      # Robust directory detection
      chain_dirs <- grep("CH[0-9]+", list.files(main_path), value = TRUE)
      
      if(length(chain_dirs) == 0) {
        showNotification("No Chain directories (CH*) found in the specified path.", type = "error")
        return()
      }
      
      raw_chains <- list()
      
      incProgress(0.2, detail = "Reading chain files...")
      
      for(d in chain_dirs) {
        full_d <- file.path(main_path, d)
        mld_file <- list.files(full_d, pattern = "MLD", full.names = TRUE)
        
        if(length(mld_file) > 0) {
          # Read CSV
          dat <- read.csv(mld_file[1])
          raw_chains[[d]] <- dat
        }
      }
      
      if(length(raw_chains) == 0) {
        showNotification("No MLD files found in chain directories.", type = "error")
        return()
      }
      
      # Process Chains (Robust Windowing)
      nsamples <- sapply(raw_chains, nrow)
      min_n <- min(nsamples)
      
      # Assuming first few cols might be junk, matching Analysis script logic
      # Remove col 1 (beta.behave / index) and col 7 (Density) if present
      # We'll try to be smart: if names match "beta" or "Density", drop them.
      
      clean_chains <- lapply(raw_chains, function(x) {
        # Basic cleanup indices based on original script: c(-1, -7)
        # We will check if ncol is sufficient
        if(ncol(x) >= 7) {
          x <- x[, c(-1, -7)]
        } else {
          x <- x[, -1] # At least drop index
        }
        # Window to min length
        x <- x[1:min_n, ]
        return(as.mcmc(x))
      })
      
      values$mcmc_list <- as.mcmc.list(clean_chains)
      
      # Extract parameter names
      values$params <- colnames(clean_chains[[1]])
      
      # Update Dropdown
      updateSelectInput(session, "param_select", choices = values$params)
      
      showNotification(paste("Loaded", length(clean_chains), "chains with", min_n, "iterations."), type = "message")
      
      incProgress(1, detail = "Done")
    })
  })
  
  # -- Trace Plot --
  output$trace_plot <- renderPlot({
    req(values$mcmc_list, input$param_select)
    
    # Convert to DF for ggplot
    df_list <- lapply(seq_along(values$mcmc_list), function(i) {
      d <- as.data.frame(values$mcmc_list[[i]])
      d$Iteration <- 1:nrow(d)
      d$Chain <- as.factor(paste("Chain", i))
      d
    })
    
    big_df <- do.call(rbind, df_list)
    
    # Filter for selected param
    p_name <- input$param_select
    if(!p_name %in% names(big_df)) return(NULL)
    
    ggplot(big_df, aes(x = Iteration, y = .data[[p_name]], color = Chain)) +
      geom_line(alpha = 0.8) +
      theme_minimal() +
      labs(title = paste("Trace Plot:", p_name), y = "Value") +
      theme(legend.position = "bottom")
  })
  
  # -- Density Plot --
  output$density_plot <- renderPlot({
    req(values$mcmc_list, input$param_select)
    
    df_list <- lapply(seq_along(values$mcmc_list), function(i) {
      d <- as.data.frame(values$mcmc_list[[i]])
      d$Chain <- as.factor(paste("Chain", i))
      d
    })
    big_df <- do.call(rbind, df_list)
    p_name <- input$param_select
    
    ggplot(big_df, aes(x = .data[[p_name]], fill = Chain)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(title = paste("Posterior Density:", p_name), x = "Value")
  })
  
  # -- Estimates Table --
  output$est_table <- renderTable({
    req(values$mcmc_list)
    
    # Calculate summary for all params
    # Combine all chains
    combined <- do.call(rbind, values$mcmc_list)
    
    means <- apply(combined, 2, mean)
    sds <- apply(combined, 2, sd)
    q2.5 <- apply(combined, 2, quantile, probs=0.025)
    q50 <- apply(combined, 2, quantile, probs=0.5)
    q97.5 <- apply(combined, 2, quantile, probs=0.975)
    
    res <- data.frame(
      Parameter = names(means),
      Mean = means,
      SD = sds,
      `2.5%` = q2.5,
      Median = q50,
      `97.5%` = q97.5
    )
    
    # Filter if user wants ONLY the selected one, or show all? 
    # Usually better to highlight selected, but showing all is good for "Estimates" tab.
    # Let's show selected first, then others.
    
    if(!is.null(input$param_select)) {
      idx <- which(res$Parameter == input$param_select)
      if(length(idx)>0) {
        # Move selected to top
        res <- rbind(res[idx,], res[-idx,])
      }
    }
    res
  })
  
  # -- Gelman Diagnostics --
  output$gelman_print <- renderPrint({
    req(values$mcmc_list)
    if(length(values$mcmc_list) < 2) {
      cat("Need at least 2 chains for Gelman-Rubin diagnostics.")
    } else {
      gelman.diag(values$mcmc_list, multivariate = FALSE)
    }
  })
  
  output$gelman_plot <- renderPlot({
    req(values$mcmc_list)
    if(length(values$mcmc_list) >= 2) {
      gelman.plot(values$mcmc_list)
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
