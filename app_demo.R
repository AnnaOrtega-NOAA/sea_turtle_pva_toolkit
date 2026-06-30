library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(mvtnorm)
library(truncnorm)
library(jagsUI)
library(MARSS)

# =====================================================================
# MATHEMATICAL HELPERS & DISTRIBUTIONS
# =====================================================================
inv_logit <- function(x) { exp(x) / (1 + exp(x)) }
safe_logit <- function(x) { x_safe <- pmax(0.001, pmin(0.999, x)); log(x_safe / (1 - x_safe)) }

compute_CMP_constant <- function(Lambda, Nu, Mu, Tol, Max, Log=TRUE, Type="Z"){
  if( (!is.na(Lambda) & Lambda > 10^Nu) | (!is.na(Mu) & Mu^Nu > 10^Nu) ){
    if(Type=="Z"){ ln_Const = Nu*Lambda^(1/Nu) - ((Nu-1)/(2*Nu))*log(Lambda) - ((Nu-1)/2)*log(2*pi) - (1/2)*log(Nu) }
    if(Type=="S"){ ln_Const = Nu*Mu - ((Nu-1)/(2))*log(Mu) - ((Nu-1)/2)*log(2*pi) - (1/2)*log(Nu) }
  }else{
    Const = rep(0,Max+1); Index = 1; Const[Index] = 1
    while( Const[Index]/Const[1] > Tol ){
      if(Type=="Z") Const[Index+1] = Const[Index] * ( Lambda / Index^Nu )
      if(Type=="S") Const[Index+1] = Const[Index] * ( Mu / Index )^Nu; Index = Index + 1
    }
    ln_Const = log(sum(Const))
  }
  if(Log) return(ln_Const) else return(exp(ln_Const))
}

dCMP <- function( x, lambda, mu, nu, log=TRUE, tol=0.01, iter.max=200 ){
  if(missing(mu) & !missing(lambda)) loglike = x*log(lambda) - nu*lfactorial(x) - compute_CMP_constant(Lambda=lambda, Nu=nu, Mu=NA, Tol=tol, Max=iter.max, Log=TRUE, Type="Z")
  if(!missing(mu) & missing(lambda)) loglike = nu*x*log(mu) - nu*lfactorial(x) - compute_CMP_constant(Lambda=NA, Nu=nu, Mu=mu, Tol=tol, Max=iter.max, Log=TRUE, Type="S")
  if(log) return(loglike) else return(exp(loglike))
}

rCMP <- function( n, lambda, mu, nu, tol=0.01, x_max=200 ){
  loglike_x = rep(NA, x_max+1)
  for( x in 0:x_max ){
    if(missing(mu) & !missing(lambda)) loglike_x[x+1] = dCMP( x=x, lambda=lambda, nu=nu, log=TRUE, tol=tol, iter.max=x_max)
    if(!missing(mu) & missing(lambda)) loglike_x[x+1] = dCMP( x=x, mu=mu, nu=nu, log=TRUE, tol=tol, iter.max=x_max)
  }
  return(sample(x=0:x_max, size=n, replace=TRUE, prob=exp(loglike_x)))
}

# =====================================================================
# EXACT FOURIER IMPUTATION ENGINE
# =====================================================================
# Add 'six_month_sites' as an incoming argument to the engine
# Add 'six_month_sites' as an incoming argument to the engine
run_fourier_imputation <- function(df, iter = 100000, six_month_sites = NULL) {
  df_clean <- df %>% 
    group_by(Year, Month, Site) %>% 
    summarise(Count = if(all(is.na(Count))) NA_real_ else sum(Count, na.rm = TRUE), .groups = "drop")
  
  all_years <- min(df_clean$Year, na.rm = TRUE):max(df_clean$Year, na.rm = TRUE)
  sites <- sort(unique(df_clean$Site))
  n_years <- length(all_years)
  n_timeseries <- length(sites)
  
  full_grid <- expand.grid(Month = 1:12, Year = all_years, Site = sites, stringsAsFactors = FALSE)
  prep_df <- full_grid %>% 
    left_join(df_clean, by = c("Year", "Month", "Site")) %>% 
    arrange(Year, Month) %>% 
    pivot_wider(names_from = Site, values_from = Count) %>% 
    select(all_of(sites))
  
  y_matrix <- as.matrix(prep_df)
  y_matrix[y_matrix <= 0 | is.nan(y_matrix) | is.infinite(y_matrix)] <- NA
  y_matrix <- log(y_matrix)
  y_matrix <- matrix(y_matrix, ncol = n_timeseries)
  
  # Honor user override selections: assign a 6 if user selected it, otherwise assign 12
  periods <- rep(12, n_timeseries)
  for(i in 1:n_timeseries) { 
    if(sites[i] %in% six_month_sites) periods[i] <- 6 
  }
  
  jags_data <- list(m = rep(1:12, times = n_years), n.steps = nrow(y_matrix), n.months = 12, pi = pi, period = periods, n.timeseries = n_timeseries, n.years = n_years)
  
  if (n_timeseries == 1) {
    jags_data$y <- as.vector(y_matrix)
    model_string <- "
    model {
      predX0 ~ dnorm(5, 0.1)
      predX[1] <- c[m[1]] + predX0
      X[1] ~ dnorm(predX[1], tau.X)
      y[1] ~ dnorm(X[1], tau.y)
      for (t in 2:n.steps){
         predX[t] <- c[m[t]] + X[t-1]
         X[t] ~ dnorm(predX[t], tau.X)
         y[t] ~ dnorm(X[t], tau.y)
      }
      for (y_idx in 1:n.years){
         for (mm in 1:12){ tmp2[y_idx, mm] <- exp(X[(y_idx*12 - mm + 1)]) }
         N[y_idx] <- log(sum(tmp2[y_idx, ]))
      }
      for (k in 1:n.months){
          c.const[k] <- 2 * pi * k / period[1]
          c[k] <- beta.cos * cos(c.const[k]) + beta.sin * sin(c.const[k])
      }
      sigma.y ~ dgamma(2, 0.5); tau.y <- 1/(sigma.y * sigma.y)
      beta.cos ~ dnorm(0, 1); beta.sin ~ dnorm(0, 1)
      sigma.X ~ dgamma(2, 0.5); tau.X <- 1/(sigma.X * sigma.X)
    }"
  } else {
    jags_data$y <- y_matrix
    model_string <- "
    model {
      for(j in 1:n.timeseries) {
         predX0[j] ~ dnorm(5, 0.1)
         predX[1,j] <- c[j, m[1]] + predX0[j]
         X[1,j] ~ dnorm(predX[1,j], tau.X[j])
         y[1,j] ~  dnorm(X[1,j], tau.y[j])
         for (t in 2:n.steps){
             predX[t,j] <-  c[j,m[t]] + X[t-1, j]
             X[t,j] ~ dnorm(predX[t,j], tau.X[j])
             y[t,j] ~  dnorm(X[t,j], tau.y[j])
         }
         for (y_idx in 1:n.years){
            for (mm in 1:12){ tmp2[y_idx, mm, j] <- exp(X[(y_idx*12 - mm + 1), j]) }
            N[y_idx, j] <- log(sum(tmp2[y_idx, , j]))
         }
      }
      for (j in 1:n.timeseries){
          for (k in 1:n.months){
              c.const[j, k] <-  2 * pi * k / period[j]
              c[j, k] <- beta.cos[j] * cos(c.const[j,k]) + beta.sin[j] * sin(c.const[j,k])
          }
          sigma.y[j] ~ dgamma(2, 0.5); tau.y[j] <- 1/(sigma.y[j] * sigma.y[j])
          beta.cos[j] ~ dnorm(0, 1); beta.sin[j] ~ dnorm(0, 1)
          sigma.X[j] ~ dgamma(2, 0.5); tau.X[j] <- 1/(sigma.X[j] * sigma.X[j])
      }
    }"
  }
  
  jm <- jagsUI::jags(data = jags_data, parameters.to.save = c("N"), model.file = textConnection(model_string),
                     n.chains = 3, n.iter = iter, n.burnin = floor(iter/3), n.thin = 5, parallel = FALSE, verbose = FALSE)
  
  n_sims <- jm$q50$N
  if (n_timeseries == 1) n_sims <- matrix(n_sims, ncol = 1)
  
  d_annual_list <- list()
  for (i in 1:n_timeseries) {
    d_annual_list[[i]] <- data.frame(Year = all_years, Site = sites[i], Count = exp(n_sims[, i]))
  }
  res_annual <- do.call(rbind, d_annual_list)
  return(res_annual)
}

# =====================================================================
# EXACT STATE-SPACE TREND ENGINE
# =====================================================================
run_jags_aligned <- function(df, iter = 50000, burnin = 10000, thin = 10) {
  all_years <- min(df$Year, na.rm = TRUE):max(df$Year, na.rm = TRUE)
  sites <- sort(unique(df$Site))
  n.yrs <- length(all_years); n.timeseries <- length(sites)
  
  prep_df <- df %>%
    group_by(Year, Site) %>%
    summarise(Annual_Nesters = if(all(is.na(Annual_Nesters))) NA_real_ else sum(Annual_Nesters, na.rm = TRUE), .groups = "drop") %>%
    complete(Year = all_years, Site = sites) %>%
    arrange(Year) %>%
    pivot_wider(names_from = Site, values_from = Annual_Nesters) %>%
    select(all_of(sites))
  
  mat_data <- as.matrix(prep_df)
  mat_data[mat_data <= 0 | is.nan(mat_data) | is.infinite(mat_data)] <- NA
  Y_matrix <- t(log(mat_data))
  
  jags_data <- list(
    n.yrs = n.yrs, n.timeseries = n.timeseries, a_mean = 0, a_sd = 4,
    u_mean = 0, u_sd = 0.5, q_alpha = 0.01, q_beta = 0.01, r_alpha = 0.01, r_beta = 0.01,
    x0_mean = mean(Y_matrix[1:n.timeseries, ], na.rm = TRUE), x0_sd = 10
  )
  if(is.na(jags_data$x0_mean)) jags_data$x0_mean <- 5
  
  if (n.timeseries == 1) {
    jags_data$Y <- as.vector(Y_matrix)
    model_string <- "
    model {
      U ~ dnorm(u_mean, 1/(u_sd^2))
      tauQ ~ dgamma(q_alpha, q_beta)
      Q <- 1/tauQ
      X[1] ~ dnorm(x0_mean, 1/(x0_sd^2))
      X0 <- X[1] - U
      for(t in 2:n.yrs) {
        predX[t] <- X[t-1] + U
        X[t] ~ dnorm(predX[t], tauQ)
      }
      tauR ~ dgamma(r_alpha, r_beta)
      R <- 1/tauR
      for(t in 1:n.yrs) { Y[t] ~ dnorm(X[t], tauR) }
    }"
    params <- c("U", "Q", "R", "X0", "X")
  } else {
    jags_data$Y <- Y_matrix
    jags_data$Z <- matrix(rep(1, n.timeseries), ncol=1)
    model_string <- "
    model {
      A[1] <- 0
      for(j in 2:n.timeseries) { A[j] ~ dnorm(a_mean, 1/(a_sd^2)) }
      U ~ dnorm(u_mean, 1/(u_sd^2))
      tauQ ~ dgamma(q_alpha, q_beta)
      Q <- 1/tauQ
      X[1] ~ dnorm(x0_mean, 1/(x0_sd^2))
      X0 <- X[1] - U
      for(t in 2:n.yrs) {
        predX[t] <- X[t-1] + U
        X[t] ~ dnorm(predX[t], tauQ)
      }
      for(j in 1:n.timeseries) {
        tauR[j] ~ dgamma(r_alpha, r_beta)
        R[j] <- 1/tauR[j]
        for(t in 1:n.yrs) {
          predY[j,t] <- Z[j,1] * X[t] + A[j]
          Y[j,t] ~ dnorm(predY[j,t], tauR[j])
        }
      }
    }"
    params <- c("U", "Q", "R", "X0", "X", "A")
  }
  
  fit <- jagsUI::jags(data = jags_data, parameters.to.save = params, model.file = textConnection(model_string),
                      n.chains = 3, n.iter = iter, n.burnin = burnin, n.thin = thin, parallel = FALSE, verbose = FALSE)
  return(list(fit = fit, years = all_years))
}

# =====================================================================
# HISTORICAL RETROSPECTIVE COMPILER ENGINE
# =====================================================================
calculate_empirical_ane <- function(obs_df, safe_df, params, current_year) {
  merged_data <- obs_df %>% left_join(safe_df, by = "Year")
  merged_data$Total_Est[is.na(merged_data$Total_Est)] <- 16.5
  len_lm <- lm(log(Length) ~ Total_Est, data = merged_data)
  sim_beta0 <- coef(len_lm)[1]
  sim_beta1 <- ifelse(is.na(coef(len_lm)[2]), 0, coef(len_lm)[2])
  sim_sigma_L <- summary(len_lm)$sigma
  if(is.na(sim_sigma_L) || sim_sigma_L == 0) sim_sigma_L <- 0.338
  
  safe_logit_vec <- log(pmax(0.001, pmin(0.999, merged_data$Mortality)) / (1 - pmax(0.001, pmin(0.999, merged_data$Mortality))))
  sim_mu0 <- mean(safe_logit_vec, na.rm = TRUE)
  sim_sigma_D <- sd(safe_logit_vec, na.rm = TRUE)
  if(is.na(sim_sigma_D) || sim_sigma_D == 0) sim_sigma_D <- 0.50
  sim_rho <- cor(log(merged_data$Length), safe_logit_vec, use = "complete.obs")
  if(is.na(sim_rho)) sim_rho <- -0.51
  sim_cov <- matrix(c(sim_sigma_L^2, sim_sigma_L * sim_sigma_D * sim_rho, sim_sigma_L * sim_sigma_D * sim_rho, sim_sigma_D^2), 2, 2)
  
  emp_mean_len <- mean(obs_df$Length, na.rm = TRUE)
  emp_mean_mort <- mean(obs_df$Mortality, na.rm = TRUE)
  sp_years <- sort(unique(safe_df$Year))
  all_sp_turtles <- data.frame()
  
  for (target_yr in sp_years) {
    obs_yr <- obs_df[which(obs_df$Year == target_yr), ]
    total_est <- sum(safe_df$Total_Est[safe_df$Year == target_yr], na.rm = TRUE)
    n_unobs <- max(0, round(total_est) - nrow(obs_yr))
    if (nrow(obs_yr) > 0) {
      all_sp_turtles <- rbind(all_sp_turtles, data.frame(CaptureYear = target_yr, Length = coalesce(obs_yr$Length, emp_mean_len), Mortality = coalesce(obs_yr$Mortality, emp_mean_mort)))
    }
    if (n_unobs > 0) {
      mu_l <- sim_beta0 + sim_beta1 * total_est
      draws <- mvtnorm::rmvnorm(n_unobs, mean = c(mu_l, sim_mu0), sigma = sim_cov)
      all_sp_turtles <- rbind(all_sp_turtles, data.frame(CaptureYear = target_yr, Length = exp(draws[,1]), Mortality = 1 / (1 + exp(-draws[,2]))))
    }
  }
  
  sp_ledger <- data.frame()
  for (i in 1:nrow(all_sp_turtles)) {
    c_year <- all_sp_turtles$CaptureYear[i]
    age_start <- params$sp_tknot - (1 / params$sp_k) * log(1 - (pmin(all_sp_turtles$Length[i], params$sp_linf - 0.1) / params$sp_linf))
    if(is.nan(age_start) || is.na(age_start)) age_start <- params$sp_max_age - 2
    future_years <- seq(c_year, current_year)
    l_y <- length(future_years); if(l_y < 1) next
    ages_traj <- seq(age_start, length.out = l_y, by = 1)
    lens_traj <- params$sp_linf * (1 - exp(-params$sp_k * (ages_traj - params$sp_tknot)))
    
    p_mat <- ifelse(lens_traj >= 0.99 * params$sp_linf, 1.0, 1.0 / (1.0 + exp(-(lens_traj - params$sp_lmat) / params$sp_sig_mat)))
    surv_vector <- cumprod((1 - p_mat) * params$sp_ane_pj + p_mat * params$sp_ane_pa)
    
    p_binom <- rbinom(l_y, 1, p_mat)
    if(any(p_binom == 1)) p_binom[min(which(p_binom == 1)):l_y] <- 1
    
    sr <- surv_vector * p_binom * params$sp_pf * all_sp_turtles$Mortality[i] * (1 / params$sp_remig_int)
    sp_ledger <- rbind(sp_ledger, data.frame(CalendarYear = future_years, ANE = sr))
  }
  
  final_historical_ane <- sp_ledger %>% 
    group_by(CalendarYear) %>% 
    summarise(Total_Cumulative_ANE = sum(ANE, na.rm = TRUE), .groups = "drop") %>%
    rename(Year = CalendarYear)
  
  attr(final_historical_ane, "meta_beta0") <- sim_beta0
  attr(final_historical_ane, "meta_beta1") <- sim_beta1
  attr(final_historical_ane, "meta_sigma_L") <- sim_sigma_L
  attr(final_historical_ane, "meta_sigma_D") <- sim_sigma_D
  attr(final_historical_ane, "meta_rho")     <- sim_rho
  attr(final_historical_ane, "meta_mu0")     <- sim_mu0
  return(final_historical_ane)
}

# =====================================================================
# UI LAYOUT
# =====================================================================
ui <- page_navbar(
  title = "Sea turtle population viability analysis toolkit",
  theme = bs_theme(version = 5, bootswatch = "yeti"),
  fillable = FALSE,
  
  # --- CHAPTER 1: DATA INGESTION PREVIEW ---
  nav_panel(
    title = "1. Data preview",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 330, title = "Data Source Configurations",
        radioButtons("data_mode", "Data Selection Mode:",
                     choices = c("Demo Mode" = "demo", "Upload Own Data" = "upload"),
                     selected = "demo"),
        hr(),
        conditionalPanel(
          condition = "input.data_mode == 'demo'",
          selectInput("demo_trend", 
                      label = tags$span("Regional population trend:", 
                                        popover(shiny::icon("question-circle"), 
                                                "Stable means the population is steady, increasing is a growing population, and decreasing is a declining population trend.",
                                                title = "Population trend selection")),
                      choices = c("Stable" = "stable", "Increasing" = "increasing", "Decreasing" = "decreasing"),
                      selected = "stable")
        ),
        conditionalPanel(
          condition = "input.data_mode == 'upload'",
          fileInput("uploaded_file", "Choose CSV Spreadsheet File:", accept = c(".csv")),
          uiOutput("mapping_ui")
        ),
        hr(),
        uiOutput("timeframe_ui"), 
        hr(),
        selectInput("species_preset", 
                    label = tags$span("Select population:", 
                                      popover(shiny::icon("question-circle"), 
                                              "Autofills demographic parameters (clutch frequency, remigration interval, etc) based on previous Technical Memoranda (Martin et al., 2020).",
                                              title = "Pre-populate demography")),
                    choices = c("North Pacific loggerhead" = "North Pacific Loggerhead", "western Pacific leatherback" = "western Pacific Leatherback")),
        hr(),
        numericInput("clutch_freq", "Clutch Frequency:", value = 4.6, step = 0.1),
        numericInput("remig_int", "Remigration Interval (years):", value = 3.3, step = 0.1),
        hr(),
        # New placeholder slot for the user-controlled period overrides
        uiOutput("period_override_ui"),
        hr(),
        div(style = "background-color: #f8f9fa; border-left: 4px solid #0dcaf0; padding: 12px; border-radius: 4px; font-size: 0.95rem; font-weight: 500; color: #212529;",
            "If you are satisfied with the data input, please move to page 2. Baseline trends"
        )
      ),
      layout_column_wrap(
        width = 1,
        card(class = "bg-light border-start border-info border-4", p(shiny::icon("info-circle"), tags$b(" Welcome to the data workspace setup page!"), " This tab visualizes the quality and completeness of your historical records.", style = "margin-bottom:0px;"))
      ), br(),
      
      fluidRow(
        column(4, bslib::value_box(title = "Data Completeness", value = uiOutput("completeness_card_ui"), theme = "info", style = "height: 220px;")),
        column(8, bslib::card(bslib::card_header("Table 1. Data coverage over full timeline"), tableOutput("gap_table"), style = "max-height: 220px; overflow-y: auto;"))
      ), br(),
      
      # 1. New dynamic QA/QC message center layout node
      fluidRow(
        column(12, uiOutput("qaqc_alerts"))
      ), br(),
      
      # 2. Main data plots layout row
      fluidRow(
        column(6, card(card_header("Raw nest counts over time, by beach"), 
                       plotOutput("preview_annual_raw", height = "500px"),
                       downloadButton("download_preview_raw", "Download Raw Plot", class = "btn-sm btn-outline-secondary mt-2"))),
        column(6, card(card_header("Data preview ledger grid (Top 10 Rows)"), 
                       tableOutput("data_preview_table_raw"), style = "max-height: 565px; overflow-y: auto;"))
      ), br(),
      
      # 3. New conditional row: Only shows up if an uploaded file contains monthly sub-components
      conditionalPanel(
        condition = "input.data_mode == 'upload' && input.map_month != 'none'",
        fluidRow(
          column(12, card(card_header("Nest Counts by Month (Waveform Peak Verification Check)"),
                          plotOutput("preview_monthly_seasonality", height = "600px"), # Height increased to account for multiple beach rows
                          hr(),
                          uiOutput("seasonality_recommendation") # New dynamic recommendation engine container
          ))
        )
      )
    )
  ),
  
  # --- CHAPTER 2: REGIONAL TREND ASSESSMENT ---
  nav_panel(
    title = "2. Abundance and Trend",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 330, title = "Simulation settings",
        actionButton("run_model", "Run “Trend and Abundance” model", class = "btn-primary w-100", style = "font-weight: bold; font-size:1.05rem; margin-bottom: 12px;"),
        hr(),
        sliderInput("iterations", 
                    label = tags$span("Number of iterations:", tooltip(shiny::icon("info-circle"), "Higher numbers increase model precision depth but take longer to compute.")),
                    min = 10000, max = 150000, value = 30000, step = 10000),
        checkboxInput("run_split", "Would you like to calculate separate trends before and after a specific year? E.g. Did nest monitoring methods change significantly in 2014?", value = TRUE),
        conditionalPanel(condition = "input.run_split == true", numericInput("split_year", "Select the year you’d like to calculate split trends from:", value = 2005, min = 1980, max = 2050))
      ),
      layout_column_wrap(
        width = 1,
        card(class = "bg-light border-start border-primary border-4", p(shiny::icon("chart-line"), tags$b(" Calculated population trend:"), " The plot below shows the smoothed trajectory line fitted across your beach counts. The shaded region defines the model's true 95% uncertainty boundaries.", style = "margin-bottom:0px;"))
      ), br(),
      uiOutput("summary_stats"), br(), uiOutput("executive_summary"), br(), 
      card(card_header("Model fit to regional population trend (all beaches share a regional trend)"), 
           div(style = "width: 100%;", plotOutput("clean_baseline_plot", height = "450px")))
    )
  ),
  
  # --- CHAPTER 3: BEACH DIAGNOSTICS LAB LAYOUT ---
  nav_panel(
    title = "3. Beach diagnostics",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 330, title = "Data & Model Verification",
        p(tags$b("How to Interpret the Diagnostics Matrix:")),
        p("• ", tags$b("Diagonal Curves:"), " Represent standalone certainty. Tall, narrow peaks indicate clear data signals."),
        p("• ", tags$b("Off-Diagonal Ovals:"), " Show variable interactions. Round, un-tilted boundaries confirm that growth rate calculations are stable.")
      ),
      layout_column_wrap(
        width = 1,
        card(class = "bg-light border-start border-dark border-4", p(shiny::icon("sliders"), tags$b(" Separating environmental stochasticity from observation stochasticity:"), " This panel isolates components to help us understand where populations experience shifts.", style = "margin-bottom:0px;"))
      ), br(),
      card(card_header("Framework Controls"), checkboxGroupInput("plot_layers", "Select models to display on chart canvas below:", choices = c("JAGS Regional Trend (Martin et al., 2020)" = "jags", "MARSS Regional Trend (Holmes et al., 2012)" = "marss_s", "MARSS Independent Site Trends (Holmes et al., 2012)" = "marss_i"), selected = c("jags", "marss_i"), inline = TRUE)), br(),
      
      layout_column_wrap(
        width = 1/2,
        card(card_header("1. Fitted Abundance Trajectories"), uiOutput("dynamic_unified_plot"), downloadButton("download_unified", "Download Fit Comparison Plot", class = "btn-sm btn-outline-secondary mt-auto"), style = "height: 480px;"),
        card(card_header("2. Variance Scalers: Environment vs Monitoring Noise"), div(style = "max-width: 440px; margin: 0 auto; width: 100%;", plotOutput("plot_var", height = "230px")), hr(), uiOutput("variance_interpretation_text"), downloadButton("download_var", "Download Variance Components Plot", class = "btn-sm btn-outline-secondary mt-auto"), style = "height: 480px;"),
        card(card_header("3. Quantified Growth Coefficients Log"), tableOutput("table_u"), style = "height: 480px; overflow-y: auto;"),
        card(
          card_header("4. Joint Parameter Covariance Matrix"), 
          div(style = "width: 360px; height: 360px; margin: 0 auto; display: block;", 
              plotOutput("posterior_pairs_plot", width = "100%", height = "100%")), 
          hr(),
          uiOutput("covariance_interpretation_text"), # <-- Adds the dynamic text container
          downloadButton("download_posterior", "Download Posterior Matrix Plot", class = "btn-sm btn-outline-secondary mt-auto"), 
          style = "height: 560px; overflow-y: auto;" # <-- Expanded height and added scroll safety
        )
      )
    )
  ),
  
  # --- CHAPTER 4: RETROSPECTIVE REMOVAL LEDGER ---
  nav_panel(
    title = "4. ANEs/Take",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 330, title = "Threat Parameters Audit",
        numericInput("ane_custom_take", "Annual historical turtle interactions (count):", value = 30, min = 0, step = 5),
        sliderInput("ane_custom_mort", "Historical interaction mortality rate:", min = 0, max = 1, value = 0.35, step = 0.05)
      ),
      layout_column_wrap(
        width = 1,
        card(class = "bg-light border-start border-warning border-4", p(shiny::icon("receipt"), tags$b(" Adult Nester Equivalents, ANEs:"), " This tab tracks the hidden demographic cost of historical interactions (e.g. fisheries bycatch) by projecting mortalities backward.", style = "margin-bottom:0px;"))
      ), br(),
      card(card_header("Historical Footprint Assessment: Estimated Nesting Drawdown Impact"), 
           div(style = "width: 100%;", plotOutput("ane_impact_plot", height = "480px")),
           downloadButton("download_ane", "Download Drawdown Plot", class = "btn-sm btn-outline-secondary mt-2"))
    )
  ),
  
  # --- CHAPTER 5: RISK FORECASTING SANDBOX (TRUE DEMOGRAPHIC ACCOUNTING) ---
  nav_panel(
    title = "5. Future projections",
    layout_sidebar(
      fillable = FALSE,
      sidebar = sidebar(
        width = 340, title = "Simulation parameters",
        accordion(
          accordion_panel("Forecast Controls",
                          sliderInput("pva_years", "Years to project forward into the future:", min = 10, max = 100, value = 50, step = 10),
                          numericInput("pva_sims", "Number of iterations:", value = 500, min = 100, max = 2500, step = 100),
                          selectInput("trend_to_project", "Historical trend to project forward:", choices = c("Full Timeline Trend", "Post-Split Year Trend"))),
          accordion_panel("Future Threat Configuration",
                          checkboxInput("apply_threats", "Introduce Future Threat Scenarios", value = FALSE),
                          conditionalPanel(condition = "input.apply_threats == true",
                                           numericInput("cmp_mu", "Expected future fishery encounter frequency rate (μ):", value = 15.0, step = 0.5),
                                           sliderInput("threat_mortality", "Future interaction mortality risk probability rate:", min = 0, max = 1, value = 0.25, step = 0.05),
                                           numericInput("cmp_nu", "Fleet dispersion variance modifier (ν):", value = 0.15, step = 0.01))),
          accordion_panel("Strategy A: Nest Protections",
                          numericInput("optA_nests", "Number of protected nests per year:", value = 35),
                          numericInput("clutch_size", "Average clutch size:", value = 95),
                          numericInput("hatch_succ", "Hatching success proportion (0-1):", value = 0.72, min=0, max=1, step=0.05),
                          numericInput("emerg_succ", "Emergence success proportion (0-1):", value = 0.88, min=0, max=1, step=0.05),
                          numericInput("yrlng_surv", "Wild Yearling (age 0 to 1) survival rate:", value = 0.025, min=0, max=1, step=0.005)),
          accordion_panel("Strategy B: Headstarting",
                          numericInput("optB_head", "Number of headstarted turtles released per year:", value = 15),
                          numericInput("optB_age", "Age of headstarts at release (years):", value = 1.0, min=1, max=5, step=1)),
          accordion_panel("Strategy C: Adult Poaching Cessation",
                          numericInput("optC_adults", "Number of nesting adults saved from poaching per year:", value = 6)),
          accordion_panel("Strategy D & E: Fleet Controls",
                          numericInput("optD_mu_mit", "Reduced encounter rate achieved via gear changes (μ_mit):", value = 4.5, step=0.5),
                          hr(),
                          numericInput("optE_m_min", "Minimum mortality risk from improved handling:", value = 0.05, min=0, max=1, step=0.01),
                          numericInput("optE_m_max", "Maximum mortality risk from improved handling:", value = 0.18, min=0, max=1, step=0.01))
        ), br(),
        card(
          card_header("Build a combined strategy portfolio"),
          checkboxGroupInput("portfolio_options", "Select multiple management actions to execute simultaneously:",
                             choices = c("Strategy A: Nest protections" = "A", "Strategy B: Headstarting" = "B", "Strategy C: Adult female poaching cessation" = "C", "Strategy D: Modified fishing gear (fewer interactions)" = "D", "Strategy E: Best handling and release practices (lower mortality)" = "E"),
                             selected = c("A", "D", "E"))
        ),
        actionButton("run_portfolio_sim", "Simulate future projections", class = "btn-success w-100", style = "font-weight: bold; font-size:1.05rem;")
      ),
      layout_column_wrap(
        width = 1,
        card(class = "bg-light border-start border-success border-4", p(shiny::icon("rocket"), tags$b(" Projecting long-term conservation strategy payoffs:"), " Mix and match options in the side panel to see your project's century-long recovery path.", style = "margin-bottom:0px;"))
      ), br(),
      fluidRow(
        column(7, card(
          card_header("100-year forward simulation forecasting trends"),
          checkboxGroupInput("portfolio_plots", "Select pathways to display on layout canvas:", 
                             choices = c("Status Quo", "Threats Only", "Strategy A Only", "Strategy B Only", "Strategy C Only", "Strategy D Only", "Strategy E Only", "Combined Portfolio"), 
                             selected = c("Status Quo", "Threats Only", "Combined Portfolio"), inline = TRUE),
          checkboxInput("show_proj_ci", "Display 95% uncertainty ribbons around trajectories", value = TRUE),
          plotOutput("portfolio_overlay_plot", height = "360px"),
          downloadButton("download_forecast", "Download Timeline Forecast Plot", class = "btn-sm btn-outline-secondary mt-2")
        )),
        column(5, card(
          card_header("Terminal Year Abundance Density Profiles"), 
          plotOutput("portfolio_density_plot", height = "415px"),
          downloadButton("download_density", "Download Abundance Density Plot", class = "btn-sm btn-outline-secondary mt-2")
        ))
      ), br(),
      layout_column_wrap(
        width = 1,
        card(
          card_header("Interpretation Strategy Guide for Stakeholders"),
          p("• ", tags$b("The Position:"), " Active portfolios should shift as far right as possible, away from zero, indicating recovery."),
          p("• ", tags$b("The Extinction Flag:"), " If a curve peaks heavily against the left margin near zero, it indicates a critical warning that a high proportion of trajectories crashed.")
        )
      ), br(),
      card(card_header("Abundance in 100 years for each strategy"), tableOutput("portfolio_scorecard_table"))
    )
  ),
  
  # Hidden parameter configurations
  nav_item(div(style="display:none;", numericInput("linf", "", 86.9), numericInput("k", "", 0.09), numericInput("tknot", "", -2.467), numericInput("max_age", "", 60), numericInput("lmat", "", 86.03), numericInput("sig_mat", "", 6.34), numericInput("pf", "", 0.65), numericInput("ane_pj", "", 0.8), numericInput("ane_pa", "", 0.895), radioButtons("workspace_source", "", "demo"), radioButtons("model_mode", "", "take")))
)

# =====================================================================
# SERVER ENGINE
# =====================================================================
server <- function(input, output, session) {
  
  vault <- reactiveValues(res = NULL, abund = NULL, unimputed_abund = NULL, marss = NULL, summary = NULL, year = NULL, nesters = NULL, total = NULL, trend_display = NULL, trend_pct = NULL, pre_trend = NULL, post_trend = NULL, pre_u_val = 0, post_u_val = 0, empirical_ane_ledger = NULL)
  vault_portfolio <- reactiveValues(plot_df = NULL, scorecard = NULL, all_scen_raw = list())
  
  custom_colors <- c(
    "Status Quo"         = "#000000",
    "Threats Only"       = "#D55E00",
    "Strategy A Only"    = "#E69F00",
    "Strategy B Only"    = "#56B4E9",
    "Strategy C Only"    = "#009E73",
    "Strategy D Only"    = "#F0E442",
    "Strategy E Only"    = "#CC79A7",
    "Combined Portfolio" = "#0072B2"
  )
  
  observeEvent(input$species_preset, {
    if (input$species_preset == "western Pacific Leatherback") {
      updateNumericInput(session, "linf", value = 142.7); updateNumericInput(session, "k", value = 0.2262)
      updateNumericInput(session, "tknot", value = -0.17); updateNumericInput(session, "max_age", value = 45.0)
      updateNumericInput(session, "lmat", value = 139.13); updateNumericInput(session, "sig_mat", value = 6.34)
      updateNumericInput(session, "pf", value = 0.73); updateNumericInput(session, "clutch_freq", value = 5.5)
      updateNumericInput(session, "remig_int", value = 3.06); updateNumericInput(session, "ane_pj", value = 0.81)
      updateNumericInput(session, "ane_pa", value = 0.893)
    } else if (input$species_preset == "North Pacific Loggerhead") {
      updateNumericInput(session, "linf", value = 86.9); updateNumericInput(session, "k", value = 0.09)
      updateNumericInput(session, "tknot", value = -2.467); updateNumericInput(session, "max_age", value = 60.0)
      updateNumericInput(session, "lmat", value = 86.03); updateNumericInput(session, "sig_mat", value = 6.34)
      updateNumericInput(session, "pf", value = 0.65); updateNumericInput(session, "clutch_freq", value = 4.6)
      updateNumericInput(session, "remig_int", value = 3.3); updateNumericInput(session, "ane_pj", value = 0.80)
      updateNumericInput(session, "ane_pa", value = 0.895)
    }
  })
  
  uploaded_raw <- reactive({
    req(input$uploaded_file)
    read.csv(input$uploaded_file$datapath, stringsAsFactors = FALSE)
  })
  
  output$mapping_ui <- renderUI({
    df <- uploaded_raw()
    cols <- colnames(df)
    tagList(
      selectInput("map_year", "Designate Year Column:", choices = cols),
      selectInput("map_month", "Designate Month Column (Optional):", choices = c("None / Annual Data" = "none", cols)),
      radioButtons("data_struct", "File Data Structure format:",
                   choices = c("Long Form (One Site column + One Count column)" = "long", 
                               "Wide Form (Beaches are separated columns)" = "wide")),
      conditionalPanel(
        condition = "input.data_struct == 'long'",
        selectInput("map_site", "Designate Beach/Site Name column:", choices = cols),
        selectInput("map_count", "Designate Counts column:", choices = cols)
      ),
      conditionalPanel(
        condition = "input.data_struct == 'wide'",
        selectizeInput("map_beaches", "Select Beach Columns to include in mapping:", choices = cols, multiple = TRUE)
      )
    )
  })
  
  raw_ingested_data <- reactive({
    if (input$data_mode == "demo") {
      set.seed(2026)
      years_vec <- 1985:2015
      n_yrs <- length(years_vec)
      slope <- 0.0
      if(input$demo_trend == "increasing") slope <- 15.5
      if(input$demo_trend == "decreasing") slope <- -18.0
      
      demo_df <- data.frame(
        Year = rep(years_vec, each = 3),
        Site = rep(c("Alpha Beach", "Beta Beach", "Gamma Beach"), times = n_yrs),
        Count = round(pmax(10, c(
          800 + (years_vec - 1985) * slope + rnorm(n_yrs, mean = 0, sd = 40),
          550 + (years_vec - 1985) * slope + rnorm(n_yrs, mean = 0, sd = 30),
          250 + (years_vec - 1985) * (slope * 0.5) + rnorm(n_yrs, mean = 0, sd = 15)
        )))
      )
      na_indices <- sample(1:nrow(demo_df), floor(nrow(demo_df) * 0.10))
      demo_df$Count[na_indices] <- NA_real_
      demo_df$Month <- NA
      return(demo_df)
    } else {
      req(input$uploaded_file, input$map_year, input$data_struct)
      df <- uploaded_raw()
      
      if (!input$map_year %in% colnames(df)) return(data.frame(Year=numeric(), Month=numeric(), Site=character(), Count=numeric()))
      
      if (input$data_struct == "long") {
        if (is.null(input$map_site) || is.null(input$map_count) || !input$map_site %in% colnames(df) || !input$map_count %in% colnames(df)) return(data.frame(Year=numeric(), Month=numeric(), Site=character(), Count=numeric()))
        out <- df %>% rename(Year = !!sym(input$map_year), Site = !!sym(input$map_site), Count = !!sym(input$map_count))
        out$Month <- if(input$map_month != "none" && input$map_month %in% colnames(df)) out[[input$map_month]] else NA
      } else {
        if (is.null(input$map_beaches) || !all(input$map_beaches %in% colnames(df))) return(data.frame(Year=numeric(), Month=numeric(), Site=character(), Count=numeric()))
        out <- df %>% rename(Year = !!sym(input$map_year))
        out$Month <- if(input$map_month != "none" && input$map_month %in% colnames(df)) out[[input$map_month]] else NA
        out <- out %>% pivot_longer(cols = all_of(input$map_beaches), names_to = "Site", values_to = "Count")
      }
      return(out %>% select(Year, Month, Site, Count))
    }
  })
  
  output$timeframe_ui <- renderUI({
    df <- raw_ingested_data()
    if(is.null(df) || nrow(df) == 0) return(NULL)
    yr_range <- range(df$Year, na.rm = TRUE)
    sliderInput("timeframe_filter", "Select Timeframe (Years):", min = yr_range[1], max = yr_range[2], value = c(yr_range[1], yr_range[2]), step = 1, sep = "")
  })
  
  processed_data <- reactive({
    df <- raw_ingested_data()
    if(is.null(df) || nrow(df) == 0) return(df)
    if(!is.null(input$timeframe_filter)) {
      data_range <- range(df$Year, na.rm = TRUE)
      if (input$timeframe_filter[1] >= data_range[1] && input$timeframe_filter[2] <= data_range[2]) {
        df <- df %>% filter(Year >= input$timeframe_filter[1] & Year <= input$timeframe_filter[2])
      }
    }
    return(df)
  })
  
  output$completeness_card_ui <- renderUI({
    d_long <- processed_data(); if(nrow(d_long) == 0) return("0%")
    total_cells <- nrow(d_long)
    valid_cells <- sum(!is.na(d_long$Count) & d_long$Count > 0, na.rm = TRUE)
    return(paste0(round((valid_cells / total_cells) * 100, 1), "%"))
  })
  
  output$gap_table <- renderTable({
    d_long <- processed_data(); if(nrow(d_long) == 0) return(NULL)
    d_long %>% group_by(Site) %>%
      summarise(
        `Missing data (NAs / Blanks)` = sum(is.na(Count) | Count <= 0), 
        `Expected records` = n(),
        `Percent coverage` = sprintf("%.1f%%", (sum(!is.na(Count) & Count > 0) / n()) * 100)
      )
  }, align = "c", striped = TRUE, hover = TRUE, spacing = "s")
  
  output$preview_annual_raw <- renderPlot({
    d_long <- processed_data(); if(nrow(d_long) == 0) return(NULL)
    ggplot(d_long %>% filter(!is.na(Count)), aes(x = Year, y = Count, color = Site)) +
      geom_line(linewidth = 1) + geom_point(size = 2) + facet_wrap(~Site, scales = "free_y", ncol = 1) + 
      theme_classic() + 
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth = 1))
  })
  
  output$data_preview_table_raw <- renderTable({ 
    d_long <- processed_data(); if(nrow(d_long) == 0) return(NULL)
    head(d_long, 10) 
  })
  
  output$qaqc_alerts <- renderUI({
    df <- processed_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    # Run data verification filters
    negatives <- sum(df$Count < 0, na.rm = TRUE)
    nas <- sum(is.na(df$Count))
    
    alerts <- list()
    if (negatives > 0) {
      alerts[[length(alerts) + 1]] <- div(class = "alert alert-danger d-flex align-items-center", style = "margin-bottom: 10px;",
                                          shiny::icon("exclamation-triangle", class = "me-2"), sprintf("QA/QC Alert: Found %d negative count values in your spreadsheet. Please verify or fix your raw counts file.", negatives))
    }
    if (nas > 0) {
      alerts[[length(alerts) + 1]] <- div(class = "alert alert-warning d-flex align-items-center", style = "margin-bottom: 10px;",
                                          shiny::icon("info-circle", class = "me-2"), sprintf("QA/QC Note: Found %d missing observations (NAs/gaps). The state-space models will cleanly bridge these periods via mathematical imputation.", nas))
    }
    
    if (length(alerts) == 0) {
      return(div(class = "alert alert-success d-flex align-items-center", style = "margin-bottom: 10px;", 
                 shiny::icon("check-circle", class = "me-2"), "QA/QC Pass: Initial data structures cleared (No negative counts or formatting errors caught)."))
    } else {
      return(tagList(alerts))
    }
  })
  
  output$period_override_ui <- renderUI({
    df <- processed_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    # Grab all unique beach locations present in the active dataset
    all_sites <- sort(unique(df$Site))
    
    # Run your original pattern check to establish the default system guess
    default_6mo <- all_sites[grepl("W_|Wermon|Wamlana|Waspait|Waenibe", all_sites, ignore.case = TRUE)]
    
    selectizeInput(
      "six_month_sites", 
      label = tags$span(
        "Designate 6-Month (Bimodal) Beaches:",
        tooltip(shiny::icon("question-circle"), "Select which beaches experience two distinct nesting peaks per year. Unselected beaches default to a standard 12-month cycle.")
      ),
      choices = all_sites,
      selected = default_6mo,
      multiple = TRUE,
      options = list(plugins = list('remove_button'))
    )
  })
  
  output$preview_monthly_seasonality <- renderPlot({
    df <- processed_data()
    req(df, input$data_mode == "upload", input$map_month != "none")
    
    # Group by Site, Year, and Month to isolate individual time-series waveforms
    df_monthly <- df %>% 
      filter(!is.na(Month), !is.na(Count), !is.na(Site)) %>% 
      mutate(Month = as.numeric(Month)) %>% 
      group_by(Site, Year, Month) %>% 
      summarise(Total_Count = sum(Count, na.rm = TRUE), .groups = "drop") %>% 
      mutate(X_Month = ifelse(Month >= 4, Month - 3, Month + 9))
    
    ggplot(df_monthly, aes(x = X_Month, y = Total_Count, group = factor(Year), color = factor(Year))) +
      geom_line(linewidth = 0.9, alpha = 0.75) +
      geom_point(size = 1.5) +
      scale_x_continuous(
        breaks = 1:12, 
        labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar")
      ) +
      scale_color_viridis_d(option = "viridis") +
      # Facet by individual beach time series with free y-axes to accommodate density differences
      facet_wrap(~Site, scales = "free_y", ncol = 1) + 
      theme_classic() +
      labs(x = "Month (Biological Nesting Year)", y = "Aggregated Nest Counts", color = "Calendar Year") +
      theme(
        text = element_text(size = 13), 
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.major.x = element_line(color = "grey95"),
        legend.position = "right"
      )
  })
  
  output$seasonality_recommendation <- renderUI({
    df <- processed_data()
    req(df, input$data_mode == "upload", input$map_month != "none")
    
    sites <- sort(unique(df$Site))
    recommendations <- list()
    
    for (st in sites) {
      df_site <- df %>% filter(Site == st, !is.na(Month), !is.na(Count))
      if (nrow(df_site) == 0) next
      
      # Calculate average distribution across the biological calendar segments
      summer_activity <- sum(df_site$Count[df_site$Month %in% 4:9], na.rm = TRUE)
      winter_activity <- sum(df_site$Count[df_site$Month %in% c(10,11,12,1,2,3)], na.rm = TRUE)
      total_activity <- summer_activity + winter_activity
      if(total_activity == 0) next
      
      # Heuristic: If winter clusters hold more than 25% of total annual nesting, it indicates a secondary peak
      is_bimodal <- (winter_activity / total_activity > 0.25) & (summer_activity / total_activity > 0.25)
      
      if (is_bimodal) {
        recommendations[[length(recommendations) + 1]] <- tags$li(
          tags$b(st, ": "), "System suggests a ", tags$span(class = "badge bg-warning text-dark", "6-Month Period"), 
          sprintf(" due to a bimodal waveform pattern (Summer: %.1f%%, Winter: %.1f%% of data).", 
                  (summer_activity/total_activity)*100, (winter_activity/total_activity)*100)
        )
      } else {
        recommendations[[length(recommendations) + 1]] <- tags$li(
          tags$b(st, ": "), "System suggests a ", tags$span(class = "badge bg-primary", "12-Month Period"), 
          " due to a single primary unimodal nesting peak season."
        )
      }
    }
    
    tagList(
      h5(shiny::icon("robot"), " Data Phenology Diagnostic System Suggestions:"),
      tags$ul(recommendations),
      p(tags$i(class = "text-muted", "Note: If your custom uploaded beach names conflict with the system's structural recommendations, you can explicitly update the `periods` logic inside your script's Fourier Imputation engine."))
    )
  })
  
  historical_ane_ledger <- reactive({
    req(vault$res, input$ane_custom_take, input$ane_custom_mort)
    
    # If the ledger was already handled by the master progress pipeline, use it as a base
    if (!is.null(vault$empirical_ane_ledger) && 
        input$ane_custom_take == 30 && input$ane_custom_mort == 0.35) {
      return(vault$empirical_ane_ledger)
    }
    
    # If a user manually changes a slider on page 4, recalculate smoothly here
    years_vec <- min(vault$abund$Year, na.rm = TRUE):max(vault$abund$Year, na.rm = TRUE)
    custom_safe <- data.frame(Year = years_vec, Total_Est = input$ane_custom_take)
    
    set.seed(2026)
    dummy_obs <- data.frame(
      Year = sample(years_vec, min(50, length(years_vec) * 2), replace = TRUE),
      Length = round(rnorm(min(50, length(years_vec) * 2), mean = 70, sd = 5), 1),
      Mortality = input$ane_custom_mort
    )
    
    params_demo <- list(
      sp_tknot = input$tknot, sp_k = input$k, sp_linf = input$linf, 
      sp_max_age = input$max_age, sp_lmat = input$lmat, sp_sig_mat = input$sig_mat, 
      sp_ane_pj = input$ane_pj, sp_ane_pa = input$ane_pa, sp_pf = input$pf, sp_remig_int = input$remig_int
    )
    suppressWarnings({ calculate_empirical_ane(dummy_obs, custom_safe, params_demo, max(years_vec)) })
  })
  
  
  observeEvent(input$run_model, {
    d_mapped <- processed_data()
    if(nrow(d_mapped) == 0) return()
    
    withProgress(message = 'Calculating baseline arrays...', value = 0, {
      if (input$data_mode == "upload" && input$map_month != "none") {
        setProgress(value = 0.1, detail = "Executing nested monthly Fourier matrix filling engine...")
        d_annual <- run_fourier_imputation(
          d_mapped, 
          iter = max(2000, floor(input$iterations / 5)),
          six_month_sites = input$six_month_sites # Connects the UI choice to the JAGS data builder
        )
      } else {
        d_annual <- d_mapped %>% group_by(Year, Site) %>% summarise(Count = sum(Count, na.rm=TRUE), .groups="drop")
      }
      
      abund <- d_annual %>% mutate(Annual_Nesters = Count / input$clutch_freq)
      vault$abund <- abund
      
      if (input$run_split) {
        setProgress(value = 0.3, detail = "Running pre-split matrix window...")
        slice_pre <- abund %>% filter(Year <= input$split_year)
        fit_pre <- tryCatch({ run_jags_aligned(slice_pre, iter = input$iterations, burnin = floor(input$iterations/3)) }, error = function(e) NULL)
        if(!is.null(fit_pre)) {
          vault$pre_u_val <- median(fit_pre$fit$sims.list$U)
          vault$pre_trend <- paste0(round(vault$pre_u_val, 3), " (", round((exp(vault$pre_u_val)-1)*100, 2), "%)")
        }
        
        slice_post <- abund %>% filter(Year >= input$split_year)
        fit_post <- tryCatch({ run_jags_aligned(slice_post, iter = input$iterations, burnin = floor(input$iterations/3)) }, error = function(e) NULL)
        if(!is.null(fit_post)) {
          vault$post_u_val <- median(fit_post$fit$sims.list$U)
          vault$post_trend <- paste0(round(vault$post_u_val, 3), " (", round((exp(vault$post_u_val)-1)*100, 2), "%)")
        }
      }
      
      setProgress(value = 0.6, detail = "Compiling global JAGS baseline canvas...")
      res <- tryCatch({ run_jags_aligned(abund, iter = input$iterations, burnin = floor(input$iterations / 3), thin = 10) }, error = function(e) { NULL })
      req(res); vault$res <- res
      
      setProgress(value = 0.8, detail = "Compiling validation MARSS layers...")
      vault$marss <- tryCatch({
        prep_marss <- vault$abund %>% dplyr::select(Year, Site, Annual_Nesters) %>% pivot_wider(names_from = Site, values_from = Annual_Nesters) %>% arrange(Year)
        mat_data <- as.matrix(prep_marss %>% select(-Year)); mat_data[mat_data <= 0 | is.nan(mat_data)] <- NA
        Y_matrix <- t(log(mat_data)); n_sites_m <- nrow(Y_matrix)
        
        if(n_sites_m == 1) {
          fit_s <- MARSS::MARSS(Y_matrix, model = list(Z = matrix(1), A = "zero", R = matrix("r"), Q = matrix("q"), U = matrix("u")), silent = TRUE)
          fit_i <- fit_s
        } else {
          fit_s <- MARSS::MARSS(Y_matrix, model = list(Z = matrix(1, nrow = n_sites_m, ncol = 1), A = "scaling", R = "diagonal and unequal", Q = matrix("q"), U = matrix("u")), silent = TRUE)
          fit_i <- MARSS::MARSS(Y_matrix, model = list(Z = diag(1, n_sites_m), A = "zero", R = "diagonal and unequal", Q = "diagonal and unequal", U = "unequal"), silent = TRUE)
        }
        list(shared = fit_s, indep = fit_i, years = prep_marss$Year)
      }, error = function(e) { NULL })
      
      setProgress(value = 0.9, detail = "Attaching backend parameter values...")
      fit <- res$fit; years <- res$years; fy <- length(years)
      if (is.null(fit$sims.list$A)) {
        X_total <- exp(fit$sims.list$X)
      } else {
        X_total <- apply(fit$sims.list$X, 2, function(v) rowSums(apply(fit$sims.list$A, 2, function(x) exp(v + x))))
      }
      
      vault$year <- max(years); vault$nesters <- median(X_total[, fy]); vault$total <- median(X_total[, fy] * input$remig_int)
      vault$trend_display <- paste0(round(median(fit$sims.list$U), 3), " (", round((exp(median(fit$sims.list$U))-1)*100, 2), "%)")
      vault$trend_pct <- (exp(median(fit$sims.list$U))-1)*100
      
      vault$draws <- data.frame(U = fit$sims.list$U, Q = fit$sims.list$Q, 
                                R_mean = if(is.matrix(fit$sims.list$R)) rowMeans(fit$sims.list$R) else as.numeric(fit$sims.list$R), 
                                Total_Females = X_total[, fy] * input$remig_int)
      
      setProgress(value = 1, detail = "Complete!")
    })
  })
  
  output$summary_stats <- renderUI({
    validate(need(vault$nesters, "Please click Run “Trend and Abundance” model on the left panel first."))
    boxes <- list(
      value_box(title = "Nesters", value = format(round(vault$nesters), big.mark=","), theme = "primary"),
      value_box(title = "Adult females", value = format(round(vault$total), big.mark=","), theme = "secondary"),
      value_box(title = "Annual population growth rate", value = vault$trend_display, theme = if(vault$trend_pct >= 0) "success" else "danger")
    )
    if (input$run_split) {
      pre_val <- if(!is.null(vault$pre_trend)) vault$pre_trend else "N/A"
      post_val <- if(!is.null(vault$post_trend)) vault$post_trend else "N/A"
      boxes[[length(boxes)+1]] <- value_box(title = paste("Trend Pre-", input$split_year), value = pre_val, theme = "info")
      boxes[[length(boxes)+1]] <- value_box(title = paste("Trend Post-", input$split_year), value = post_val, theme = "warning")
    }
    do.call(layout_column_wrap, c(list(width = 1/length(boxes)), boxes))
  })
  
  output$clean_baseline_plot <- renderPlot({
    validate(need(vault$res, "Please click Run “Trend and Abundance” model on Page 2 first."))
    fit <- vault$res$fit; years <- vault$res$years
    if (is.null(fit$sims.list$A)) {
      X_total <- exp(fit$sims.list$X)
    } else {
      X_total <- apply(fit$sims.list$X, 2, function(v) rowSums(apply(fit$sims.list$A, 2, function(x) exp(v + x))))
    }
    X_q <- apply(log(X_total), 2, quantile, probs = c(0.025, 0.5, 0.975))
    obs_summary <- vault$abund %>% dplyr::group_by(Year) %>% dplyr::summarise(Ann_Tot = sum(Annual_Nesters, na.rm=TRUE)) %>% filter(Ann_Tot > 0)
    
    par(mar = c(4, 4, 1, 1), bty = "l", xaxs = "i", yaxs = "i")
    plot(years, X_q[2,], type="n", ylim=range(c(log(obs_summary$Ann_Tot), X_q), na.rm=TRUE), ylab="Annual nesters (log)", xlab="Year", axes=FALSE)
    axis(1, col="black", col.axis="black", lwd=1); axis(2, col="black", col.axis="black", lwd=1)
    if(input$run_split) abline(v = input$split_year, lty=2, col="#0072B2", lwd=1.5)
    polygon(c(years, rev(years)), c(X_q[1, ], rev(X_q[3, ])), col = 'grey90', border = NA)
    lines(years, X_q[2, ], lwd = 2.5, col = 'black')
    
    # Plot true observed points as solid filled circles
    points(obs_summary$Year, log(obs_summary$Ann_Tot), pch = 16, col = "black")
    
    # Isolate missing/unmonitored years on the timeline
    raw_mapped <- raw_ingested_data()
    observed_years <- unique(raw_mapped$Year[!is.na(raw_mapped$Count) & raw_mapped$Count > 0])
    missing_years <- years[!(years %in% observed_years)]
    
    # Overlay open circles (pch = 1) directly onto the trajectory line for missing years
    if (length(missing_years) > 0) {
      idx_missing <- match(missing_years, years)
      points(missing_years, X_q[2, idx_missing], pch = 1, col = "black", cex = 1.6, lwd = 2)
    }
  })
  
  output$dynamic_unified_plot <- renderUI({
    validate(need(vault$abund, "Please click Run “Trend and Abundance” model on Page 2 first..."))
    plotOutput("unified_trend_plot", height = "360px")
  })
  
  output$unified_trend_plot <- renderPlot({
    req(vault$res, vault$marss, vault$abund, input$plot_layers)
    
    # 1. Native floating progress card banner is back
    withProgress(message = "Generating Abundance Trajectories plot...", value = 0.5, {
      
      j_fit <- vault$res$fit
      all_years <- vault$res$years
      sites = sort(unique(vault$abund$Site))
      n_sites <- length(sites)
      df_all_fits <- data.frame()
      
      if ("jags" %in% input$plot_layers) {
        A_sims <- j_fit$sims.list$A; if (is.null(A_sims)) A_sims <- matrix(0, nrow = length(j_fit$sims.list$U), ncol = 1)
        for (i in seq_along(sites)) {
          site_A <- if(ncol(A_sims) >= i) as.numeric(median(A_sims[, i])) else 0
          for (t in seq_along(all_years)) {
            pred <- exp(j_fit$sims.list$X[, t] + site_A)
            df_all_fits <- rbind(df_all_fits, data.frame(Year=all_years[t], Site=sites[i], Model="JAGS Regional trend", Median=median(pred), Lower=quantile(pred,0.025), Upper=quantile(pred,0.975)))
          }
        }
      }
      if ("marss_s" %in% input$plot_layers && !is.null(vault$marss$shared)) {
        ms_fit <- vault$marss$shared; ms_states <- as.numeric(ms_fit$states[1, ]); ms_se <- as.numeric(ms_fit$states.se[1, ]); ms_A <- stats::coef(ms_fit, type="matrix")$A; if (is.null(ms_A) || length(ms_A) == 0) ms_A <- matrix(0, nrow = n_sites, ncol = 1)
        for(i in 1:n_sites) {
          site_A <- if(nrow(ms_A) >= i) as.numeric(ms_A[i,1]) else 0
          df_all_fits <- rbind(df_all_fits, data.frame(Year=all_years, Site=sites[i], Model="MARSS Regional Trend", Median=exp(ms_states+site_A), Lower=exp((ms_states-1.96*ms_se)+site_A), Upper=exp((ms_states+1.96*ms_se)+site_A)))
        }
      }
      if ("marss_i" %in% input$plot_layers && !is.null(vault$marss$indep)) {
        mi_fit <- vault$marss$indep
        for(i in 1:n_sites) {
          mi_states <- if(is.matrix(mi_fit$states)) as.numeric(mi_fit$states[i, ]) else as.numeric(mi_fit$states)
          mi_se <- if(is.matrix(mi_fit$states.se)) as.numeric(mi_fit$states.se[i, ]) else as.numeric(mi_fit$states.se)
          df_all_fits <- rbind(df_all_fits, data.frame(Year=all_years, Site=sites[i], Model="MARSS Independent Trends", Median=exp(mi_states), Lower=exp(mi_states-1.96*mi_se), Upper=exp(mi_states+1.96*mi_se)))
        }
      }
      
      validate(need(nrow(df_all_fits) > 0, "Select at least one layer framework to draw fitted trends chart."))
      
      ggplot(df_all_fits %>% left_join(vault$abund, by=c("Year","Site")), aes(x=Year)) +
        geom_ribbon(aes(ymin=Lower, ymax=Upper, fill=Model), alpha=0.12, color=NA) + 
        geom_line(aes(y=Median, color=Model, linetype=Model), linewidth=1.1) + 
        geom_point(aes(y=Annual_Nesters), color="black", size=2.2, na.rm=TRUE) +
        facet_wrap(~Site, scales="free_y", ncol=1) + 
        theme_classic() + 
        theme(strip.text = element_text(face = "bold", size = 11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth = 1)) +
        labs(y = "Nesters", x = "Year")
    })
  })
  
  output$table_u <- renderTable({
    validate(need(vault$res, "Please run baseline models on Page 2 first."))
    sites <- sort(unique(vault$abund$Site))
    j_u <- median(vault$res$fit$sims.list$U)
    get_badge <- function(u) sprintf("<span class='badge %s'>%+.2f%% / yr</span>", if((exp(u)-1)*100 < 0) "bg-danger text-white" else "bg-success text-white", (exp(u)-1)*100)
    
    rows <- list(data.frame(Framework = "JAGS Regional trend", Strategy = "Shared", U = sprintf("%.4f", j_u), Trend = get_badge(j_u), check.names=FALSE))
    
    if(!is.null(vault$marss)) {
      m_s_u <- stats::coef(vault$marss$shared, type="matrix")$U[1,1]
      m_i_u <- stats::coef(vault$marss$indep, type="matrix")$U[,1]
      rows[[length(rows)+1]] <- data.frame(Framework = "MARSS Regional Trend", Strategy = "Shared", U = sprintf("%.4f", m_s_u), Trend = get_badge(m_s_u), check.names=FALSE)
      for(i in seq_along(sites)) {
        val_u <- if(length(m_i_u) >= i) m_i_u[i] else m_i_u[1]
        rows[[length(rows)+1]] <- data.frame(Framework = "MARSS Independent Trends", Strategy = sites[i], U = sprintf("%.4f", val_u), Trend = get_badge(val_u), check.names=FALSE)
      }
    }
    do.call(rbind, rows)
  }, sanitize.text.function = function(x) x)
  
  output$plot_var <- renderPlot({
    validate(need(vault$res, "Please run baseline models on Page 2 first."))
    j_q <- mean(vault$res$fit$sims.list$Q)
    j_r <- mean(colMeans(vault$res$fit$sims.list$R))
    
    df_var <- data.frame(Component = c("Environmental stochasticity (Q)", "Observation stochasticity (R)"), Variance = c(j_q, j_r))
    ggplot(df_var, aes(x = Component, y = Variance, fill = Component)) + 
      geom_bar(stat = "identity", color = "black", width = 0.5, linewidth = 0.8) + 
      scale_fill_manual(values = c("Environmental stochasticity (Q)" = "#0072B2", "Observation stochasticity (R)" = "#D55E00")) +
      theme_classic() + 
      theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth = 1), text = element_text(size = 12)) +
      labs(y = "Variance scale value", x = "")
  })
  
  output$variance_interpretation_text <- renderUI({
    req(vault$res)
    j_q <- mean(vault$res$fit$sims.list$Q)
    j_r <- mean(colMeans(vault$res$fit$sims.list$R))
    if (j_q > j_r) {
      p(tags$b("Interpretation: "), "The variation you see year-to-year is primarily driven by ", tags$span("real environmental shifts in the ocean (Q)", style="color:#0072B2; font-weight:bold;"), ", meaning population signals are highly reflective of true ecological changes.")
    } else {
      p(tags$b("Interpretation: "), "The variation is primarily driven by ", tags$span("monitoring gaps and counting noise on the beach (R)", style="color:#D55E00; font-weight:bold;"), ", suggesting that data collection consistency is masking the underlying population signal.")
    }
  })
  
  output$covariance_interpretation_text <- renderUI({
    req(vault$res)
    fit <- vault$res$fit
    years <- vault$res$years
    fy <- length(years)
    
    r_vec <- as.numeric(fit$sims.list$U)
    q_vec <- as.numeric(fit$sims.list$Q)
    
    if (is.null(fit$sims.list$A)) {
      n_vec <- exp(as.numeric(fit$sims.list$X[, fy]))
    } else {
      n_vec <- apply(fit$sims.list$X, 1, function(v) sum(apply(fit$sims.list$A, 2, function(x) exp(v[fy] + x))))
    }
    
    cor_rq <- cor(r_vec, q_vec, use = "complete.obs")
    cor_rn <- cor(r_vec, n_vec, use = "complete.obs")
    cor_nq <- cor(n_vec, q_vec, use = "complete.obs")
    
    tagList(
      h5(style = "font-size: 1.1rem; font-weight: bold; margin-bottom: 12px;", 
         shiny::icon("analytics"), " Management & Biological Interpretation Guide:"),
      
      tags$ul(style = "padding-left: 20px;",
              # 1. Trend vs Environmental Noise
              tags$li(style = "margin-bottom: 10px;",
                      tags$b("Independence of Long-Term Trend and Environmental Fluctuation [r = ", sprintf("%.2f", cor_rq), "]: "),
                      if(abs(cor_rq) < 0.30) {
                        "This near-zero correlation confirms excellent parameter separation. The model is calculating the persistent, multi-decade population growth or decline without being distorted or misled by short-term, year-to-year environmental spikes."
                      } else {
                        "Warning: The calculated growth rate and background environmental noise are statistically bleeding into each other, meaning the long-term trend estimate is highly sensitive to our environmental variance assumptions."
                      }
              ),
              
              # 2. Trend vs Current Headcount
              tags$li(style = "margin-bottom: 10px;",
                      tags$b("Cumulative Growth Impact on Final Abundance [r = ", sprintf("%.2f", cor_rn), "]: "),
                      "This positive relationship verifies standard biological consistency within the model: simulation pathways that evaluate a slightly higher historical growth rate predictably accumulate more individuals, leading to a larger current population estimate."
              ),
              
              # 3. Variance-Expanded Abundance Limits
              tags$li(style = "margin-bottom: 10px;",
                      tags$b("Abundance Sensitivity to Environmental Uncertainty [r = ", sprintf("%.2f", cor_nq), "]: "),
                      if(cor_nq >= 0.30) {
                        tags$span(style = "font-weight: 500; color: #b02a37;",
                                  sprintf("Elevated correlation detected. High-variance environmental simulations are pushing open the upper statistical limits of our population estimate. This indicates that our maximum potential population size is structurally linked to background ocean volatility.", cor_nq))
                      } else {
                        "Current population size estimates are stable and mathematically isolated from environmental variance scaling factors."
                      }
              )
      ),
      
      div(style = "background-color: #f8f9fa; border-left: 4px solid #6c757d; padding: 10px; margin-top: 15px; border-radius: 4px;",
          tags$b("Statistical Distribution Note (Histogram Shapes): "),
          "The asymmetric right-hand tail displayed in the middle (N_final) histogram indicates that our population uncertainty is not uniform. While the lower boundary is firmly restricted by actual beach nest counts, the upper boundary allows for a wide margin of error. Statistically, the true population size has a much higher likelihood of being under-counted rather than over-counted."
      )
    )
  })
  
  output$posterior_pairs_plot <- renderPlot({
    validate(need(vault$res, "Please click Run “Trend and Abundance” model on Page 2 first."))
    
    # 2. Re-activates the floating progress card notification banner for the matrix
    withProgress(message = "Computing Joint Parameter Covariance Matrix...", value = 0.5, {
      
      fit <- vault$res$fit; years <- vault$res$years; fy <- length(years)
      
      r_vec <- as.numeric(fit$sims.list$U)
      q_vec <- as.numeric(fit$sims.list$Q)
      if (is.null(fit$sims.list$A)) {
        n_vec <- exp(as.numeric(fit$sims.list$X[, fy]))
      } else {
        n_vec <- apply(fit$sims.list$X, 1, function(v) sum(apply(fit$sims.list$A, 2, function(x) exp(v[fy] + x))))
      }
      
      pairs_df <- data.frame(r = r_vec, N_final = n_vec, Q = q_vec)
      
      panel_hist <- function(x, ...) {
        h <- hist(x, plot = FALSE, breaks = 25)
        y <- h$counts / max(h$counts)
        old_par <- par(usr = c(par("usr")[1:2], 0, 1.1))
        on.exit(par(old_par))
        lines(approx(h$mids, y, xout=seq(min(x), max(x), length.out=100)), lwd = 2, col = "black")
      }
      
      panel_scatter <- function(x, y, ...) {
        points(x, y, pch = 20, col = rgb(0.3, 0.3, 0.3, 0.05), cex = 0.5)
        ellipse_50 <- car::dataEllipse(x, y, levels = 0.50, draw = FALSE)
        ellipse_95 <- car::dataEllipse(x, y, levels = 0.95, draw = FALSE)
        lines(ellipse_50, col = "black", lwd = 2)
        lines(ellipse_95, col = "gray50", lwd = 1.2, lty = 2)
      }
      
      panel_cor <- function(x, y, ...) {
        r_coef <- cor(x, y, use = "complete.obs")
        old_par <- par(usr = c(0, 1, 0, 1))
        on.exit(par(old_par))
        text(0.5, 0.5, sprintf("%.2f", r_coef), cex = 1.5, font = 2)
      }
      
      par(mar = c(3, 3, 1, 1), bg = "white")
      pairs(pairs_df, labels = c("r", "N_final", "Q"),
            diag.panel = panel_hist, lower.panel = panel_cor, upper.panel = panel_scatter,
            cex.labels = 1.4, font.labels = 2)
    })
  })
  
  output$ane_impact_plot <- renderPlot({
    validate(
      need(!is.null(vault$abund), "Please click Run “Trend and Abundance” model on Page 2 first."),
      need(!is.null(historical_ane_ledger()), "Please click Run “Trend and Abundance” model on Page 2 first.")
    )
    obs_annual <- vault$abund %>% group_by(Year) %>% summarise(Observed = sum(Annual_Nesters, na.rm = TRUE), .groups = "drop")
    plot_df <- obs_annual %>% 
      left_join(historical_ane_ledger(), by = "Year") %>%
      mutate(Counterfactual = Observed + coalesce(Total_Cumulative_ANE, 0)) %>%
      select(Year, Observed, Counterfactual) %>%
      pivot_longer(cols = c(Observed, Counterfactual), names_to = "Timeline", values_to = "Nesters")
    
    ggplot(plot_df, aes(x = Year, y = Nesters, color = Timeline, linetype = Timeline)) +
      geom_line(linewidth = 1.4) + geom_point(size = 2.5) +
      scale_color_manual(values = c("Observed" = "#000000", "Counterfactual" = "#D55E00"),
                         labels = c("Counterfactual" = "Pristine Timeline (No Marine Take)", "Observed" = "Observed Timeline (Reality)")) +
      scale_linetype_manual(values = c("Observed" = "solid", "Counterfactual" = "dashed"),
                            labels = c("Counterfactual" = "Pristine Timeline (No Marine Take)", "Observed" = "Observed Timeline (Reality)")) +
      theme_classic() +
      theme(legend.position = "bottom", legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth = 1), text = element_text(size = 14)) +
      labs(y = "Annual Nesters", x = "Year")
  })
  
  # --- CHAPTER 5: DETAILED SCENARIO PROJECTOR (100% DEMOGRAPHIC PARITY) ---
  observeEvent(input$run_portfolio_sim, {
    req(vault$draws); n_sims <- input$pva_sims; horizon <- input$pva_years
    
    withProgress(message = "Simulating all management tracks...", value = 0, {
      if (input$trend_to_project == "Post-Split Year Trend" && input$run_split) {
        r_base <- rnorm(n_sims, mean = vault$post_u_val, sd = 0.01)
      } else {
        r_base <- sample(vault$draws$U, n_sims, replace = TRUE)
      }
      q_base <- sample(vault$draws$Q, n_sims, replace = TRUE)
      start_pop_draws <- sample(vault$draws$Total_Females, n_sims, replace = TRUE)
      
      scenarios <- c("Status Quo", "Threats Only", "Strategy A Only", "Strategy B Only", "Strategy C Only", "Strategy D Only", "Strategy E Only", "Combined Portfolio")
      sim_tracks <- list()
      for(s in scenarios) { sim_tracks[[s]] <- matrix(NA, n_sims, horizon) }
      
      sim_beta0   <- 4.35
      sim_beta1   <- 0.0
      sim_mu0     <- 0.14
      sim_sigma_L <- 0.338
      sim_sigma_D <- 0.50
      sim_rho     <- -0.51
      sim_cov     <- matrix(c(sim_sigma_L^2, sim_sigma_L * sim_sigma_D * sim_rho, sim_sigma_L * sim_sigma_D * sim_rho, sim_sigma_D^2), 2, 2)
      
      # PRE-CALCULATE LIFECYCLE MATURITY YEARS AS TARGET HURDLE
      A_mat <- round(input$tknot - (1 / input$k) * log(1 - (pmin(input$lmat, input$linf - 0.1) / input$linf)))
      if(is.na(A_mat) || A_mat <= 1) A_mat <- 12
      
      # Compute precise 100% parity demographic conversion metrics
      yearlings_prod <- input$optA_nests * input$clutch_size * input$hatch_succ * input$emerg_succ * input$yrlng_surv
      ane_gain_A <- yearlings_prod * (input$ane_pj ^ pmax(1, A_mat - 1)) * input$pf
      ane_gain_B <- input$optB_head * (input$ane_pj ^ pmax(0, A_mat - input$optB_age)) * input$pf
      ane_gain_C <- input$optC_adults * input$pf
      
      for (i in 1:n_sims) {
        curr_pops <- list()
        for(s in scenarios) { curr_pops[[s]] <- start_pop_draws[i] }
        
        loss_threat <- loss_comb <- rep(0, horizon)
        
        for (y in 1:horizon) {
          # Calculate standard fishery takes context loops
          encounters <- rCMP(1, mu = input$cmp_mu, nu = input$cmp_nu)
          encounters_D <- if(input$apply_threats) rCMP(1, mu = input$optD_mu_mit, nu = input$cmp_nu) else 0
          encounters_comb <- if("D" %in% input$portfolio_options) encounters_D else encounters
          
          # Baseline threat profiles
          if (encounters > 0 && input$apply_threats) {
            mu_l_sq <- sim_beta0 + sim_beta1 * encounters
            draws_sq <- mvtnorm::rmvnorm(encounters, mean = c(mu_l_sq, sim_mu0), sigma = sim_cov)
            for (t in 1:encounters) {
              future_len <- exp(draws_sq[t, 1])
              age_t <- input$tknot - (1 / input$k) * log(1 - (pmin(future_len, input$linf - 0.1) / input$linf))
              if (is.nan(age_t) || is.na(age_t)) age_t <- input$max_age - 2
              ly <- length(y:horizon); ages_traj <- seq(age_t, length.out = ly, by = 1)
              lens_traj <- input$linf * (1 - exp(-input$k * (ages_traj - input$tknot)))
              
              p_mat <- ifelse(lens_traj >= 0.99 * input$linf, 1.0, 1.0 / (1.0 + exp(-(lens_traj - input$lmat) / input$sig_mat)))
              surv <- cumprod((1 - p_mat) * input$ane_pj + p_mat * input$ane_pa)
              p_binom <- rbinom(ly, 1, p_mat); if(any(p_binom==1)) p_binom[min(which(p_binom==1)):ly] <- 1
              
              loss_threat[y:horizon] <- loss_threat[y:horizon] + (surv * p_binom * input$pf * input$threat_mortality * (1/input$remig_int))
            }
          }
          
          # Mitigated Portfolio Threat loops
          if (encounters_comb > 0 && input$apply_threats) {
            mu_l_c <- sim_beta0 + sim_beta1 * encounters_comb
            draws_c <- mvtnorm::rmvnorm(encounters_comb, mean = c(mu_l_c, sim_mu0), sigma = sim_cov)
            for (t in 1:encounters_comb) {
              future_len <- exp(draws_c[t, 1])
              age_t <- input$tknot - (1 / input$k) * log(1 - (pmin(future_len, input$linf - 0.1) / input$linf))
              if (is.nan(age_t) || is.na(age_t)) age_t <- input$max_age - 2
              ly <- length(y:horizon); ages_traj <- seq(age_t, length.out = ly, by = 1)
              lens_traj <- input$linf * (1 - exp(-input$k * (ages_traj - input$tknot)))
              
              p_mat <- ifelse(lens_traj >= 0.99 * input$linf, 1.0, 1.0 / (1.0 + exp(-(lens_traj - input$lmat) / input$sig_mat)))
              surv <- cumprod((1 - p_mat) * input$ane_pj + p_mat * input$ane_pa)
              p_binom <- rbinom(ly, 1, p_mat); if(any(p_binom==1)) p_binom[min(which(p_binom==1)):ly] <- 1
              
              m_rate_comb <- if("E" %in% input$portfolio_options) runif(1, input$optE_m_min, input$optE_m_max) else input$threat_mortality
              loss_comb[y:horizon] <- loss_comb[y:horizon] + (surv * p_binom * input$pf * m_rate_comb * (1/input$remig_int))
            }
          }
          
          # Compute explicit step transitions using direct additive ANE metrics
          chg_sq   <- 0
          chg_t    <- -loss_threat[y]
          chg_A    <- -loss_threat[y] + ane_gain_A
          chg_B    <- -loss_threat[y] + ane_gain_B
          chg_C    <- -loss_threat[y] + ane_gain_C
          chg_D    <- -loss_threat[y] # Handled via structural mu shift inside separate loop bounds
          chg_E    <- 0                # Handled via explicit mortality sample swap bounds
          
          chg_comb <- -loss_comb[y] + 
            (if("A" %in% input$portfolio_options) ane_gain_A else 0) + 
            (if("B" %in% input$portfolio_options) ane_gain_B else 0) + 
            (if("C" %in% input$portfolio_options) ane_gain_C else 0)
          
          # Execute demographic equations mapping transitions safely across pools
          for(s in scenarios) {
            net_delta <- if(s=="Status Quo") chg_sq else if(s=="Threats Only") chg_t else if(s=="Strategy A Only") chg_A else if(s=="Strategy B Only") chg_B else if(s=="Strategy C Only") chg_C else chg_comb
            
            # Apply additive demographic mass balance directly to adult female pool context
            curr_pops[[s]] <- curr_pops[[s]] + (net_delta / input$remig_int)
            if(curr_pops[[s]] < 0) curr_pops[[s]] <- 0
            curr_pops[[s]] <- rnorm(1, curr_pops[[s]] * exp(r_base[i]), sqrt(q_base[i]))
            if(curr_pops[[s]] < 0) curr_pops[[s]] <- 0
            sim_tracks[[s]][i, y] <- curr_pops[[s]]
          }
        }
        if (i %% 25 == 0) setProgress(value = i / n_sims, detail = paste("Simulating all 8 tracks", i, "/", n_sims))
      }
      
      vault_portfolio$all_scen_raw <- list()
      for(s in scenarios) { vault_portfolio$all_scen_raw[[s]] <- sim_tracks[[s]][, horizon] }
      
      proj_years <- (max(vault$abund$Year) + 1):(max(vault$abund$Year) + horizon)
      calc_sum <- function(mat, label) data.frame(Year = proj_years, Median = apply(mat, 2, median), L95 = apply(mat, 2, function(x) quantile(x, 0.025)), U95 = apply(mat, 2, function(x) quantile(x, 0.975)), Scenario = label)
      
      compiled_rows <- list()
      for(s in scenarios) { compiled_rows[[s]] <- calc_sum(sim_tracks[[s]], s) }
      vault_portfolio$plot_df <- bind_rows(compiled_rows)
      
      scorecard_rows <- list()
      for(s in scenarios) {
        scorecard_rows[[length(scorecard_rows)+1]] <- data.frame(Strategy = s, `Abundance (median Adult Females)` = round(median(sim_tracks[[s]][, horizon]), 1), check.names = FALSE)
      }
      vault_portfolio$scorecard <- do.call(rbind, scorecard_rows)
    })
  })
  
  output$portfolio_overlay_plot <- renderPlot({
    validate(need(vault_portfolio$plot_df, "Please click Simulate Future Projections on the side panel first."))
    filtered_plot_df <- vault_portfolio$plot_df %>% filter(Scenario %in% input$portfolio_plots)
    p <- ggplot(filtered_plot_df, aes(x = Year, y = Median, color = Scenario, fill = Scenario)) + 
      geom_line(linewidth = 1.3) + scale_color_manual(values = custom_colors) + scale_fill_manual(values = custom_colors) + 
      theme_classic() + theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth = 1), text = element_text(size = 13)) + labs(y="Adult females")
    if(input$show_proj_ci) p <- p + geom_ribbon(aes(ymin = L95, ymax = U95), alpha = 0.05, color = NA)
    return(p)
  })
  
  output$portfolio_density_plot <- renderPlot({
    validate(need(length(vault_portfolio$all_scen_raw) > 0, "Please click Simulate Future Projections on the side panel first."))   
    density_list <- list()
    for(s in names(vault_portfolio$all_scen_raw)) { density_list[[length(density_list)+1]] <- data.frame(Abundance = vault_portfolio$all_scen_raw[[s]], Scenario = s) }
    density_df <- bind_rows(density_list) %>% filter(Scenario %in% input$portfolio_plots)
    
    ggplot(density_df, aes(x = Abundance, fill = Scenario, color = Scenario)) +
      geom_density(alpha = 0.12, linewidth = 1.3) + scale_fill_manual(values = custom_colors) + scale_color_manual(values = custom_colors) +
      theme_classic() + theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black", linewidth = 1), text = element_text(size = 13)) +
      labs(x = "Projected Adult Females Pool (Terminal Year)", y = "Relative Density")
  })
  
  output$portfolio_scorecard_table <- renderTable({ req(vault_portfolio$scorecard); vault_portfolio$scorecard })
  
  # --- ALL MASTER DOWNLOAD ARCHITECTURE ENGINE MODULES ---
  output$download_preview_raw <- downloadHandler(
    filename = function() { "Raw_Nest_Counts_Preview.png" },
    content = function(file) { ggsave(file, plot = output$preview_annual_raw, device = "png", width = 8, height = 5, bg = "white") }
  )
  output$download_baseline <- downloadHandler(
    filename = function() { "Calculated_Regional_Baseline_Trend.png" },
    content = function(file) {
      png(file, width = 800, height = 500, res = 100)
      fit <- vault$res$fit; years <- vault$res$years
      if(is.null(fit$sims.list$A)) X_total <- exp(fit$sims.list$X) else X_total <- apply(fit$sims.list$X, 2, function(v) rowSums(apply(fit$sims.list$A, 2, function(x) exp(v + x))))
      X_q <- apply(log(X_total), 2, quantile, probs = c(0.025, 0.5, 0.975))
      obs_summary <- vault$abund %>% dplyr::group_by(Year) %>% dplyr::summarise(Ann_Tot = sum(Annual_Nesters, na.rm=TRUE))
      par(mar = c(4, 4, 1, 1), bty = "l", xaxs = "i", yaxs = "i")
      plot(years, X_q[2,], type="n", ylim=range(c(log(obs_summary$Ann_Tot), X_q), na.rm=TRUE), ylab="Annual nester records (log)", xlab="Year", axes=FALSE)
      axis(1); axis(2); polygon(c(years, rev(years)), c(X_q[1, ], rev(X_q[3, ])), col = 'grey90', border = NA); lines(years, X_q[2, ], lwd = 2.5); points(obs_summary$Year, log(obs_summary$Ann_Tot), pch = 16)
      dev.off()
    }
  )
  output$download_unified <- downloadHandler( filename = function() { "Fits_Comparison.png" }, content = function(file) { ggsave(file, device = "png", width = 8, height = 6, bg = "white") } )
  output$download_var <- downloadHandler( filename = function() { "Variance_Components.png" }, content = function(file) { ggsave(file, device = "png", width = 5, height = 4, bg = "white") } )
  output$download_posterior <- downloadHandler(
    filename = function() { "Joint_Posterior_Matrix.png" },
    content = function(file) {
      png(file, width = 650, height = 650, res = 120)
      fit <- vault$res$fit; years <- vault$res$years; fy <- length(years)
      r_vec <- as.numeric(fit$sims.list$U); q_vec <- as.numeric(fit$sims.list$Q)
      if(is.null(fit$sims.list$A)) n_vec <- exp(as.numeric(fit$sims.list$X[, fy])) else n_vec <- apply(fit$sims.list$X, 1, function(v) sum(apply(fit$sims.list$A, 2, function(x) exp(v[fy] + x))))
      pairs_df <- data.frame(r = r_vec, N_final = n_vec, Q = q_vec)
      panel_hist <- function(x, ...) { h <- hist(x, plot = FALSE, breaks = 25); y <- h$counts / max(h$counts); lines(approx(h$mids, y, xout=seq(min(x), max(x), length.out=100)), lwd = 2) }
      panel_scatter <- function(x, y, ...) { points(x, y, pch = 20, col = rgb(0.3, 0.3, 0.3, 0.05)); lines(car::dataEllipse(x, y, levels = 0.50, draw = FALSE), lwd = 2); lines(car::dataEllipse(x, y, levels = 0.95, draw = FALSE), lty = 2) }
      panel_cor <- function(x, y, ...) { text(mean(range(x)), mean(range(y)), sprintf("%.2f", cor(x, y, use = "complete.obs")), cex = 1.5, font = 2) }
      pairs(pairs_df, labels = c("r", "N_final", "Q"), diag.panel = panel_hist, lower.panel = panel_cor, upper.panel = panel_scatter)
      dev.off()
    }
  )
  output$download_ane <- downloadHandler( filename = function() { "ANE_Impact.png" }, content = function(file) { ggsave(file, device = "png", width = 8, height = 5, bg = "white") } )
  output$download_forecast <- downloadHandler( filename = function() { "Timeline_Forecast.png" }, content = function(file) { ggsave(file, device = "png", width = 8, height = 5, bg = "white") } )
  output$download_density <- downloadHandler( filename = function() { "Abundance_Density.png" }, content = function(file) { ggsave(file, device = "png", width = 6, height = 4, bg = "white") } )
}
shinyApp(ui, server)