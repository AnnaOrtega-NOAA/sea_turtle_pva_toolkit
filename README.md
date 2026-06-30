# Sea Turtle Population Viability Analysis (PVA) Toolkit

Welcome to the **Sea Turtle Population Viability Analysis (PVA) Toolkit**—a comprehensive R Shiny workspace designed for demographic accounting, state-space population trend estimation, historical take auditing, and stochastic recovery forecasting. 

This platform bridges the gap between raw beach monitoring metrics and robust, publishable population assessments, supporting regional turtle conservation programs (e.g., North Pacific Loggerheads, western Pacific Leatherbacks).

---

## System Dependencies & Requirements

To run this application successfully, you need the R language environment along with a localized compiled software library for Bayesian Gibbs sampling.

### 1. External System Requirements
The state-space assessment and monthly imputation engines rely on **JAGS (Just Another Gibbs Sampler)**. You must install this external binary on your machine before compiling the R packages:
* **Windows/macOS/Linux:** Download and install the latest stable version of JAGS from [SourceForge](https://sourceforge.net/projects/mcmc-jags/).

### 2. Required R Packages
Ensure your local R library contains the following active dependencies. You can install them collectively by executing this block in your R console:

```R
install.packages(c(
  "shiny",    # Interactive framework architecture
  "bslib",    # Modern Yeti layout UI theme structures
  "dplyr",    # Data manipulation & pipe processing
  "tidyr",    # Wide-to-long form matrix operations
  "ggplot2",  # Trajectory and density canvas plotting
  "purrr",    # Structural loop vectorization
  "mvtnorm",  # Multivariate normal distributions for somatic models
  "truncnorm",# Truncated normal demographic wrappers
  "jagsUI",   # Interface to JAGS via parallel/R pipelines
  "MARSS"     # Multivariate Autoregressive State-Space modeling engine
))
```

---

## App Architecture & Core Chapters

The application is structured sequentially into five core chapters, tracking a rigorous data-to-decision pipeline:

### 1. Data Ingestion Workspace
* **Functionality:** Ingests raw nesting data via built-in simulation demo channels or user-supplied spreadsheet `.csv` uploads. Supports both **Long Form** and **Wide Form** data structures.
* **Automated QA/QC Scanning:** Instantly evaluates data health upon upload, flagging illegal negative counts, identifying outliers, and mapping missing observation gaps (`NA`s) before modeling.
* **Phenology Diagnostics & Overrides:** Visualizes monthly nesting waveforms (shifted to an April–March biological calendar) to detect unimodal (12-month) vs. bimodal (6-month) nesting peaks. Includes an automated heuristic recommendation engine to suggest cycle periodicities, paired with manual user override controls.

### 2. Abundance and Trend Assessment
* **Fourier Monthly Imputation:** If data includes a monthly sub-component with standard data gaps, an automated monthly JAGS Fourier loop fills missing matrices while accommodating the distinct reproductive periodicities assigned in Chapter 1.
* **JAGS State-Space Trend Engine:** Fits a joint log-linear trend across beaches. It separates true biological process variance (Q) from observation counting errors (R) and supports structural break testing (e.g., assessing if a shift in monitoring protocols in a specific year altered calculated growth coefficients).
* **Missing Data Visibility:** Explicitly maps completely unmonitored years as open circles on the final trajectory line so stakeholders can see exactly where the model relies on statistical imputation.

### 3. Beach Diagnostics Lab
* **Dual-Framework Validation:** Superimposes your unified Bayesian JAGS results against independent maximum likelihood frameworks (MARSS Shared vs. Independent Beach trend models).
* **Variance Scaler Isolator:** Explicitly visualizes the exact ratio of environmental noise vs observation error, giving managers objective visibility into whether counts are clouding underlying signals.
* **Joint Parameter Covariance Matrix:** Generates an interactive Pairs Plot showing the model's posterior distributions, Pearson Correlation logs, and joint 50%/95% confidence data ellipses to confirm complete parameter separation. Includes an automated, non-technical translation guide to help stakeholders understand variance-driven abundance ceilings and asymmetrical uncertainty buffers.

### 4. Retrospective Removal Ledger (ANEs)
* **Adult Nester Equivalents (ANE):** Converts raw historic counts of turtle interactions (e.g., longline bycatch) into reproductive adult nester equivalents (ANEs).
* **Demographic Hindcasting:** Projects takes backward using custom von Bertalanffy somatic growth curves, logistic sexual maturity milestones, remigration intervals, and stage-specific multi-year survival coefficients to reconstruct a "Pristine counterfactual timeline" showing where population levels would stand today if historical takes had not happened.

### 5. Future Projections Sandbox
* **Stochastic Forecaster:** Projects population pools up to 100 years into the future by sampling directly from the joint uncertainty profiles computed in the baseline diagnostics.
* **Conway-Maxwell-Poisson (CMP) Threat Profiles:** Injects customized future fishery hazard events modeling over- or under-dispersed fleet encounter dynamics.
* **Strategy Portfolio Builder:** Simulates standalone or combined recovery portfolios across five management levers:
  * *Strategy A:* Localized nest protection enclosures (clutch size, hatchling emergence multipliers).
  * *Strategy B:* Headstarting programs (custom release age and count).
  * *Strategy C:* Adult female poaching cessation targets.
  * *Strategy D:* Gear modifications to actively limit fleet interaction rates.
  * *Strategy E:* Advanced handling protocols to reduce post-release interaction mortality.

---

## Getting Started

To launch the app workspace locally from your R command line console, execute the following script:

```R
library(shiny)
# Path pointing to where your script is saved
runApp("C:/Users/YOUR_USER/PATH_TO_APP/app_demo.R")
```
