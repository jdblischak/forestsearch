# =============================================================================
# Execute Vignette Code with Synthetic Dataset
# =============================================================================

library(survival)
# library(randomizr) - implementing block_ra locally due to installation issues

# Local implementation of block_ra from randomizr package
# Performs block (stratified) random assignment
block_ra <- function(blocks, prob = 0.5) {
  blocks <- as.character(blocks)
  unique_blocks <- unique(blocks)
  assignment <- rep(NA, length(blocks))
  
  for (block in unique_blocks) {
    block_indices <- which(blocks == block)
    n_block <- length(block_indices)
    n_treat <- round(n_block * prob)
    
    # Randomly assign within block
    treat_indices <- sample(block_indices, n_treat)
    assignment[block_indices] <- 0
    assignment[treat_indices] <- 1
  }
  
  return(assignment)
}

# Try to load data.table, use base R if not available
has_dt <- requireNamespace("data.table", quietly = TRUE)
if (has_dt) {
  library(data.table)
} else {
  # Provide a simple copy function as fallback
  data.table <- list(
    copy = function(x) { 
      if (is.data.frame(x)) return(as.data.frame(x)) 
      else return(x)
    },
    data.table = function(..., keep.rownames = FALSE) {
      df <- data.frame(...)
      if (keep.rownames && !is.null(rownames(list(...)[[1]]))) {
        df <- cbind(rn = rownames(list(...)[[1]]), df)
      }
      return(df)
    }
  )
}

cat("=" , rep("=", 70), "\n", sep="")
cat("EXECUTING VIGNETTE WITH SYNTHETIC DATA\n")
cat("=", rep("=", 70), "\n\n", sep="")

# =============================================================================
# STEP 1: Create Synthetic Dataset Matching Expected Structure
# =============================================================================

cat("STEP 1: Creating synthetic dataset...\n")

set.seed(12345)
N <- 859  # Match the original dataset size mentioned in filename

# Create regional indicators (mutually exclusive)
# Approximate realistic clinical trial regional distribution
region_probs <- c(EU = 0.35, Asia = 0.25, US = 0.25, RoW = 0.15)
regions <- sample(names(region_probs), N, replace = TRUE, prob = region_probs)

# Create the synthetic case study dataset
dfcase <- data.frame(
  id = 1:N,
  regneu = ifelse(regions == "EU", "EU", "Other"),
  regasia = ifelse(regions == "Asia", "Asia", "Other"),
  regnus = ifelse(regions == "US", "US", "Other"),
  combo = rbinom(N, 1, 0.5),  # Treatment assignment (will be re-randomized)
  stratar = paste0("S", sample(1:4, N, replace = TRUE))  # Randomization strata
)

# Generate survival times using Weibull distribution
# Parameters chosen to give reasonable event rates
shape <- 1.2
scale_control <- 24  # months
scale_treat <- 28    # months (treatment benefit)

# Generate times based on treatment
dfcase$os_time <- ifelse(
  dfcase$combo == 1,
  rweibull(N, shape = shape, scale = scale_treat),
  rweibull(N, shape = shape, scale = scale_control)
)

# Add some regional variation
dfcase$os_time <- dfcase$os_time * ifelse(regions == "Asia", 1.1, 1.0)
dfcase$os_time <- dfcase$os_time * ifelse(regions == "US", 0.95, 1.0)

# Generate censoring times (administrative censoring around 36 months)
censor_time <- runif(N, 30, 42)

# Observed time and event indicator
dfcase$os_event <- as.integer(dfcase$os_time <= censor_time)
dfcase$os_time <- pmin(dfcase$os_time, censor_time)

cat("  Dataset created with", N, "subjects\n")
cat("  Regional distribution:\n")
print(table(regions))
cat("  Event rate:", round(mean(dfcase$os_event), 3), "\n")
cat("  Median follow-up:", round(median(dfcase$os_time), 1), "months\n\n")

# =============================================================================
# STEP 2: Define Functions from Vignette
# =============================================================================

cat("STEP 2: Defining vignette functions...\n")

# Function 1: get_dgm_stratified
get_dgm_stratified <- function(df, log.hrs = log(c(0.75, 0.75, 0.75, 0.75)), strata_tte = NULL) {
  
  loghr.0 <- log.hrs[1]
  loghr.z1 <- log.hrs[2]
  loghr.z2 <- log.hrs[3]
  loghr.z3 <- log.hrs[4]
  
  # Create interactions with subgroup indicators z1,z2,z3
  dfa2 <- within(df, {
    z1.treat <- z1 * treat
    z2.treat <- z2 * treat
    z3.treat <- z3 * treat
  })
  
  # strata 
  if (!is.null(strata_tte)) {
    aa <- paste("strata(", eval(strata_tte), ")")
    bb <- c("Surv(tte,event) ~ treat + z1 + z1.treat + z2 + z2.treat + z3 + z3.treat +")
    weib.formula <- as.formula(paste(bb, aa))
  }
  
  if (is.null(strata_tte)) {
    weib.formula <- as.formula("Surv(tte,event) ~ treat + z1 + z1.treat + z2 + z2.treat + z3 + z3.treat")
  }
  
  fit.weibk <- survreg(weib.formula, dist = 'weibull', data = dfa2)
  
  # Censoring model independent of covariates
  fitC.weib <- survreg(Surv(tte, 1 - event) ~ 1, dist = 'weibull', data = dfa2)
  tauC <- c(fitC.weib$scale)
  muC <- c(coef(fitC.weib)[1])
  
  mu <- c(coef(fit.weibk)[1])
  gamma <- c(coef(fit.weibk)[c(-1)])
  tau <- c(fit.weibk$scale)
  
  if (!is.null(strata_tte)) {
    strataO <- dfa2[, c(strata_tte)]
    aa <- names(tau)
    tau_ids <- unlist(lapply(strataO, function(x) { grep(x, aa) }))
    tau.strataO <- tau[tau_ids]
    if (length(tau.strataO) != nrow(dfa2)) stop("strata_tte not uniquely identified via matching")
    tau.approx <- median(tau.strataO)
    dfa2$tau.strataO <- tau.strataO
  }
  
  if (is.null(strata_tte)) {
    strataO <- "All"
    tau.strataO <- tau
    tau.approx <- tau
    dfa2$tau.strataO <- tau
  }
  
  # Re-define to satisfy log(hrs) pattern 
  b0 <- c(-gamma) / tau.approx
  b0[1] <- loghr.0
  b0[3] <- (loghr.z1 - b0[1])
  b0[5] <- (loghr.z2 - b0[1])
  b0[7] <- (loghr.z3 - b0[1])
  
  gamma.true <- -b0 * tau.approx
  beta.true <- (-1) * gamma.true / tau.approx
  
  return(list(
    df_super = dfa2,
    gamma.true = gamma.true,
    beta.true = beta.true,
    mu = mu,
    tau = tau,
    muC = muC,
    tauC = tauC,
    strata_tte = strata_tte,
    tau.approx = tau.approx
  ))
}


# Function 2: draw_sim_stratified
draw_sim_stratified <- function(dgm, ss = 1, Ndraw = nrow(dgm$df_super), 
                                 strata_rand = c("stratar"), details = TRUE) {
  
  df_super <- dgm$df_super
  
  var_names <- c(strata_rand, c("z1", "z2", "z3", "z1.treat", "z2.treat", "z3.treat"))
  if (all(var_names %in% names(df_super)) != TRUE) {
    stop("strata_rand and region variables not in dgm$df_super")
  }
  
  gamma.true <- dgm$gamma.true
  mu <- dgm$mu
  tau <- dgm$tau
  muC <- dgm$muC
  tauC <- dgm$tauC
  strata_tte <- dgm$strata_tte
  
  set.seed(8316951 + ss * 1000)
  
  if (Ndraw != nrow(df_super)) {
    # Use base R copy if data.table not available
    if (has_dt) {
      dfNew <- data.table::copy(df_super)
    } else {
      dfNew <- as.data.frame(df_super)
    }
    id_sample <- sample(c(1:nrow(dfNew)), size = Ndraw, replace = TRUE)
    df_super <- dfNew[id_sample, ]
  }
  
  tau.strataO <- df_super$tau.strataO
  N <- nrow(df_super)
  zmat <- as.matrix(df_super[, c("treat", "z1", "z1.treat", "z2", "z2.treat", "z3", "z3.treat")])
  dfsim <- df_super
  strataR <- df_super[, c(strata_rand)]
  
  if (!is.null(strata_tte)) strataO <- df_super[, c(strata_tte)]
  if (is.null(strata_tte)) strataO <- "All"
  
  # Set treatment to 1
  zmat.1 <- zmat
  zmat.1[, "treat"] <- 1.0
  zmat.1[, "z1.treat"] <- zmat.1[, "z1"]
  zmat.1[, "z2.treat"] <- zmat.1[, "z2"]
  zmat.1[, "z3.treat"] <- zmat.1[, "z3"]
  
  # Set treatment to 0
  zmat.0 <- zmat
  zmat.0[, "treat"] <- 0.0
  zmat.0[, "z1.treat"] <- 0.0
  zmat.0[, "z2.treat"] <- 0.0
  zmat.0[, "z3.treat"] <- 0.0
  
  epsilon <- log(rexp(N))
  eta1 <- mu + c(zmat.1 %*% gamma.true)
  phi1 <- (-1) * c(zmat.1 %*% gamma.true) / tau.strataO
  log.Y1 <- eta1 + tau.strataO * epsilon
  
  eta0 <- mu + c(zmat.0 %*% gamma.true)
  log.Y0 <- eta0 + tau.strataO * epsilon
  phi0 <- (-1) * c(zmat.0 %*% gamma.true) / tau.strataO
  
  loghr.po <- phi1 - phi0
  blocks <- strataR
  Zr <- block_ra(blocks = blocks)
  log.Yr <- Zr * log.Y1 + (1 - Zr) * log.Y0
  Yr <- exp(log.Yr)
  
  epsilonC <- log(rexp(N))
  log.YC <- muC + tauC * epsilonC
  log.YrC <- Zr * log.YC + (1 - Zr) * log.YC
  YrC <- exp(log.YrC)
  
  dfsim$event.sim <- ifelse(Yr <= YrC, 1, 0)
  dfsim$y.sim <- pmin(Yr, YrC)
  dfsim$treat.sim <- Zr
  dfsim$strata.simR <- strataR
  dfsim$strata.simO <- strataO
  dfsim$loghr.po <- loghr.po
  dfsim$log.Y1 <- log.Y1
  dfsim$log.Y0 <- log.Y0
  
  if (details & ss <= 10) cat("  % censored =", round(mean(1 - dfsim$event.sim), 3), "\n")
  
  if (details) {
    cat("  Stratification parm (taus) df_super:", round(tau, 4), "\n")
    
    aa <- paste("strata(", eval("strata.simO"), ")")
    bb <- c("Surv(y.sim,event.sim) ~ treat.sim + z1 + z1.treat + z2 + z2.treat + z3 + z3.treat +")
    weib.formula <- as.formula(paste(bb, aa))
    fitit <- survreg(weib.formula, dist = 'weibull', data = dfsim)
    fittau <- c(fitit$scale)
    cat("  Stratification parm (taus) simulated:", round(fittau, 4), "\n")
    
    dcheck <- loghr.po - (log.Y0 - log.Y1) / tau.strataO
    cat("  Max |loghr.po - (log.Y0-log.Y1)/tau| =", max(abs(round(dcheck, 12))), "\n")
    
    bhat.weib <- -(1) * coef(fitit)[c(-1)] / fittau
    fit.cox <- coxph(weib.formula, data = dfsim)
    fits <- cbind(bhat.weib, coef(fit.cox))
    rownames(fits) <- c("treat", "z1", "z1.treat", "z2", "z2.treat", "z3", "z3.treat")
    colnames(fits) <- c("Weibull", "Cox")
    if (has_dt) {
      fits <- data.table::data.table(fits, keep.rownames = TRUE)
    } else {
      fits <- data.frame(rn = rownames(fits), fits)
    }
    cat("\n  Parameter estimates:\n")
    print(fits)
    
    ahr_empirical <- with(dfsim, exp(mean(loghr.po)))
    cat("\n  Overall AHR =", round(ahr_empirical, 4), "\n")
    
    ahr_z1.1 <- with(subset(dfsim, z1 == 1), exp(mean(loghr.po)))
    ahr_z1.0 <- with(subset(dfsim, z1 == 0), exp(mean(loghr.po)))
    cat("  AHR Z1=1 (EU), Z1=0:", round(c(ahr_z1.1, ahr_z1.0), 4), "\n")
    
    ahr_z2.1 <- with(subset(dfsim, z2 == 1), exp(mean(loghr.po)))
    ahr_z2.0 <- with(subset(dfsim, z2 == 0), exp(mean(loghr.po)))
    cat("  AHR Z2=1 (Asia), Z2=0:", round(c(ahr_z2.1, ahr_z2.0), 4), "\n")
    
    ahr_z3.1 <- with(subset(dfsim, z3 == 1), exp(mean(loghr.po)))
    ahr_z3.0 <- with(subset(dfsim, z3 == 0), exp(mean(loghr.po)))
    cat("  AHR Z3=1 (US), Z3=0:", round(c(ahr_z3.1, ahr_z3.0), 4), "\n")
    
    ahr_z4.1 <- with(subset(dfsim, z4 == 1), exp(mean(loghr.po)))
    ahr_z4.0 <- with(subset(dfsim, z4 == 0), exp(mean(loghr.po)))
    cat("  AHR Z4=1 (RoW), Z4=0 (non-RoW):", round(c(ahr_z4.1, ahr_z4.0), 4), "\n")
    
    aa <- paste("strata(", eval("strata.simR"), ")")
    bb <- c("Surv(y.sim,event.sim) ~ treat.sim")
    coxmod1 <- as.formula(bb)
    bb <- c("Surv(y.sim,event.sim) ~ treat.sim +")
    coxmod2 <- as.formula(paste(bb, aa))
    fit1 <- summary(coxph(coxmod1, data = dfsim))$conf.int
    fit2 <- summary(coxph(coxmod2, data = dfsim))$conf.int
    cat("  Cox ITT: Un-adjusted =", round(fit1[1], 4), ", Stratified =", round(fit2[1], 4), "\n")
  }
  
  return(dfsim)
}

cat("  Functions defined successfully\n\n")

# =============================================================================
# STEP 3: Prepare Data (mimic loadData chunk)
# =============================================================================

cat("STEP 3: Preparing data with regional indicators...\n")

dfcase <- within(dfcase, {
  z1 <- ifelse(regneu == "EU", 1, 0)
  z2 <- ifelse(regasia == "Asia", 1, 0)
  z3 <- ifelse(regnus == "US", 1, 0)
  z4 <- 1 - (z1 + z2 + z3)
  treat <- combo
  tte <- os_time
  event <- os_event
})

cat("  Regional subgroup sizes:\n")
cat("    EU (z1=1):", sum(dfcase$z1), "\n")
cat("    Asia (z2=1):", sum(dfcase$z2), "\n")
cat("    US (z3=1):", sum(dfcase$z3), "\n")
cat("    RoW (z4=1):", sum(dfcase$z4), "\n\n")

# =============================================================================
# STEP 4: Test Weibull vs Cox Comparison
# =============================================================================

cat("STEP 4: Comparing Weibull vs Cox parameter estimates...\n")

fit.weib_ex <- survreg(Surv(tte, event) ~ treat + z1 + z2 + z3, dist = 'weibull', data = dfcase)
tauhat <- fit.weib_ex$scale
bhat.weib <- -(1) * coef(fit.weib_ex)[-c(1)] / tauhat

fit.cox_ex <- coxph(Surv(tte, event) ~ treat + z1 + z2 + z3, data = dfcase)
res <- data.frame(Weibull = round(bhat.weib, 4), Cox = round(coef(fit.cox_ex), 4))
cat("\n")
print(res)
cat("\n")

# =============================================================================
# STEP 5: Test DGM Functions with Uniform Effects
# =============================================================================

cat("STEP 5: Testing get_dgm_stratified with uniform HR=0.75...\n")

dgm <- get_dgm_stratified(df = dfcase, log.hrs = log(c(0.75, 0.75, 0.75, 0.75)))

cat("  Interaction parameters (should be ~0 for uniform effects):\n")
cat("    z1.treat:", round(dgm$gamma.true["z1.treat"], 6), "\n")
cat("    z2.treat:", round(dgm$gamma.true["z2.treat"], 6), "\n")
cat("    z3.treat:", round(dgm$gamma.true["z3.treat"], 6), "\n")

cat("  Overall treatment effect exp(beta[1]):", round(exp(dgm$beta.true[1]), 4), "(target: 0.75)\n\n")

# =============================================================================
# STEP 6: Draw Simulated Data
# =============================================================================

cat("STEP 6: Drawing simulated dataset (N=10000)...\n\n")

df_example <- draw_sim_stratified(dgm = dgm, ss = 123, Ndraw = 10000, strata_rand = "stratar", details = TRUE)

# =============================================================================
# STEP 7: Test Differential Effects
# =============================================================================

cat("\n", "=", rep("=", 70), "\n", sep="")
cat("STEP 7: Testing differential effects scenarios...\n")
cat("=", rep("=", 70), "\n\n", sep="")

# Scenario: z1 (EU) has weaker effect (HR=0.85 vs 0.75)
cat("Scenario A: EU (z1) has HR=0.85, others HR=0.75\n")
dgm_diff <- get_dgm_stratified(df = dfcase, log.hrs = log(c(0.75, 0.85, 0.75, 0.75)))
df_diff <- draw_sim_stratified(dgm = dgm_diff, ss = 456, Ndraw = 5000, strata_rand = "stratar", details = FALSE)

cat("  Empirical AHRs:\n")
cat("    Overall:", round(with(df_diff, exp(mean(loghr.po))), 4), "\n")
cat("    EU (z1=1):", round(with(subset(df_diff, z1 == 1), exp(mean(loghr.po))), 4), "(target: 0.85)\n")
cat("    Non-EU (z1=0):", round(with(subset(df_diff, z1 == 0), exp(mean(loghr.po))), 4), "\n\n")

# =============================================================================
# STEP 8: Run Mini-Simulation (reduced from 2000 to 50 for speed)
# =============================================================================

cat("=", rep("=", 70), "\n", sep="")
cat("STEP 8: Running mini-simulation (50 iterations instead of 2000)...\n")
cat("=", rep("=", 70), "\n\n", sep="")

sims <- 50
hr_uniform <- 0.7

hr_itt <- rep(NA, sims)
max_HRs <- rep(NA, sims)
threshold1 <- 0.80
threshold2 <- 0.90
threshold3 <- 1.0
count_thresholds1 <- rep(NA, sims)
count_thresholds2 <- rep(NA, sims)
count_thresholds3 <- rep(NA, sims)
pcensors <- rep(NA, sims)

est_sgs <- matrix(NA, nrow = sims, ncol = 8)
colnames(est_sgs) <- c("z1=1", "z1=0", "z2=1", "z2=0", "z3=1", "z3=0", "RoW", "non-RoW")

coxmod1 <- as.formula("Surv(y.sim, event.sim) ~ treat.sim")

dgm <- get_dgm_stratified(df = dfcase, log.hrs = log(c(hr_uniform, hr_uniform, hr_uniform, hr_uniform)))

cat("Target uniform HR:", hr_uniform, "\n")
cat("Running", sims, "simulations...\n\n")

pb_interval <- ceiling(sims / 10)

for (ss in 1:sims) {
  
  if (ss %% pb_interval == 0) cat("  Completed", ss, "of", sims, "\n")
  
  dfsim <- draw_sim_stratified(dgm = dgm, ss = ss, strata_rand = "stratar", details = FALSE)
  
  pcensors[ss] <- mean(1 - dfsim$event.sim)
  
  # Overall ITT
  fit <- summary(coxph(coxmod1, data = dfsim))$conf.int
  hr_itt[ss] <- c(fit[1])
  
  # Subgroup analyses
  fit <- summary(coxph(coxmod1, data = subset(dfsim, z1 == 1)))$conf.int
  est_sgs[ss, "z1=1"] <- c(fit[1])
  fit <- summary(coxph(coxmod1, data = subset(dfsim, z1 == 0)))$conf.int
  est_sgs[ss, "z1=0"] <- c(fit[1])
  
  fit <- summary(coxph(coxmod1, data = subset(dfsim, z2 == 1)))$conf.int
  est_sgs[ss, "z2=1"] <- c(fit[1])
  fit <- summary(coxph(coxmod1, data = subset(dfsim, z2 == 0)))$conf.int
  est_sgs[ss, "z2=0"] <- c(fit[1])
  
  fit <- summary(coxph(coxmod1, data = subset(dfsim, z3 == 1)))$conf.int
  est_sgs[ss, "z3=1"] <- c(fit[1])
  fit <- summary(coxph(coxmod1, data = subset(dfsim, z3 == 0)))$conf.int
  est_sgs[ss, "z3=0"] <- c(fit[1])
  
  fit <- summary(coxph(coxmod1, data = subset(dfsim, z4 == 1)))$conf.int
  est_sgs[ss, "RoW"] <- c(fit[1])
  
  fit <- summary(coxph(coxmod1, data = subset(dfsim, z4 == 0)))$conf.int
  est_sgs[ss, "non-RoW"] <- c(fit[1])
  
  max_HRs[ss] <- max(c(est_sgs[ss, ]))
  
  count_thresholds1[ss] <- c(max_HRs[ss] >= threshold1)
  count_thresholds2[ss] <- c(max_HRs[ss] >= threshold2)
  count_thresholds3[ss] <- c(max_HRs[ss] >= threshold3)
}

# =============================================================================
# STEP 9: Simulation Results
# =============================================================================

cat("\n", "=", rep("=", 70), "\n", sep="")
cat("SIMULATION RESULTS\n")
cat("=", rep("=", 70), "\n\n", sep="")

cat("Number of simulations:", sims, "\n")
cat("Average % censoring:", round(mean(pcensors) * 100, 1), "%\n\n")

cat("--- Threshold Analysis ---\n")
cat("Threshold 1 =", threshold1, "\n")
cat("  % of ITT HR >= threshold:", round(mean(ifelse(hr_itt >= threshold1, 1, 0)) * 100, 1), "%\n")
cat("  % of Max(subgroup HRs) >= threshold:", round(mean(ifelse(max_HRs >= threshold1, 1, 0)) * 100, 1), "%\n\n")

cat("Threshold 2 =", threshold2, "\n")
cat("  % of ITT HR >= threshold:", round(mean(ifelse(hr_itt >= threshold2, 1, 0)) * 100, 1), "%\n")
cat("  % of Max(subgroup HRs) >= threshold:", round(mean(ifelse(max_HRs >= threshold2, 1, 0)) * 100, 1), "%\n\n")

cat("Threshold 3 =", threshold3, "\n")
cat("  % of ITT HR >= threshold:", round(mean(ifelse(hr_itt >= threshold3, 1, 0)) * 100, 1), "%\n")
cat("  % of Max(subgroup HRs) >= threshold:", round(mean(ifelse(max_HRs >= threshold3, 1, 0)) * 100, 1), "%\n\n")

cat("--- Quantiles ---\n")
cat("ITT HR quantiles:\n")
print(round(quantile(hr_itt), 4))
cat("\nMax(subgroup HRs) quantiles:\n")
print(round(quantile(max_HRs), 4))

cat("\n--- Subgroup-Specific Results ---\n")
cat("Mean HR by subgroup:\n")
print(round(colMeans(est_sgs), 4))

cat("\n", "=", rep("=", 70), "\n", sep="")
cat("EXECUTION COMPLETED SUCCESSFULLY\n")
cat("=", rep("=", 70), "\n", sep="")
