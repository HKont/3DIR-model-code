# ==========================================
# Power–Duration (PD) curve fitting: 3CP (Morton model; k=1) & 3CPmod (Xert model; k=2)
# - Input: "PD data.csv"
#   Example format:
#   Month,5,60,180,360,720,1200
#   Sep-24,680,399.4,298.4,262.4,239.4,230.6
#
# - Choose durations to use via durations_use.
# - Outputs:
#   * pd_fit_summary.csv   (per Month x Model row with params, SEs, CIs, RMSE)
#   * pd_fit_<Month>.png   (overlay plot with both models + 95% CI ribbons)
#
# Base R only.
# ==========================================

## ---- User controls ----
input_csv      <- "PD data.csv"
output_summary <- "pd_fit_summary.csv"
plot_prefix    <- "pd_fit_"

# Durations (seconds) to fit — must match column names in the CSV
durations_use  <- c(5, 30, 60, 120, 150, 180, 360, 480, 600, 720, 1200)

# Bounds for parameters in nls(port)
bounds <- list(
  CP    = c(  150,  300),   # W
  Wprime= c( 2000, 80000),  # J
  Pmax  = c( 500, 2500)    # W
)

# Numerical-diff step for CI of predicted curves
eps_grad <- 1e-4

## ---- Helpers: model equations ----
# 3CP (k=1): from t = W'/(P-CP) - W'/(Pmax-CP) => P(t) = CP + W' / ( t + W'/(Pmax-CP) )
morton_P_of_t <- function(CP, Wp, Pmax, t) {
  CP + Wp / ( t + Wp / pmax(1e-9, (Pmax - CP)) )
}

# 3CPmod (k=2): solve quadratic in y=P-CP; see your earlier derivation
xert_P_of_t <- function(CP, Wp, Pmax, t) {
  A <- Pmax; B <- CP
  c <- (A - B) * (t^2) / (Wp^2)
  disc <- 1 + 4 * c * (A - B)
  y <- ifelse(c > 0, (-1 + sqrt(disc)) / (2 * c), (A - B))  # t -> 0 => P -> Pmax
  pmax(B, pmin(A, B + y))
}

## ---- Starting values heuristics ----
# Use the longest duration as a CP proxy; W' via mid-duration; Pmax via short sprint
make_start_vals <- function(t, P) {
  o <- order(t); t <- t[o]; P <- P[o]
  CP0   <- max(50, min(P[length(P)] - 10, 0.98 * P[length(P)]))
  Pmax0 <- max(P[1] * 1.05, CP0 + 50)
  # Rough W': use 60–180s region if available, else mid point
  idx_mid <- which.min(abs(t - 120))
  if (length(idx_mid) == 0 || is.na(idx_mid)) idx_mid <- ceiling(length(t)/2)
  Wp0 <- max(500, (P[idx_mid] - CP0) * t[idx_mid])
  list(CP=CP0, Wprime=Wp0, Pmax=Pmax0)
}

## ---- Fit wrappers using nls(port) with bounds ----
fit_morton <- function(t, P) {
  sv <- make_start_vals(t, P)
  lower <- c(CP=bounds$CP[1], Wprime=bounds$Wprime[1], Pmax=bounds$Pmax[1])
  upper <- c(CP=bounds$CP[2], Wprime=bounds$Wprime[2], Pmax=bounds$Pmax[2])
  
  # guard: require Pmax > CP
  if (sv$Pmax <= sv$CP + 5) sv$Pmax <- sv$CP + 50
  
  nls(P ~ CP + Wprime / ( t + Wprime / (Pmax - CP) ),
      start = list(CP=sv$CP, Wprime=sv$Wprime, Pmax=sv$Pmax),
      algorithm = "port", lower = lower, upper = upper,
      control = nls.control(maxiter=500, warnOnly=TRUE))
}

fit_xert <- function(t, P) {
  sv <- make_start_vals(t, P)
  lower <- c(CP=bounds$CP[1], Wprime=bounds$Wprime[1], Pmax=bounds$Pmax[1])
  upper <- c(CP=bounds$CP[2], Wprime=bounds$Wprime[2], Pmax=bounds$Pmax[2])
  
  # We can't write xert_P_of_t directly in nls formula, so wrap as a self-start
  # Use a user-defined function returning predicted P given params & t
  xert_fun <- function(CP, Wprime, Pmax, t) xert_P_of_t(CP, Wprime, Pmax, t)
  
  nls(P ~ xert_fun(CP, Wprime, Pmax, t),
      start = list(CP=sv$CP, Wprime=sv$Wprime, Pmax=sv$Pmax),
      algorithm = "port", lower = lower, upper = upper,
      control = nls.control(maxiter=500, warnOnly=TRUE))
}

## ---- Prediction with CI via delta method (numeric gradient) ----
pred_with_ci <- function(f_P_of_t, par_vec, vcov_mat, tgrid, eps=1e-4) {
  # f_P_of_t: function(par, t) -> P
  Pfit <- f_P_of_t(par_vec, tgrid)
  
  grad_num <- function(par, t) {
    g <- numeric(length(par))
    for (j in seq_along(par)) {
      step <- rep(0, length(par)); step[j] <- eps * (abs(par[j]) + 1)
      p_up <- par + step; p_dn <- par - step
      f_up <- f_P_of_t(p_up, t); f_dn <- f_P_of_t(p_dn, t)
      g[j] <- (f_up - f_dn) / (2 * step[j])
    }
    g
  }
  
  se <- sapply(tgrid, function(tt) {
    g <- grad_num(par_vec, tt)
    sqrt( max(0, as.numeric(t(g) %*% vcov_mat %*% g)) )
  })
  
  list(P = Pfit,
       lo = Pfit - 1.96 * se,
       hi = Pfit + 1.96 * se,
       se = se)
}

## ---- Plot per Month (both models + CI ribbons) ----
plot_month <- function(month_lab, t, P,
                       morton_fit, xert_fit,
                       tgrid, morton_ci, xert_ci,
                       outfile) {
  png(outfile, width = 1200, height = 700)
  par(mar = c(4.5, 5.5, 3.5, 2.0))
  yr <- range(c(P, morton_ci$lo, morton_ci$hi, xert_ci$lo, xert_ci$hi), na.rm = TRUE)
  xr <- range(c(1, max(tgrid)))
  
  # empty plot
  plot(NA, NA, xlim = xr, ylim = yr, xlab = "Duration t (s)", ylab = "Power (W)",
       main = sprintf("PD fits: %s", month_lab))
  
  # ribbons
  polygon(c(tgrid, rev(tgrid)), c(morton_ci$lo, rev(morton_ci$hi)),
          col = "#1f77b433", border = NA)
  polygon(c(tgrid, rev(tgrid)), c(xert_ci$lo, rev(xert_ci$hi)),
          col = "#ff7f0e33", border = NA)
  
  # curves
  lines(tgrid, morton_ci$P, lwd = 2, col = "#1f77b4") # Morton (k=1)
  lines(tgrid, xert_ci$P,   lwd = 2, col = "#ff7f0e") # Xert (k=2)
  
  # observed points
  points(t, P, pch = 19, cex = 1.2)
  text(t, P, labels = paste0(t, "s"), pos = 3, cex = 0.8)
  
  legend("topright", bty = "n",
         lwd = c(2, 2, NA),
         pch = c(NA, NA, 19),
         col = c("#1f77b4", "#ff7f0e", "black"),
         legend = c("Morton (k=1)", "Xert (k=2)", "Observed"))
  grid()
  dev.off()
}

## ---- Main: read, loop over rows (months), fit, summarize & plot ----
raw <- read.csv(input_csv, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Check durations exist as columns
dur_labels <- as.character(durations_use)
missing_cols <- setdiff(dur_labels, colnames(raw))
if (length(missing_cols)) stop("These duration columns are missing in PD data.csv: ", paste(missing_cols, collapse=", "))

# Prepare output summary collector
sum_rows <- list()

for (ri in seq_len(nrow(raw))) {
  month_lab <- as.character(raw$Month[ri])
  t <- durations_use
  P <- as.numeric(raw[ri, dur_labels])
  
  # Drop NA durations
  keep <- is.finite(P)
  t <- t[keep]; P <- P[keep]
  if (length(t) < 3) {
    warning("Row ", ri, " (", month_lab, "): fewer than 3 data points; skipping.")
    next
  }
  
  # Fit Morton
  fit1 <- try(fit_morton(t, P), silent = TRUE)
  # Fit Xert
  fit2 <- try(fit_xert(t, P),   silent = TRUE)
  
  # Functions for prediction & names
  f1 <- function(par, tt) morton_P_of_t(par[1], par[2], par[3], tt)
  f2 <- function(par, tt) xert_P_of_t(par[1], par[2], par[3], tt)
  
  # Build plot grid
  tgrid <- seq(1, max(t) * 1.15, by = 1)
  
  # Aggregate results
  add_result <- function(model_name, fit_obj, f_pred) {
    if (inherits(fit_obj, "try-error")) {
      return(data.frame(
        Month = month_lab, Model = model_name,
        CP = NA, CP_SE = NA, CP_L = NA, CP_U = NA,
        Wprime = NA, Wprime_SE = NA, Wprime_L = NA, Wprime_U = NA,
        Pmax = NA, Pmax_SE = NA, Pmax_L = NA, Pmax_U = NA,
        RMSE = NA, stringsAsFactors = FALSE))
    }
    sm <- summary(fit_obj)
    co <- coef(sm)
    par_hat <- coef(fit_obj)
    names(par_hat) <- names(coef(fit_obj))
    
    # Extract SE and 95% CI
    se <- sqrt(diag(vcov(fit_obj)))
    ciL <- par_hat - 1.96 * se
    ciU <- par_hat + 1.96 * se
    
    # RMSE
    rmse <- sqrt(mean(residuals(fit_obj)^2, na.rm = TRUE))
    
    # Prediction with 95% CI (for plotting)
    vc <- vcov(fit_obj)
    ci_pred <- pred_with_ci(function(par, tt) f_pred(par, tt),
                            par_vec = par_hat,
                            vcov_mat = vc,
                            tgrid = tgrid,
                            eps = eps_grad)
    
    list_row <- data.frame(
      Month = month_lab, Model = model_name,
      CP = par_hat["CP"], CP_SE = se["CP"], CP_L = ciL["CP"], CP_U = ciU["CP"],
      Wprime = par_hat["Wprime"], Wprime_SE = se["Wprime"], Wprime_L = ciL["Wprime"], Wprime_U = ciU["Wprime"],
      Pmax = par_hat["Pmax"], Pmax_SE = se["Pmax"], Pmax_L = ciL["Pmax"], Pmax_U = ciU["Pmax"],
      RMSE = rmse, stringsAsFactors = FALSE)
    
    attr(list_row, "ci_pred") <- ci_pred
    attr(list_row, "par_hat") <- par_hat
    list_row
  }
  
  r1 <- add_result("Morton (k=1)", fit1, f1)
  r2 <- add_result("Xert (k=2)",   fit2, f2)
  
  # Save plot if at least one fit succeeded
  if (!all(is.na(c(r1$CP, r2$CP)))) {
    morton_ci <- if (!is.null(attr(r1, "ci_pred"))) attr(r1, "ci_pred") else list(P=numeric(), lo=numeric(), hi=numeric())
    xert_ci   <- if (!is.null(attr(r2, "ci_pred"))) attr(r2, "ci_pred") else list(P=numeric(), lo=numeric(), hi=numeric())
    outfile   <- paste0(plot_prefix, gsub("[^A-Za-z0-9_-]", "_", month_lab), ".png")
    plot_month(month_lab, t, P, fit1, fit2, tgrid, morton_ci, xert_ci, outfile)
  }
  
  sum_rows[[length(sum_rows)+1]] <- r1
  sum_rows[[length(sum_rows)+1]] <- r2
}

# Write summary CSV
if (length(sum_rows)) {
  summary_df <- do.call(rbind, sum_rows)
  write.csv(summary_df, output_summary, row.names = FALSE)
  message("Saved: ", output_summary)
} else {
  stop("No rows could be fitted. Check your data file and chosen durations.")
}

# =========================
# End of script
# =========================
