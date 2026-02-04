# ==========================================
# CP ~ (SS | TSS | TRIMP) via 2-EWMA gain/fatigue model
# Independent fits with identical bounds; compare SSE/RMSE/AIC
# Requires: minpack.lm, readxl
# Input Excel file columns: t,	CP,	SS, TSS, TRIMP
# ==========================================

# Packages
library(minpack.lm)
library(readxl)

# ---------- User settings ----------
in_xlsx <- "Example data for fitting CP with SS TSS TRIMP.xlsx"
sheet   <- "Sheet1"

# Baseline (untrained) CP for the performance equation p(t) = p0 + k1*G - k2*H
p0_cp <- 150   # W (adjust if you have a better prior)

# Initial state for the EWMA states (load at t=0 before records start)
SS0     <- 90
TSS0    <- 70
TRIMP0  <- 110

# Bounds & starts (common to all 3 fits); order: (k1, k2, tau1, tau2)
# tau1 = "fitness" time-constant, tau2 = "fatigue" time-constant (days or sample steps)
lb  <- c(0.2, 0.1, 15,  3)
ub  <- c(5.0, 4.0, 52, 10)
st  <- c(1.0, 0.5, 45,  5)

# ---------- Load data ----------
dat <- read_excel(in_xlsx, sheet = sheet)

num <- function(x) suppressWarnings(as.numeric(x))
t_col     <- if ("t" %in% names(dat)) num(dat$t) else seq_len(nrow(dat)) - 1
CP_obs    <- num(dat$CP)
SS        <- num(dat$SS)
TSS       <- num(dat$TSS)
TRIMP     <- num(dat$TRIMP)

# Clean NA loads to zero (sensible default for cumulative indices)
SS[!is.finite(SS)]         <- 0
TSS[!is.finite(TSS)]       <- 0
TRIMP[!is.finite(TRIMP)]   <- 0

# Observation indices for CP
idx_cp <- which(is.finite(CP_obs)) - 1
y_cp   <- CP_obs[is.finite(CP_obs)]
if (length(y_cp) < 3) stop("Need at least 3 CP observations to estimate parameters.")

tseq <- seq_len(nrow(dat)) - 1

# ---------- EWMA & model ----------
compute_ewma <- function(w, tau, init = 0) {
  n <- length(w)
  out <- numeric(n)
  a <- 1 - exp(-1 / tau)
  if (!is.finite(a) || a <= 0) a <- 1e-8
  prev <- init
  for (i in 1:n) {
    prev <- prev * (1 - a) + w[i] * a
    out[i] <- prev
  }
  out
}

# Performance trajectory for a single metric series
# par5 = (p0, k1, k2, tau1, tau2)
perf_from_params <- function(par5, load_series, init_load = 0) {
  p0   <- par5[1]; k1 <- par5[2]; k2 <- par5[3]; tau1 <- par5[4]; tau2 <- par5[5]
  G <- compute_ewma(load_series, tau1, init = init_load)  # fitness
  H <- compute_ewma(load_series, tau2, init = init_load)  # fatigue
  p <- p0 + k1 * G - k2 * H
  list(p = p, G = G, H = H)
}

# ---------- Boxed residual (unconstrained z -> [L,U]) ----------
map_box <- function(z, L, U) L + (U - L) / (1 + exp(-z))

# theta order: (k1, k2, tau1, tau2)
resid_boxed_cp <- function(theta_z, load_series, idx_obs, y_obs, init_load, p0,
                           L, U, penalize_tau = TRUE) {
  theta <- map_box(theta_z, L, U)
  par5  <- c(p0, theta)
  pred  <- perf_from_params(par5, load_series, init_load = init_load)$p
  # optional gentle regularizer encouraging tau2 <= tau1 (fatigue faster than fitness)
  pen <- if (penalize_tau) 1e-3 * max(0, theta[4] - theta[3]) else 0
  residuals <- pred[idx_obs + 1] - y_obs
  residuals + pen
}

# Robust z-start builder
to_zstarts <- function(st, L, U, name = "block") {
  if (length(st) != 4L || length(L) != 4L || length(U) != 4L) {
    stop(sprintf("'%s' starts/bounds must be length 4. Got: st=%d, L=%d, U=%d",
                 name, length(st), length(L), length(U)))
  }
  st_clipped <- pmin(pmax(st, L + 1e-9), U - 1e-9)
  p <- (st_clipped - L) / (U - L)
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  as.numeric(qlogis(p))
}

st_z <- to_zstarts(st, lb, ub, "CP")

# ---------- Fit helper ----------
fit_one_metric <- function(metric_name, load_series, init_load) {
  fit <- minpack.lm::nls.lm(
    par   = st_z,
    fn    = resid_boxed_cp,
    # boxing inside fn(): do NOT pass lower/upper here
    load_series = load_series,
    idx_obs = idx_cp, y_obs = y_cp,
    init_load = init_load, p0 = p0_cp,
    L = lb, U = ub, penalize_tau = TRUE,
    control = minpack.lm::nls.lm.control(maxiter = 500, ftol = 1e-10, ptol = 1e-10)
  )
  theta_hat <- map_box(fit$par, lb, ub)         # (k1,k2,tau1,tau2) in-bounds
  par5      <- setNames(c(p0_cp, theta_hat), c("p0","k1","k2","tau1","tau2"))
  traj      <- perf_from_params(par5, load_series, init_load = init_load)
  pred_obs  <- traj$p[idx_cp + 1]
  sse       <- sum((pred_obs - y_cp)^2)
  n         <- length(y_cp)
  k_params  <- 4  # (k1,k2,tau1,tau2); p0 is fixed here
  rmse      <- sqrt(sse / n)
  # AIC for least-squares with Gaussian errors: n*log(SSE/n) + 2*k
  aic       <- n * log(sse / n) + 2 * k_params
  
  list(
    metric = metric_name,
    par    = par5,
    traj   = traj,
    SSE    = sse,
    RMSE   = rmse,
    AIC    = aic,
    N      = n
  )
}

# ---------- Run all three independent fits ----------
res_SS    <- fit_one_metric("SS",    SS,    SS0)
res_TSS   <- fit_one_metric("TSS",   TSS,   TSS0)
res_TRIMP <- fit_one_metric("TRIMP", TRIMP, TRIMP0)

# ---------- Summaries ----------
sum_df <- rbind(
  data.frame(metric = res_SS$metric,    t(res_SS$par),    SSE = res_SS$SSE,    RMSE = res_SS$RMSE,    AIC = res_SS$AIC,    N = res_SS$N,    row.names = NULL, check.names = FALSE),
  data.frame(metric = res_TSS$metric,   t(res_TSS$par),   SSE = res_TSS$SSE,   RMSE = res_TSS$RMSE,   AIC = res_TSS$AIC,   N = res_TSS$N,   row.names = NULL, check.names = FALSE),
  data.frame(metric = res_TRIMP$metric, t(res_TRIMP$par), SSE = res_TRIMP$SSE, RMSE = res_TRIMP$RMSE, AIC = res_TRIMP$AIC, N = res_TRIMP$N, row.names = NULL, check.names = FALSE)
)

write.csv(sum_df, "cp_vs_load_fits.csv", row.names = FALSE)

# ---------- Time-series export ----------
ts_out <- data.frame(
  t      = t_col,
  CP_obs = CP_obs,
  CP_SS_fit    = res_SS$traj$p,
  CP_TSS_fit   = res_TSS$traj$p,
  CP_TRIMP_fit = res_TRIMP$traj$p,
  SS    = SS,    TSS = TSS,    TRIMP = TRIMP,
  SS_G  = res_SS$traj$G,  SS_H  = res_SS$traj$H,
  TSS_G = res_TSS$traj$G, TSS_H = res_TSS$traj$H,
  TRIMP_G = res_TRIMP$traj$G, TRIMP_H = res_TRIMP$traj$H
)
write.csv(ts_out, "cp_vs_load_timeseries.csv", row.names = FALSE)

## ---------- Plot overlay (interactive, not saved) ----------
op <- par(mar = c(4.5, 5.5, 3.5, 2.0))
on.exit(par(op), add = TRUE)

yr <- range(c(CP_obs, res_SS$traj$p, res_TSS$traj$p, res_TRIMP$traj$p), na.rm = TRUE)
plot(tseq, res_SS$traj$p, type = "l", lwd = 2, col = "green",
     xlab = "Time (days)", ylab = "CP (W)",
     ylim = yr, main = "")
lines(tseq, res_TSS$traj$p,   lwd = 2, col = "#1f77b4")
lines(tseq, res_TRIMP$traj$p, lwd = 2, col = "red")
points(tseq[idx_cp + 1], y_cp, pch = 19, cex = 1.0)
legend("bottomright", bty = "n",
       col = c("green", "#1f77b4", "red", "black"),
       lwd = c(2,2,2,NA), pch = c(NA,NA,NA,19),
       legend = c(
         sprintf("SS  (RMSE %.1f)",  res_SS$RMSE),
         sprintf("TSS (RMSE %.1f)",  res_TSS$RMSE),
         sprintf("TRIMP (RMSE %.1f)",res_TRIMP$RMSE),
         "Observed CP"))

cat("Saved:\n - cp_vs_load_fits.csv\n - cp_vs_load_timeseries.csv\n - cp_vs_load_plot.png\n")
