# ==========================================
# CP, W', Pmax ~ via 2-EWMA gain/fatigue model
# Independent fits; compare SSE/RMSE/AIC
# Requires: minpack.lm, readxl
# Input Excel file columns: t, CP, Wprime, Pmax,	SSCP,	SSW, SSPmax
# ==========================================

# Packages
library(minpack.lm)
library(readxl)

SSCP0   <- 90   # baseline for CP load
SSW0    <- 12    # baseline for W′ load
SSPmax0 <- 2     # baseline for Pmax load

p0_cp   <- 150     # untrained CP (W)
p0_w    <- 7500    # untrained W′ (J)
p0_pmax <- 600     # untrained Pmax (W)

# --- Load data ---
data <- read_excel("Example data for fitting 3D.xlsx", sheet = "Sheet1")

# --- Bounds & starts per outcome, order: c(k1, k2, tau1, tau2) ---
lb_cp    <- c(0.2, 0.1,  15,  3)
ub_cp    <- c(5,   4,   52, 10)
st_cp    <- c(1, 0.5, 45,  5)

lb_w     <- c(2000,   2000,    3,  3)
ub_w     <- c(6000,6000, 40, 14)
st_w     <- c(4000,4000,  30,  5)

lb_pmax  <- c(0,   0,    7,  1)
ub_pmax  <- c(500, 400,  10, 10)
st_pmax  <- c(30,  6,   20,  5)

# Coerce numeric & guard
num <- function(x) suppressWarnings(as.numeric(x))
CP     <- num(data$CP)
Wprime <- num(data$Wprime)
Pmax   <- num(data$Pmax)

SSCP    <- num(data$SSCP);    SSCP[!is.finite(SSCP)]     <- 0
SSW     <- num(data$SSW);     SSW[!is.finite(SSW)]       <- 0
SSPmax  <- num(data$SSPmax);  SSPmax[!is.finite(SSPmax)] <- 0

tseq <- seq_len(nrow(data)) - 1

# Observation indices (zero-based for p[t] addressing)
idx_cp    <- which(is.finite(CP))     - 1
idx_w     <- which(is.finite(Wprime)) - 1
idx_pmax  <- which(is.finite(Pmax))   - 1

obs_cp    <- CP[is.finite(CP)]
obs_w     <- Wprime[is.finite(Wprime)]
obs_pmax  <- Pmax[is.finite(Pmax)]

# Weights to balance scales (avoid one metric dominating)
w_cp    <- if (length(obs_cp)   > 1 && sd(obs_cp)   > 0) 1/sd(obs_cp)   else 1
w_w     <- if (length(obs_w)    > 1 && sd(obs_w)    > 0) 1/sd(obs_w)    else 1
w_pmax  <- if (length(obs_pmax) > 1 && sd(obs_pmax) > 0) 1/sd(obs_pmax) else 1

# --- EWMA with init ---
compute_ewma <- function(w, tau, init = 0) {
  n <- length(w); ew <- numeric(n)
  a <- 1 - exp(-1 / tau); if (!is.finite(a) || a <= 0) a <- 1e-8
  prev <- init
  for (i in 1:n) { prev <- prev * (1 - a) + w[i] * a; ew[i] <- prev }
  ew
}

# --- Helper: one metric trajectory ---
perf_from_params <- function(par5, ss, init_ss = 0) {
  p0 <- par5[1]; k1 <- par5[2]; k2 <- par5[3]; tau1 <- par5[4]; tau2 <- par5[5]
  g <- compute_ewma(ss, tau1, init = init_ss)
  h <- compute_ewma(ss, tau2, init = init_ss)
  p <- p0 + k1 * g - k2 * h
  list(p = p, g = g, h = h)
}

# ---------- Bounds-enforcing residual (boxed reparameterization) ----------
# map z (R) -> [L,U] with logistic; guarantees staying inside box
map_box <- function(z, L, U) L + (U - L) / (1 + exp(-z))  # note: exp(-z), not -x

# theta_z are unconstrained; we map them to [L,U] each call
# theta order: c(k1, k2, tau1, tau2)
resid_one_boxed <- function(theta_z, ss, idx_obs, y_obs, w_obs, init_ss, p0,
                            L, U, penalize_tau = TRUE) {
  theta <- map_box(theta_z, L, U)                        # k1,k2,tau1,tau2 in-bounds
  par5  <- c(p0, theta)
  pred  <- perf_from_params(par5, ss, init_ss = init_ss)$p
  pen   <- if (penalize_tau) 1e-3 * max(0, theta[4] - theta[3]) else 0  # tau2 > tau1
  w_obs * (pred[idx_obs + 1] - y_obs) + pen
}

# helper to get z-start near desired mid-point in the box
mid_to_z <- function(m, L, U) {
  # (m-L)/(U-L) in (0,1); clamp just in case
  p <- (m - L) / (U - L)
  p <- min(max(p, 1e-6), 1 - 1e-6)
  qlogis(p)
}

# --- Fit (independent problems) ---

# --- Robust z-start builder (replaces mid_to_z & old st_*_z lines) ---
to_zstarts <- function(st, L, U, name = "block") {
  if (length(st) != 4L || length(L) != 4L || length(U) != 4L) {
    stop(sprintf("'%s' starts/bounds must be length 4. Got: st=%d, L=%d, U=%d",
                 name, length(st), length(L), length(U)))
  }
  # ensure starts are within the box (clip gently if needed)
  st_clipped <- pmin(pmax(st, L + 1e-9), U - 1e-9)
  p <- (st_clipped - L) / (U - L)       # should be in (0,1)
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)    # guard edges
  z <- as.numeric(qlogis(p))
  if (any(!is.finite(z))) {
    # Fallback to box mid-point (z=0) if anything went weird
    warning(sprintf("Non-finite z starts in '%s'; falling back to zeros.", name))
    z <- rep(0, 4)
  }
  z
}

# Build z-starts (length 4 each) — replaces previous st_*_z construction
st_cp_z  <- to_zstarts(st_cp,   lb_cp,   ub_cp,   "CP")
st_w_z   <- to_zstarts(st_w,    lb_w,    ub_w,    "Wprime")
st_pm_z  <- to_zstarts(st_pmax, lb_pmax, ub_pmax, "Pmax")

# Optional: quick diagnostics
cat("z-start lengths: CP=", length(st_cp_z), 
    " W'=", length(st_w_z), 
    " Pmax=", length(st_pm_z), "\n")

## ---- CP ----
fit_cp <- minpack.lm::nls.lm(
  par   = st_cp_z,
  fn    = resid_one_boxed,
  # no lower/upper here; boxing is handled in resid_one_boxed via map_box()
  ss = SSCP, idx_obs = idx_cp, y_obs = obs_cp, w_obs = w_cp,
  init_ss = SSCP0, p0 = p0_cp, L = lb_cp, U = ub_cp, penalize_tau = TRUE,
  control = minpack.lm::nls.lm.control(maxiter = 300, ftol = 1e-10, ptol = 1e-10)
)
theta_cp <- map_box(fit_cp$par, lb_cp, ub_cp)
par_cp   <- setNames(c(p0_cp, theta_cp), c("p0","k1","k2","tau1","tau2"))

## ---- W′ ----
fit_w <- minpack.lm::nls.lm(
  par   = st_w_z,
  fn    = resid_one_boxed,
  ss = SSW, idx_obs = idx_w, y_obs = obs_w, w_obs = w_w,
  init_ss = SSW0, p0 = p0_w, L = lb_w, U = ub_w, penalize_tau = TRUE,
  control = minpack.lm::nls.lm.control(maxiter = 300, ftol = 1e-10, ptol = 1e-10)
)
theta_w <- map_box(fit_w$par, lb_w, ub_w)
par_w   <- setNames(c(p0_w, theta_w), c("p0","k1","k2","tau1","tau2"))

## ---- Pmax ----
fit_pm <- minpack.lm::nls.lm(
  par   = st_pm_z,
  fn    = resid_one_boxed,
  ss = SSPmax, idx_obs = idx_pmax, y_obs = obs_pmax, w_obs = w_pmax,
  init_ss = SSPmax0, p0 = p0_pmax, L = lb_pmax, U = ub_pmax, penalize_tau = TRUE,
  control = minpack.lm::nls.lm.control(maxiter = 300, ftol = 1e-10, ptol = 1e-10)
)
theta_pm <- map_box(fit_pm$par, lb_pmax, ub_pmax)
par_pmax <- setNames(c(p0_pmax, theta_pm), c("p0","k1","k2","tau1","tau2"))


# --- Print key outputs ---
cat("CP:",    par_cp,   "\n")
cat("W′:",    par_w,    "\n")
cat("Pmax:",  par_pmax, "\n")

# Recompute trajectories
cp_res   <- perf_from_params(par_cp,   SSCP,   init_ss = SSCP0)
w_res    <- perf_from_params(par_w,    SSW,    init_ss = SSW0)
pmax_res <- perf_from_params(par_pmax, SSPmax, init_ss = SSPmax0)

# --- SSEs (unweighted at observed points) ---
sse_cp    <- sum((cp_res$p[idx_cp + 1]     - obs_cp)^2)
sse_w     <- sum((w_res$p[idx_w + 1]       - obs_w)^2)
sse_pmax  <- sum((pmax_res$p[idx_pmax + 1] - obs_pmax)^2)
sse_total <- sse_cp + sse_w + sse_pmax

cat("=== CP params ===\n");    print(par_cp)
cat("=== W' params ===\n");    print(par_w)
cat("=== Pmax params ===\n");  print(par_pmax)
cat(sprintf("SSE: CP=%.3f, W'=%.3f, Pmax=%.3f, Total=%.3f\n",
            sse_cp, sse_w, sse_pmax, sse_total))

# --- Save results ---
out_params <- rbind(
  cbind(metric = "CP",      as.data.frame(t(par_cp),   check.names = FALSE)),
  cbind(metric = "Wprime",  as.data.frame(t(par_w),    check.names = FALSE)),
  cbind(metric = "Pmax",    as.data.frame(t(par_pmax), check.names = FALSE))
)
out_params$SSE_metric <- c(sse_cp, sse_w, sse_pmax)
write.csv(out_params, "fit_results_multi.csv", row.names = FALSE)

ts_out <- data.frame(
  day = tseq,
  SSCP = SSCP, SSW = SSW, SSPmax = SSPmax,
  CP_fit = cp_res$p, CP_g = cp_res$g, CP_h = cp_res$h,
  Wprime_fit = w_res$p, Wprime_g = w_res$g, Wprime_h = w_res$h,
  Pmax_fit = pmax_res$p, Pmax_g = pmax_res$g, Pmax_h = pmax_res$h
)
write.csv(ts_out, "model_timeseries_multi.csv", row.names = FALSE)

# ===== Separate plots: CP, W′, Pmax =====
op <- par(mfrow = c(3,1), mar = c(4,4,2,1))

# 1) CP
yl_cp <- range(c(cp_res$p, obs_cp), na.rm = TRUE)
yl_cp <- yl_cp + c(-0.05, 0.05) * diff(yl_cp)
plot(tseq, cp_res$p, type = "l", col = "blue", lwd = 2,
     ylim = yl_cp, xlab = "Day", ylab = "CP (W)", main = "CP: model vs. observations")
points(idx_cp, cp_res$p[idx_cp + 1], pch = 1, col = "blue")
points(idx_cp, obs_cp, pch = 19, col = "black")
legend("bottomright", c("Model","Observed"), col = c("blue","black"), pch = c(1,19), lwd = c(2,NA))

# 2) W′
yl_w <- range(c(w_res$p, obs_w), na.rm = TRUE)
yl_w <- yl_w + c(-0.05, 0.05) * diff(yl_w)
plot(tseq, w_res$p, type = "l", col = "darkorange", lwd = 2,
     ylim = yl_w, xlab = "Day", ylab = "W′ (J)", main = "W′: model vs. observations")
points(idx_w, w_res$p[idx_w + 1], pch = 1, col = "darkorange")
points(idx_w, obs_w, pch = 19, col = "black")
legend("bottomright", c("Model","Observed"), col = c("darkorange","black"), pch = c(1,19), lwd = c(2,NA))

# 3) Pmax
yl_pm <- range(c(pmax_res$p, obs_pmax), na.rm = TRUE)
yl_pm <- yl_pm + c(-0.05, 0.05) * diff(yl_pm)
plot(tseq, pmax_res$p, type = "l", col = "purple", lwd = 2,
     ylim = yl_pm, xlab = "Day", ylab = "Pmax (W)", main = "Pmax: model vs. observations")
points(idx_pmax, pmax_res$p[idx_pmax + 1], pch = 1, col = "purple")
points(idx_pmax, obs_pmax, pch = 19, col = "black")
legend("bottomright", c("Model","Observed"), col = c("purple","black"), pch = c(1,19), lwd = c(2,NA))

par(op)
