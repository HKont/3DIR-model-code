# ---- Parameters (specific to athlete) ----
CP     <- 240    # watts
Wprime <- 18000   # joules
Pmax   <- 950     # watts
dt     <- 1.0     # seconds per sample
# ---- TRIMP options ----
GENDER <- "female"   # "male" or "female"

# ---- Libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(tools)
  library(lubridate)
  library(data.table)
})

# =========================
# Strain Score Batch Script 
# =========================
# - Expects 1 Hz power in column G (7th column) of each CSV
# - Uses Skiba differential W′ balance (spreadsheet E-update on W′ expended)
# - Creates per-ride timeseries CSV, QA CSV, and stacked SS bar PNG
# - Aggregates QA across rides + a combined stacked composition chart
# - Adds start_date, start_time to QA outputs
# --------------------------------------------------------------

ride_dir <- "Power files"     # folder with ride CSVs
outdir   <- "Outputs"   # where to write outputs
tz_local <- "America/Edmonton"

# ---- Colors for plots ----
cols <- c("CP" = "#1f77b4",    # blue
          "W′" = "#ff7f0e",    # orange
          "Pmax" = "#9467bd")  # purple

# ---- Helper: get first local start datetime from column A and convert first timestamp to local time ----
get_start_dt_local <- function(df, tz_local = "America/Edmonton") {
  # Try to find a non-empty timestamp string
  ts_chr <- as.character(df[[1]])
  ts_chr <- gsub("^[[:space:]\u00A0]+|[[:space:]\u00A0]+$", "", ts_chr, perl = TRUE)
  idx <- which(!is.na(ts_chr) & nzchar(ts_chr))
  if (!length(idx)) return(as.POSIXct(NA, tz = tz_local))
  
  s <- ts_chr[min(idx)]
  
  # Case 1: numeric timestamp (like POSIX seconds since epoch)
  if (suppressWarnings(!is.na(as.numeric(s)))) {
    t_utc <- as.POSIXct(as.numeric(s), origin = "1970-01-01", tz = "UTC")
    return(lubridate::with_tz(t_utc, tz_local))
  }
  
  # Case 2: ISO 8601 with Z or offset
  if (grepl("Z$", s) || grepl("[+-]\\d\\d:?\\d\\d$", s)) {
    t_utc <- suppressWarnings(lubridate::ymd_hms(s, tz = "UTC", quiet = TRUE))
    if (is.na(t_utc)) {
      t_utc <- suppressWarnings(as.POSIXct(s, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"))
    }
    return(lubridate::with_tz(t_utc, tz_local))
  }
  
  # Case 3: string without timezone (assume UTC, then convert)
  t_guess <- suppressWarnings(lubridate::ymd_hms(s, tz = "UTC", quiet = TRUE))
  if (is.na(t_guess))
    t_guess <- suppressWarnings(lubridate::parse_date_time(
      s,
      orders = c("Ymd HMS","mdy HMS","dmy HMS","Ymd HM","mdy HM","dmy HM"),
      tz = "UTC", quiet = TRUE
    ))
  
  if (is.na(t_guess)) return(as.POSIXct(NA, tz = tz_local))
  
  lubridate::with_tz(t_guess, tz_local)
}


# ---- Helper: pick a power column (prefers 'power', else 'power_filled', else fallback) ----
pick_power <- function(df, fallback_col_index = 7) {
  # candidate names in priority order (case-insensitive)
  candidates <- c("power", "power_filled", "p", "watts", "power_w", "power (w)", "power..w.")
  nms_low <- tolower(names(df))
  idx <- which(nms_low %in% candidates)
  
  source <- NULL
  if (length(idx) > 0) {
    # take the first match in priority order
    ord <- match(nms_low[idx], candidates)
    pick <- idx[order(ord)][1]
    vec  <- df[[pick]]
    source <- names(df)[pick]
  } else if (!is.na(fallback_col_index) && ncol(df) >= fallback_col_index) {
    vec    <- df[[fallback_col_index]]
    source <- paste0("col_", fallback_col_index, " (fallback)")
    warning(sprintf("No named power column found; using column %d as fallback.", fallback_col_index))
  } else {
    stop("No power-like column found and fallback index unavailable.")
  }
  
  # numeric coercion
  v <- suppressWarnings(readr::parse_number(as.character(vec)))
  
  # LOCF fill for any NA/invalid samples (keeps vector length unchanged)
  na_idx <- which(is.na(v))
  if (length(na_idx) > 0) {
    for (k in na_idx) v[k] <- if (k == 1) 0 else v[k - 1]
    message(sprintf("Filled %d missing/invalid power sample(s) using LOCF from '%s'.",
                    length(na_idx), source))
  }
  
  list(P = v, source = source)
}

# Compute Banister TRIMP with per-sample integration
compute_trimp <- function(hr_vec, dt_seconds, hr_rest, hr_max, gender = GENDER) {
  if (length(hr_vec) == 0 || all(is.na(hr_vec)) ||
      is.na(hr_rest) || is.na(hr_max) || hr_max <= hr_rest) return(NA_real_)
  # %HRR for each sample, clamped to [0,1]
  hrr <- (hr_vec - hr_rest) / (hr_max - hr_rest)
  hrr <- pmin(pmax(hrr, 0), 1)
  # Gender weighting
  k <- if (tolower(gender) == "male") 1.67 else 1.92
  y <- 0.64 * exp(k * hrr)
  # Sum over time (seconds -> minutes)
  sum((dt_seconds / 60) * hrr * y, na.rm = TRUE)
}

# ---- Normalized Power (NP) from second-by-second power ----
# Uses a 30-second rolling mean (time-aware via dt), ^4, mean, then 4th-root.
compute_np <- function(P, dt_seconds) {
  v <- suppressWarnings(as.numeric(P))
  if (!length(v) || all(is.na(v)) || !is.finite(dt_seconds) || dt_seconds <= 0)
    return(NA_real_)
  # number of samples in a 30s window
  k <- max(1L, round(30 / dt_seconds))
  # centered rolling mean; edges become NA (that’s OK for the later mean with na.rm=TRUE)
  rm30 <- as.numeric(stats::filter(v, rep(1/k, k), sides = 2))
  # NP definition
  m4 <- mean(rm30^4, na.rm = TRUE)
  if (!is.finite(m4)) return(NA_real_)
  m4^(1/4)
}

# ---- Helper: process one file ----
process_one <- function(path, CP, Wprime, Pmax, dt = 1, outdir = "Outputs", tz_local = "America/Edmonton") {
  df <- readr::read_csv(path, show_col_types = FALSE)
  
  # Skip files that have too few columns (malformed conversions)
  if (ncol(df) < 7) {
    cat("  ✖ Skipping", basename(path), "— only", ncol(df), "columns found (expected ≥7)\n")
    return(NULL)
  }
  
  # Parse start datetime (local) from column A
  start_dt <- get_start_dt_local(df, tz_local)
  
  start_date           <- if (!is.na(start_dt)) format(start_dt, "%Y-%m-%d") else NA
  start_time           <- if (!is.na(start_dt)) format(start_dt, "%H:%M:%S") else NA
  
  # Pick power column: prefers 'power', then 'power_filled', else falls back to column 7
  pow <- pick_power(df, fallback_col_index = 7)
  P <- pow$P
  power_source <- pow$source
  
  
  n <- length(P)
  t <- seq(0, by = dt, length.out = n)
  
  # Partition power into CP / W′ / Pmax channels
  PCP    <- pmin(P, CP)
  above  <- pmax(P - CP, 0)
  fracPM <- ifelse(P <= CP, 0, pmin(1, above / (Pmax - CP)))  # 0..1 between CP and Pmax
  PPmax  <- above * fracPM
  PWprim <- pmax(P - PCP - PPmax, 0)  # lactic (W′) channel
  
  # Skiba differential model using W′ expended:
  # if P > CP:  E_i = E_{i-1} + (P - CP)*dt
  # else:       E_i = E_{i-1} * exp(- (CP - P)*dt / W′)
  E <- numeric(n); E[1] <- 0
  for (i in 2:n) {
    if (P[i - 1] > CP) {
      E[i] <- E[i - 1] + (P[i - 1] - CP) * dt
    } else {
      E[i] <- E[i - 1] * exp(-(CP - P[i - 1]) * dt / Wprime)
    }
    if (E[i] < 0)       E[i] <- 0
    if (E[i] > Wprime)  E[i] <- Wprime
  }
  Wbal <- Wprime - E
  
  # Map W′bal -> MPA; guard against P == MPA
  MPA <- pmax(CP + (Pmax - CP) * (Wbal / Wprime), P + 1e-6)
  
  # Strain rates and Strain Scores
  k_strain <- (Pmax - MPA + CP) / (Pmax - P + CP)
  SR_total <- k_strain * P
  SR_CP    <- k_strain * PCP
  SR_Wp    <- k_strain * PWprim
  SR_Pmax  <- k_strain * PPmax
  
  # Normalize so 1 hour at CP ≈ 100 SS (when fully recovered)
  scale_factor <- (100/3600) * (Pmax / (CP^2))
  
  SS_total <- sum(SR_total) * scale_factor * dt
  SS_CP    <- sum(SR_CP)    * scale_factor * dt
  SS_Wp    <- sum(SR_Wp)    * scale_factor * dt
  SS_Pmax  <- sum(SR_Pmax)  * scale_factor * dt
  
  # Check whether recorded or estimated power was used
  power_source <- {
    if ("power" %in% names(df) && "power_filled" %in% names(df)) {
      same <- isTRUE(all.equal(df$power, df$power_filled, tolerance = 1e-8, check.attributes = FALSE))
      if (same) "estimated" else "recorded"
    } else if ("power" %in% names(df)) {
      "recorded"
    } else if ("power_filled" %in% names(df)) {
      "estimated"
    } else "unknown"
  }
  
  # Determine effective sample period (use your 'dt', or if you added auto-infer use that)
  dt_eff <- dt  # keep as-is unless you're using an inferred dt
  
  # TRIMP (Banister) using per-sample HR
  TRIMP <- if ("heart_rate" %in% names(df)) {
    # Coerce to numeric for safety
    hr_vec <- suppressWarnings(readr::parse_number(as.character(df$heart_rate)))
    compute_trimp(hr_vec, dt_eff, HR_REST, HR_MAX, gender = GENDER)
  } else {
    NA_real_
  }  
  # Effective sample period used throughout (use your existing dt; or dt_eff if you added it)
  dt_eff <- dt
  seconds_total <- n * dt_eff
  
  # ---- TSS components ----
  NP <- compute_np(P, dt_eff)
  IF <- if (is.finite(NP) && is.finite(CP) && CP > 0) NP / CP else NA_real_
  
  # TSS definition:
  # TSS = 100 * (seconds * NP * IF) / (CP * 3600)
  TSS <- if (is.finite(NP) && is.finite(IF) && is.finite(CP) && CP > 0) {
    100 * (seconds_total * NP * IF) / (CP * 3600)
  } else NA_real_
  
  # QA row
  qa <- tibble(
    file                 = basename(path),
    start_date           = start_date,
    start_time           = start_time,
    power_source         = power_source,
    duration_s           = n * dt,
    duration_h           = round(n * dt / 3600, 3),
    mean_power_W         = round(mean(P), 1),
    TRIMP                = round(TRIMP, 1), 
    NP_W                 = round(NP, 1), 
    pct_time_gt_CP       = round(100 * sum(P > CP) / n, 1),
    SS_total             = round(SS_total, 1),
    SS_CP                = round(SS_CP, 1),
    SS_Wprime            = round(SS_Wp, 1),
    SS_Pmax              = round(SS_Pmax, 1),
    IF                   = round(IF, 3),      
    TSS                  = round(TSS, 1),         
    min_Wbal_J           = round(min(Wbal), 0),
    final_Wbal_J         = round(tail(Wbal, 1), 0)
  )
  
  # Outputs
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  base <- tools::file_path_sans_ext(basename(path))
  
  # Per-second results CSV
  results <- tibble(
    time_s   = t,
    timestamp = if (!all(is.na(df[[1]]))) as.character(df[[1]]) else NA_character_,
    P        = P,
    PCP      = PCP,
    PWp      = PWprim,
    PPmax    = PPmax,
    Eexp_J   = E,
    Wbal_J   = Wbal,
    MPA_W    = MPA,
    k_strain = k_strain,
    SR_total = SR_total,
    SR_CP    = SR_CP,
    SR_Wp    = SR_Wp,
    SR_Pmax  = SR_Pmax
  )
  readr::write_csv(results, file.path(outdir, paste0(base, "_timeseries.csv")))
  readr::write_csv(qa,      file.path(outdir, paste0(base, "_QA.csv")))
  
  # ---- Stacked total bar (per ride) ----
  df_ss <- data.frame(
    component = factor(c("CP","W′","Pmax"), levels = c("CP","W′","Pmax")),
    SS = c(SS_CP, SS_Wp, SS_Pmax)
  )
  p_stacked <- ggplot(df_ss, aes(x = "Total", y = SS, fill = component)) +
    geom_col(width = 0.5) +
    scale_fill_manual(values = cols, name = NULL) +
    labs(x = NULL, y = "Strain Score (AU)",
         title = sprintf("%s — Total SS = %.1f", base, sum(df_ss$SS))) +
    theme_classic(base_size = 12)
  
  df_stack <- transform(df_ss,
                        y_mid = cumsum(SS) - SS/2,
                        pct = 100 * SS / sum(SS))
  p_stacked <- p_stacked +
    geom_text(data = df_stack,
              aes(x = "Total", y = y_mid, label = paste0(round(pct, 1), "%")),
              color = "white", size = 3.4)
  
  ggsave(file.path(outdir, paste0(base, "_SS_stacked.png")),
         p_stacked, width = 4.5, height = 4.5, dpi = 150)
  
  qa
}

# ---- Batch over folder ./Power files ----
files <- list.files(ride_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  warning(sprintf("No CSV files found in '%s'.", ride_dir))
} else {
  qa_all <- dplyr::bind_rows(lapply(
    files, function(f) process_one(f, CP, Wprime, Pmax, dt = dt, outdir = outdir, tz_local = tz_local)
  ))
  
  # Write aggregated QA (includes local date/time)
  readr::write_csv(qa_all, file.path(outdir, "strain_scores_QA_all.csv"))
  
  # Combined stacked composition chart across rides
  df_all <- qa_all |>
    select(file, SS_CP, SS_Wprime, SS_Pmax) |>
    pivot_longer(cols = starts_with("SS_"),
                 names_to = "component", values_to = "SS") |>
    mutate(component = recode(component,
                              SS_CP = "CP", SS_Wprime = "W′", SS_Pmax = "Pmax"),
           component = factor(component, levels = c("CP","W′","Pmax")))
  
  p_all <- ggplot(df_all, aes(x = file, y = SS, fill = component)) +
    geom_col() +
    scale_fill_manual(values = cols, name = NULL) +
    labs(x = "Ride", y = "Strain Score (AU)", title = "Strain Score Composition by Ride") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  width_in <- max(6, 0.2 * length(files) + 6)
  ggsave(file.path(outdir, "SS_composition_all.png"),
         p_all, width = width_in, height = 4.5, dpi = 150)
}

# =========================
# End of script
# =========================
