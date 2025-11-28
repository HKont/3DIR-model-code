# =========================================================
# Mean Maximal Power (MMP) extractor
# - Reads .csv files from ./Power files
# - Expects a column named 'power_filled' (Watts)
# - Tries to infer elapsed time; handles irregular sampling by resampling to 1 Hz
# - Outputs ./MMP_summary.csv with one row per file
# =========================================================

# ---- User-configurable durations (seconds) ----
durations_sec <- c(5, 30, 45, 60, 75, 90, 120, 150, 180, 360, 480, 600, 720, 1200, 1800, 3600)

# ---- Libraries ----
suppressPackageStartupMessages({
  library(data.table)
})

# ---- Helpers ----

# infer an "elapsed seconds" vector from common time columns; returns numeric seconds
infer_elapsed_seconds <- function(DT) {
  # Candidate time columns (common names). Add your own if needed.
  cand <- c("elapsed_s", "elapsed_sec", "sec", "secs", "time_s", "t_s", "time_sec",
            "timestamp", "time", "datetime", "record_time", "start_time", "device_time")
  
  # Which exist?
  have <- intersect(cand, names(DT))
  if (length(have)) {
    # Try numeric first
    for (nm in have) {
      v <- DT[[nm]]
      # If it's already numeric, standardize to start at 0
      if (is.numeric(v)) {
        v <- as.numeric(v)
        if (all(is.finite(v))) {
          v <- v - v[1]
          return(v)
        }
      }
    }
    # Try parsing datetimes
    parse_try <- function(x) {
      # Try as POSIXct via fasttime or base
      xchr <- as.character(x)
      # Attempt ISO-like parsing
      z <- suppressWarnings(as.POSIXct(xchr, tz = "UTC"))
      if (all(is.na(z))) {
        # Try with fasttime if available
        if (requireNamespace("fasttime", quietly = TRUE)) {
          z <- suppressWarnings(fasttime::fastPOSIXct(xchr, tz = "UTC"))
        }
      }
      if (!all(is.na(z))) as.numeric(z) - as.numeric(z[1]) else NULL
    }
    for (nm in have) {
      v <- DT[[nm]]
      out <- parse_try(v)
      if (!is.null(out)) return(out)
    }
  }
  
  # If no time-like column, assume 1 Hz in original order
  seq_along(DT[[1]]) - 1L
}

# Resample to 1 Hz with linear interpolation
# Returns a list: list(t = integer seconds, p = numeric power)
resample_to_1hz <- function(t_sec, power) {
  # Clean inputs
  ok <- is.finite(t_sec) & is.finite(power)
  t_sec <- t_sec[ok]
  power <- power[ok]
  
  if (length(t_sec) < 2) return(NULL)
  
  # Ensure strictly increasing time (remove duplicates)
  ord <- order(t_sec)
  t_sec <- t_sec[ord]
  power <- power[ord]
  keep <- c(TRUE, diff(t_sec) > 0)
  t_sec <- t_sec[keep]
  power <- power[keep]
  
  if (length(t_sec) < 2) return(NULL)
  
  # 1 Hz grid
  t0 <- ceiling(t_sec[1])
  tN <- floor(t_sec[length(t_sec)])
  if (tN <= t0) return(NULL) # not enough span
  grid <- seq.int(t0, tN, by = 1L)
  
  # Linear interpolation
  p_interp <- stats::approx(x = t_sec, y = power, xout = grid, method = "linear",
                            rule = 2, ties = "ordered")$y
  
  list(t = grid, p = as.numeric(p_interp))
}

# Compute peak rolling mean (MMP) for integer-second window on 1 Hz data using cumsum trick
peak_mean_1hz <- function(p_1hz, window_sec) {
  w <- as.integer(window_sec)
  n <- length(p_1hz)
  if (!is.finite(w) || w < 1L || n < w) return(NA_real_)
  cs <- c(0, cumsum(p_1hz))
  # Means over windows ending at indices w..n
  means <- (cs[(w + 1):(n + 1)] - cs[1:(n - w + 1)]) / w
  max(means, na.rm = TRUE)
}

# ---- Main ----

# Locate input files
in_dir <- file.path(getwd(), "Power files")
if (!dir.exists(in_dir)) {
  stop("Input folder not found: ", in_dir)
}
csv_files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
  stop("No .csv files found in: ", in_dir)
}

# Prepare output container
out_cols <- c("file")
out_cols <- c(out_cols, paste0(durations_sec, "s"))
res_list <- vector("list", length(csv_files))

for (i in seq_along(csv_files)) {
  f <- csv_files[i]
  bn <- basename(f)
  
  # Read quickly with data.table
  DT <- tryCatch(fread(f, showProgress = FALSE), error = function(e) NULL)
  if (is.null(DT) || nrow(DT) == 0) {
    # Return all NAs for this file
    res <- as.list(c(bn, rep(NA_real_, length(durations_sec))))
    names(res) <- out_cols
    res_list[[i]] <- res
    next
  }
  
  # Ensure power column exists
  if (!("power_filled" %in% names(DT))) {
    # Try common alternatives (optional)
    alt <- intersect(c("power", "watts", "Power", "powerFilled"), names(DT))
    if (length(alt)) {
      setnames(DT, alt[1], "power_filled")
    } else {
      res <- as.list(c(bn, rep(NA_real_, length(durations_sec))))
      names(res) <- out_cols
      res_list[[i]] <- res
      next
    }
  }
  
  # Coerce power to numeric
  DT[, power_filled := suppressWarnings(as.numeric(power_filled))]
  
  # Infer time, resample to 1 Hz
  t_sec <- infer_elapsed_seconds(DT)
  rs <- resample_to_1hz(t_sec, DT[["power_filled"]])
  
  if (is.null(rs)) {
    res <- as.list(c(bn, rep(NA_real_, length(durations_sec))))
    names(res) <- out_cols
    res_list[[i]] <- res
    next
  }
  
  p1 <- rs$p
  
  # Compute MMPs
  mmp_vals <- vapply(durations_sec, function(d) peak_mean_1hz(p1, d), numeric(1))
  # Optional: round to 1 W
  mmp_vals <- round(mmp_vals, 1)
  
  res <- as.list(c(bn, mmp_vals))
  names(res) <- out_cols
  res_list[[i]] <- res
}

# Bind & write
OUT <- rbindlist(lapply(res_list, as.data.table), use.names = TRUE, fill = TRUE)

# Ensure columns are ordered correctly
setcolorder(OUT, out_cols)

# Write to disk
out_path <- file.path(getwd(), "MMP_summary.csv")
fwrite(OUT, out_path)

cat("Wrote:", out_path, "\n")

# =========================
# End of script
# =========================
