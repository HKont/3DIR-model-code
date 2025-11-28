# ---- Libraries ----
suppressPackageStartupMessages({
  library(tools)        # file_path_sans_ext
  library(lubridate)    # ymd_hms
  library(FITfileR)     # reading FIT files
  library(data.table)   # rbindlist
  library(xml2)         # GPX parsing
  library(readr)        # used only for renaming block
})

#===========================
# Config
#===========================
fit_dir <- file.path(getwd(), "FIT files")
gpx_dir <- file.path(getwd(), "GPX files")

out_dir <- file.path(getwd(), "Power files")
dir.create(out_dir, showWarnings = FALSE)

#===========================
# Estimated power
#===========================

USE_FIXED_FORMULA <- FALSE

fixed_slope <- 1.96
fixed_intercept <- -107

calibration_csv <- "hr_power_calibration.csv"

cal <- read.csv(calibration_csv, check.names = FALSE)
HR_REST <- if ("HR rest" %in% names(cal)) cal$`HR rest`[1] else NA_real_
HR_MAX  <- if ("HR max"  %in% names(cal)) cal$`HR max`[1]  else NA_real_

if (is.na(HR_REST) || is.na(HR_MAX) || HR_MAX <= HR_REST) {
  warning("TRIMP: HR_rest/HR_max missing or invalid; TRIMP will be NA.")
}

build_power_estimator <- function() {
  if (USE_FIXED_FORMULA) {
    function(hr, hr_rest, hr_max) {
      hr2 <- pmin(pmax(hr, hr_rest, na.rm = TRUE), hr_max, na.rm = TRUE)
      p <- fixed_slope * hr2 + fixed_intercept
      p[p < 0] <- 0
      p
    }
  } else {
    cal <- read.csv(calibration_csv, check.names = FALSE)
    hr_rest <- if ("HR rest" %in% names(cal)) cal$`HR rest`[1] else NA_real_
    hr_max  <- if ("HR max"  %in% names(cal)) cal$`HR max`[1]  else NA_real_
    ss <- cal[, c("SteadyS HR","SteadyS Power")]
    ss <- ss[complete.cases(ss), ]
    names(ss) <- c("hr","power")
    if (nrow(ss) < 2) stop("Need at least two HR-power rows.")
    fit <- lm(power ~ hr, data = ss)
    slope <- coef(fit)[["hr"]]
    intercept <- coef(fit)[["(Intercept)"]]
    function(hr, hr_rest_in, hr_max_in) {
      hr_lo <- if (!is.na(hr_rest)) hr_rest else hr_rest_in
      hr_hi <- if (!is.na(hr_max))  hr_max  else hr_max_in
      hr2 <- hr
      if (!is.na(hr_lo)) hr2 <- pmax(hr2, hr_lo)
      if (!is.na(hr_hi)) hr2 <- pmin(hr2, hr_hi)
      p <- slope * hr2 + intercept
      p[p < 0] <- 0
      p
    }
  }
}

estimate_power_from_hr <- build_power_estimator()

#===========================
# Helpers 
#===========================
to_deg <- function(x) x * 180 / 2^31

zap_fit_invalids <- function(x) {
  if (!is.numeric(x)) return(x)
  bad_exact  <- c(255, 65535, 32767, 2147483647, -128, 127)
  x[x %in% bad_exact |
      abs(x - 42949672.95) < 1e-2 |
      abs(x - 21474836.47) < 1e-2] <- NA_real_
  x
}

haversine <- function(lat1, lon1, lat2, lon2) {
  R <- 6371000
  to_rad <- pi / 180
  phi1 <- lat1 * to_rad; phi2 <- lat2 * to_rad
  dphi <- (lat2 - lat1) * to_rad
  dlambda <- (lon2 - lon1) * to_rad
  a <- sin(dphi/2)^2 + cos(phi1)*cos(phi2)*sin(dlambda/2)^2
  2 * R * atan2(sqrt(a), sqrt(1 - a))
}

process_records_df <- function(out, out_path) {
  num_cols <- intersect(
    c("position_lat","position_long","heart_rate",
      "cadence","distance","power","temperature"),
    names(out)
  )
  out[num_cols] <- lapply(out[num_cols], function(v) suppressWarnings(as.numeric(v)))
  out[num_cols] <- lapply(out[num_cols], zap_fit_invalids)
  
  if ("timestamp" %in% names(out)) {
    if (!inherits(out$timestamp, "POSIXct")) {
      out$timestamp <- suppressWarnings(lubridate::ymd_hms(out$timestamp, quiet = TRUE))
    }
    out <- out[order(out$timestamp), ]
    out$elapsed_s <- as.numeric(difftime(out$timestamp, out$timestamp[1], units = "secs"))
  }
  
  if (!"distance" %in% names(out) &&
      all(c("position_lat","position_long") %in% names(out))) {
    n <- nrow(out)
    d_step <- numeric(n)
    if (n > 1) {
      for (i in 2:n) {
        d_step[i] <- haversine(
          out$position_lat[i-1], out$position_long[i-1],
          out$position_lat[i],   out$position_long[i]
        )
      }
    }
    out$distance <- cumsum(d_step)
  }
  
  n <- nrow(out)
  out$power_filled <- out$power
  
  if ("heart_rate" %in% names(out)) {
    est <- estimate_power_from_hr(out$heart_rate, HR_REST, HR_MAX)
    na_idx <- which(is.na(out$power) & !is.na(out$heart_rate))
    if (length(na_idx) > 0) {
      out$power_filled[na_idx] <- est[na_idx]
      out$power[na_idx] <- est[na_idx]
    }
  }
  
  if (all(is.na(out$power_filled))) {
    cat("  ✖ No usable power or HR in:", basename(out_path), "\n")
    return(FALSE)
  }
  
  final_cols <- intersect(
    c("timestamp","elapsed_s","position_lat","position_long",
      "heart_rate","cadence","distance","power","power_filled","temperature"),
    names(out)
  )
  
  write.csv(out[final_cols], out_path, row.names = FALSE)
  TRUE
}

read_gpx_as_df <- function(gpx_file) {
  x  <- xml2::read_xml(gpx_file)
  ns <- xml2::xml_ns(x)
  
  trkpts <- xml2::xml_find_all(x, ".//trkpt")
  if (length(trkpts) == 0) stop("No <trkpt> nodes found in GPX.")
  
  lat <- as.numeric(xml2::xml_attr(trkpts, "lat"))
  lon <- as.numeric(xml2::xml_attr(trkpts, "lon"))
  
  get_child_text <- function(nodes, xpath, ns = NULL) {
    vapply(nodes, function(n) {
      node <- xml2::xml_find_first(n, xpath, ns = ns)
      if (length(node) == 0) return(NA_character_)
      xml2::xml_text(node)
    }, character(1))
  }
  
  timestamp <- suppressWarnings(lubridate::ymd_hms(
    get_child_text(trkpts, ".//time")
  ))
  
  hr_txt <- get_child_text(trkpts, ".//gpxtpx:hr", ns)
  if (all(is.na(hr_txt))) hr_txt <- get_child_text(trkpts, ".//hr")
  heart_rate <- suppressWarnings(as.numeric(hr_txt))
  
  cad_txt <- get_child_text(trkpts, ".//gpxtpx:cad", ns)
  if (all(is.na(cad_txt))) cad_txt <- get_child_text(trkpts, ".//cad")
  cadence <- suppressWarnings(as.numeric(cad_txt))
  
  p_txt <- get_child_text(trkpts, ".//power")
  power <- suppressWarnings(as.numeric(p_txt))
  
  temp_txt <- get_child_text(trkpts, ".//gpxtpx:atemp", ns)
  if (all(is.na(temp_txt))) temp_txt <- get_child_text(trkpts, ".//atemp")
  temperature <- suppressWarnings(as.numeric(temp_txt))
  
  data.frame(
    timestamp     = timestamp,
    position_lat  = lat,
    position_long = lon,
    heart_rate    = heart_rate,
    cadence       = cadence,
    distance      = NA_real_,
    power         = power,
    temperature   = temperature,
    stringsAsFactors = FALSE
  )
}

#===========================
# FIT Processing
#===========================
fit_files <- list.files(path = fit_dir, pattern = "\\.fit$", full.names = TRUE)
for (f in fit_files) {
  cat("Processing FIT:", basename(f), "...\n")
  out_path <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(f)), "_record.csv"))
  
  tryCatch({
    rec <- records(readFitFile(f))
    if (is.list(rec) && !is.data.frame(rec)) rec <- rbindlist(rec, fill = TRUE)
    rec <- as.data.frame(rec)
    
    want <- c("timestamp","position_lat","position_long",
              "heart_rate","cadence","distance","power","temperature")
    keep <- intersect(want, names(rec))
    if (length(keep) == 0) stop("No expected record fields found.")
    
    out <- rec[, keep, drop = FALSE]
    
    ok <- process_records_df(out, out_path)
    if (ok) cat("  ✔ Wrote:", out_path, "\n")
  }, error = function(e) {
    cat("  ✖ Skipped due to error:", conditionMessage(e), "\n")
  })
}

#===========================
# GPX Processing
#===========================
gpx_files <- list.files(path = gpx_dir, pattern = "\\.gpx$", full.names = TRUE)

for (g in gpx_files) {
  cat("Processing GPX:", basename(g), "...\n")
  out_path <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(g)), "_record.csv"))
  
  tryCatch({
    out <- read_gpx_as_df(g)
    ok  <- process_records_df(out, out_path)
    if (ok) cat("  ✔ Wrote:", out_path, "\n")
  }, error = function(e) {
    cat("  ✖ Skipped due to error:", conditionMessage(e), "\n")
  })
}

cat("Conversion complete. Now renaming...\n")

#===========================
# CSV RENAMING BLOCK
#===========================
conv_dir <- out_dir
csv_files <- list.files(conv_dir, pattern = "\\.csv$", full.names = TRUE)

for (f in csv_files) {
  cat("Checking:", basename(f), "...\n")
  
  df <- tryCatch(
    readr::read_csv(f, show_col_types = FALSE,
                    col_types = readr::cols_only(timestamp = readr::col_character())),
    error = function(e) NULL
  )
  
  if (is.null(df) || !"timestamp" %in% names(df)) {
    cat("  ✖ No timestamp column, skipping.\n")
    next
  }
  
  ts <- suppressWarnings(lubridate::ymd_hms(df$timestamp[1], tz = "UTC"))
  
  if (is.na(ts)) {
    cat("  ✖ Invalid timestamp, skipping.\n")
    next
  }
  
  new_name <- paste0(format(ts, "%Y%m%d_%H%M"), ".csv")
  new_path <- file.path(conv_dir, new_name)
  
  if (file.exists(new_path)) {
    i <- 1
    repeat {
      candidate <- file.path(conv_dir, paste0(tools::file_path_sans_ext(new_name), "_", i, ".csv"))
      if (!file.exists(candidate)) {
        new_path <- candidate
        break
      }
      i <- i + 1
    }
  }
  
  file.rename(f, new_path)
  cat("  ✔ Renamed to:", basename(new_path), "\n")
}

cat("Done.\n")

# =========================
# End of script
# =========================
