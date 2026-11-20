## =====================================================================
##  notes_tools.R      (base R only, no packages needed)
##
##  Standard statistical formulas written out one by one, so that every
##  number of the report can be traced to the formula that produced it.
## =====================================================================


## ---- Centre, spread, moments, skewness, kurtosis ------------------------
central_moment <- function(x, r) mean((x - mean(x))^r)               # r-th central moment mu_r, divisor N
var_pop        <- function(x) central_moment(x, 2)                   # population variance sigma^2, divisor N
sd_pop         <- function(x) sqrt(var_pop(x))                       # population SD sigma
var_samp       <- function(x) sum((x - mean(x))^2) / (length(x) - 1) # sample variance s^2, divisor n - 1
mean_abs_dev   <- function(x) mean(abs(x - mean(x)))                 # mean absolute deviation from the mean
median_abs_dev <- function(x) 1.4826 * median(abs(x - median(x)))    # robust; on the SD scale for Normal data
cv_pct         <- function(x) 100 * sd_pop(x) / mean(x)              # coefficient of variation, in percent
mode_value     <- function(x) { t <- table(x); as.numeric(names(t)[which.max(t)]) }

## Quartile i (i = 1, 2, 3) by the (n + 1) rule: the ordered value at position i/4 (n + 1)
## (the two ordered values around that position, averaged, when it is not a whole number).
quartile_n1 <- function(x, i) {
  xs <- sort(x); pos <- i / 4 * (length(xs) + 1)
  mean(xs[c(floor(pos), ceiling(pos))])
}

skew_moment     <- function(x) central_moment(x, 3) / central_moment(x, 2)^1.5   # moment coefficient gamma_1 = mu_3 / sigma^3
skew_pearson    <- function(x, mo = mode_value(x)) (mean(x) - mo) / sd_pop(x)       # Pearson's first coefficient: (mean - mode) / SD
skew_pearson2   <- function(x) 3 * (mean(x) - median(x)) / sd_pop(x)                # Pearson's second coefficient (mode ~ 3 median - 2 mean)
skew_bowley     <- function(q) (q[3] + q[1] - 2 * q[2]) / (q[3] - q[1])             # Bowley's quartile coefficient; q = c(Q1, Q2, Q3)
kurtosis_moment <- function(x) central_moment(x, 4) / central_moment(x, 2)^2        # beta_2 = mu_4 / sigma^4 (3 for a Normal curve)


## ---- Grouped data: step-deviation mean and variance, grouped quantiles ---
## Class table for integer data with classes of width h: real limits 0.5, 0.5 + h, ...
grouped_table <- function(x, h = 5, start = 0.5) {
  brk <- seq(start, by = h, length.out = ceiling((max(x) + 0.5 - start) / h) + 1)
  f   <- as.vector(table(cut(x, brk)))
  data.frame(lower = brk[-length(brk)], upper = brk[-1], mid = brk[-length(brk)] + h / 2,
             f = f, cf = cumsum(f))
}
grouped_mean <- function(g, A) {                    # step-deviation method: A + h sum(f u) / N, u = (midpoint - A) / h
  h <- g$upper[1] - g$lower[1]; u <- (g$mid - A) / h
  A + h * sum(g$f * u) / sum(g$f)
}
grouped_var <- function(g, A) {                     # h^2 [sum(f u^2) / N - (sum(f u) / N)^2]
  h <- g$upper[1] - g$lower[1]; u <- (g$mid - A) / h; N <- sum(g$f)
  h^2 * (sum(g$f * u^2) / N - (sum(g$f * u) / N)^2)
}
grouped_quantile <- function(g, r) {                # Q = L + ((N r - cf_prev) / f) h
  N <- sum(g$f); j <- which(g$cf >= N * r)[1]
  cfp <- if (j > 1) g$cf[j - 1] else 0
  g$lower[j] + ((N * r - cfp) / g$f[j]) * (g$upper[j] - g$lower[j])
}


## ---- Combined mean and variance of several groups -------------------------
## n, m, v: sizes, means and variances (divisor n) of the groups.
combined_stats <- function(n, m, v) {
  N <- sum(n); mc <- sum(n * m) / N; d <- m - mc
  c(mean = mc, var = sum(n * (v + d^2)) / N)
}


## ---- Simple linear regression by hand (least squares) ---------------------
ols_by_hand <- function(x, y) {
  n <- length(x); xb <- mean(x); yb <- mean(y)
  Sxx <- sum((x - xb)^2); Syy <- sum((y - yb)^2); Sxy <- sum((x - xb) * (y - yb))
  b1 <- Sxy / Sxx; b0 <- yb - b1 * xb                    # least squares estimates
  fitted <- b0 + b1 * x; e <- y - fitted
  SST <- Syy; SSE <- sum(e^2); SSR <- SST - SSE          # analysis of variance: SST = SSR + SSE
  s2 <- SSE / (n - 2)                                    # unbiased estimator of the error variance (divisor n - 2)
  k  <- (x - xb) / Sxx                                   # the constants k_i, so that b1 = sum(k_i y_i)
  list(n = n, xbar = xb, ybar = yb, Sxx = Sxx, Syy = Syy, Sxy = Sxy,
       b0 = b0, b1 = b1, fitted = fitted, resid = e,
       SST = SST, SSR = SSR, SSE = SSE, R2 = SSR / SST, s2 = s2, s = sqrt(s2),
       r = Sxy / sqrt(Sxx * Syy), b_yx = Sxy / Sxx, b_xy = Sxy / Syy,
       var_b1 = s2 / Sxx, var_b0 = s2 * (1 / n + xb^2 / Sxx), cov_b0b1 = -xb * s2 / Sxx,
       k = k, k_checks = c(sum_k = sum(k), sum_kx = sum(k * x), sum_k2_Sxx = sum(k^2) * Sxx),
       normal_eq = c(sum_e = sum(e), sum_ex = sum(e * x)))
}

## Slope of the line that minimises perpendicular distances (orthogonal regression: first principal component).
tls_slope <- function(x, y) {
  ev <- eigen(cov(cbind(x, y)))$vectors[, 1]
  ev[2] / ev[1]
}


## ---- Simple random sampling -----------------------------------------------
## S2 = sum (Y_i - mean)^2 / (N - 1), the population variance used in sampling theory.
srs_var_mean    <- function(S2, n, N) (1 - n / N) * S2 / n      # SRSWOR: variance of the sample mean, with finite-population correction
srs_var_mean_wr <- function(sig2, n) sig2 / n                   # SRSWR: no correction, sig2 = divisor-N variance
srs_var_total   <- function(S2, n, N) N^2 * srs_var_mean(S2, n, N)   # variance of the estimator N * ybar of a population total


## ---- Categorical data -------------------------------------------------------
odds_ratio  <- function(t2) (t2[2, 2] * t2[1, 1]) / (t2[2, 1] * t2[1, 2])   # rows (reference, group), columns (no, yes)
chisq_table <- function(tab) {                                              # expected counts, chi-square statistic, degrees of freedom
  E <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  list(expected = E, chi2 = sum((tab - E)^2 / E), df = (nrow(tab) - 1) * (ncol(tab) - 1))
}
