## =====================================================================
##  analysis.R      (base R only, no packages needed)
##
##  Every number, table and figure input of the report and of its supplement.
##  Keep these files in one folder:
##    delhi_metro_lines.csv   the station table (the only hand-made input)
##    metro_network.R         network, distances, random networks
##    notes_tools.R           standard statistical formulas, written out by hand
##    random_networks_1000.rds  (optional) saved result of the 1,000 random
##                              networks; without it they are computed again
##    random_networks_1000_no_airport_express.rds  (optional) the same without
##                              the Airport Express (Section F)
## =====================================================================
source("metro_network.R")
source("notes_tools.R")


## ---- A. Data: station table -> network -> distances -> pair and station tables ----
Lraw <- read.csv("delhi_metro_lines.csv", stringsAsFactors = FALSE)         # as typed
invisible(capture.output(L <- read_lines_csv("delhi_metro_lines.csv")))     # cleaned and checked
net <- make_network(L)
ns  <- net$n                                    # number of stations
D   <- stops_matrix(net)                        # fewest stops, every pair of stations
RM  <- route_matrices(net)
I   <- RM$I                                     # fewest interchanges
Sx  <- RM$Sx                                    # fewest stops on a fewest-interchange route
ic  <- interchange_report(L)                    # stations where two or more coloured lines meet
colour   <- sub("-North$", "", L$group)         # the two Magenta sections count as one coloured line
observed <- tapply(L$station, colour, function(s) length(unique(s)))

inc    <- sapply(sort(unique(colour)), function(cn) as.integer(net$stations %in% L$station[colour == cn]))
shared <- (inc %*% t(inc)) > 0                  # do two stations share a coloured line?
p <- which(upper.tri(D), arr.ind = TRUE)        # each unordered pair once
P <- data.frame(pair_id = seq_len(nrow(p)),
                from_station = net$stations[p[, 1]], to_station = net$stations[p[, 2]],
                min_stops = D[p], min_interchanges = I[p],
                stops_at_min_interchanges = Sx[p], same_line = shared[p])
write.csv(P, "delhi_metro_pairs.csv", row.names = FALSE)

n_col <- tapply(colour, L$station, function(g) length(unique(g)))
stn <- data.frame(station = net$stations, stringsAsFactors = FALSE)
stn$lines          <- as.character(tapply(colour, L$station, function(g) paste(sort(unique(g)), collapse = " + "))[stn$station])
stn$interchange    <- as.integer(n_col[stn$station]) >= 2          # on two or more coloured lines
stn$centrality     <- rowSums(D) / (ns - 1)                        # mean stops to the other stations
stn$mean_ic        <- rowSums(I) / (ns - 1)                        # mean interchanges to the other stations
stn$to_interchange <- apply(D[, stn$station[stn$interchange], drop = FALSE], 1, min)   # stops to the nearest interchange
q_cent      <- quantile(stn$centrality, 0.25, names = FALSE)
stn$central <- stn$centrality <= q_cent                            # the most central quarter
write.csv(stn, "delhi_metro_stations.csv", row.names = FALSE)

x <- P$min_stops                                # variable of Questions 1 and 2
k <- P$min_interchanges                         # variable of Questions 2 and 3
N <- nrow(P)                                    # number of pairs (the population)


## ---- B. Question 1: distribution of the stops ----
mu    <- mean(x)                                # census mean (a parameter)
sigma <- sd_pop(x)                              # census SD, divisor N
S2    <- var_samp(x)                            # S^2, divisor N - 1 (population variance of sampling theory)
mom   <- sapply(1:4, function(r) central_moment(x, r))        # central moments mu_1 ... mu_4
qs    <- sapply(1:3, function(i) quartile_n1(x, i))               # quartiles by the (n + 1) rule
q7    <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)       # R's default rule, for comparison
iqr   <- qs[3] - qs[1]
mo    <- mode_value(x)
Fx    <- ecdf(x)
fence_hi <- qs[3] + 1.5 * iqr                   # Tukey's fences
fence_lo <- qs[1] - 1.5 * iqr
flag  <- x > fence_hi | x < fence_lo
x_in  <- x[!flag]
z     <- (x - mu) / sigma                       # z-scores
cover <- sapply(1:3, function(kk) mean(abs(z) <= kk))         # share within k SD (empirical rule)
cheb  <- 1 - 1 / (1:3)^2                                      # Chebyshev's guaranteed minimum
n_z3  <- sum(abs(z) > 3)
o_far <- order(x, decreasing = TRUE)[1:5]                     # the five farthest pairs
cnt_flag <- sort(table(c(P$from_station[flag], P$to_station[flag])), decreasing = TRUE)

## grouped data: classes of 5 stops (step-deviation mean and variance, grouped quartiles)
gd_tab  <- grouped_table(x, h = 5)
gd_A    <- 18                                                 # assumed mean: midpoint of the class 16-20
gd_mean <- grouped_mean(gd_tab, gd_A)
gd_var  <- grouped_var(gd_tab, gd_A)
gd_q    <- sapply(c(0.25, 0.5, 0.75), function(r) grouped_quantile(gd_tab, r))

## combined mean and variance: pairs inside the Yellow-Blue frame and all other pairs
yb    <- unique(L$station[L$group %in% c("Yellow", "Blue")])
in_yb <- P$from_station %in% yb & P$to_station %in% yb
x_yb  <- x[in_yb]
comb  <- combined_stats(c(sum(in_yb), sum(!in_yb)), c(mean(x[in_yb]), mean(x[!in_yb])),
                        c(var_pop(x[in_yb]), var_pop(x[!in_yb])))

## linear transformation: minutes if one stop took a = 2 minutes (an illustration, not DMRC data)
a_min <- 2
mins  <- a_min * x
trans <- c(mean = mean(mins), sd = sd_pop(mins), skew = skew_moment(mins), kurt = kurtosis_moment(mins),
           r_interchanges = cor(mins, k), r_reversed = cor(-x, k))

## ---- B2. Census against samples ----
set.seed(20261101)                                            # keep this line so results repeat
reps <- 2000
grid <- seq_len(max(x)); Fg <- Fx(grid)
ind2   <- as.numeric(k >= 2)                                  # pairs needing two or more interchanges
T_true <- sum(ind2)                                           # their number: a population total
pair_draw <- function(pool_idx, n) {                          # one SRSWOR of n pairs and its statistics
  ix <- pool_idx[sample.int(length(pool_idx), n)]
  xs <- x[ix]; ks <- k[ix]; o <- ols_by_hand(ks, xs)
  c(mean = mean(xs), median = median(xs), s2 = var_samp(xs), s2_n = var_pop(xs),
    gap = max(abs(ecdf(xs)(grid) - Fg)),                      # Kolmogorov-Smirnov distance to the census c.d.f.
    b1 = o$b1, b0 = o$b0, mse = o$s2, mse_n = o$SSE / n, se_b1 = sqrt(o$var_b1),
    total2 = N * mean(ks >= 2))
}
res <- list(`Random 30`  = t(replicate(reps, pair_draw(seq_len(N), 30))),
            `Random 100` = t(replicate(reps, pair_draw(seq_len(N), 100))),
            `Random 300` = t(replicate(reps, pair_draw(seq_len(N), 300))),
            `Convenience 300` = t(replicate(reps, pair_draw(which(in_yb), 300))))

cs   <- stn$centrality                                        # second finite population: the 243 stations
Ns   <- length(cs)
st_draw <- function(n, replace) { s <- cs[sample.int(Ns, n, replace = replace)]; c(mean = mean(s), median = median(s)) }
res_st <- list(`Stations 60 WOR`  = t(replicate(reps, st_draw(60, FALSE))),
               `Stations 120 WOR` = t(replicate(reps, st_draw(120, FALSE))),
               `Stations 60 WR`   = t(replicate(reps, st_draw(60, TRUE))),
               `Stations 120 WR`  = t(replicate(reps, st_draw(120, TRUE))))

mk <- function(scheme, r, n, Npop, truth, S2p, sig2p, wr = FALSE, gap = NA) {
  sdf <- if (wr) sqrt(srs_var_mean_wr(sig2p, n)) else sqrt(srs_var_mean(S2p, n, Npop))
  data.frame(scheme = scheme, n = n, f = n / Npop, mean_of_means = mean(r[, "mean"]),
             sd_obs = sd(r[, "mean"]), sd_formula = sdf, sd_median = sd(r[, "median"]),
             within1 = mean(abs(r[, "mean"] - truth) <= 1), gap = gap)
}
samp <- rbind(
  mk("Pairs, SRSWOR",  res[[1]], 30,  N, mu, S2, sigma^2, gap = median(res[[1]][, "gap"])),
  mk("Pairs, SRSWOR",  res[[2]], 100, N, mu, S2, sigma^2, gap = median(res[[2]][, "gap"])),
  mk("Pairs, SRSWOR",  res[[3]], 300, N, mu, S2, sigma^2, gap = median(res[[3]][, "gap"])),
  mk("Convenience (Yellow-Blue pairs)", res[[4]], 300, length(x_yb), mu, var_samp(x_yb), var_pop(x_yb),
     gap = median(res[[4]][, "gap"])),
  mk("Stations, SRSWOR", res_st[[1]], 60,  Ns, mean(cs), var_samp(cs), var_pop(cs)),
  mk("Stations, SRSWOR", res_st[[2]], 120, Ns, mean(cs), var_samp(cs), var_pop(cs)),
  mk("Stations, SRSWR",  res_st[[3]], 60,  Ns, mean(cs), var_samp(cs), var_pop(cs), wr = TRUE),
  mk("Stations, SRSWR",  res_st[[4]], 120, Ns, mean(cs), var_samp(cs), var_pop(cs), wr = TRUE))

bessel <- c(S2 = S2, sigma2 = sigma^2,                          # n = 30: divisor n - 1 against divisor n
            mean_s2 = mean(res[[1]][, "s2"]), mean_s2_n = mean(res[[1]][, "s2_n"]))
total  <- c(true = T_true, mean_est = mean(res[[3]][, "total2"]), sd_obs = sd(res[[3]][, "total2"]),
            sd_formula = sqrt(srs_var_total(var_samp(ind2), 300, N)))


## ---- C. Question 2: stops and interchanges ----
o2     <- ols_by_hand(k, x)                     # stops on interchanges, by the least squares formulas
fit2   <- lm(x ~ k)                             # the same line by lm(), as a check
lm_gap <- max(abs(coef(fit2) - c(o2$b0, o2$b1)))
r_s    <- cor(x, k, method = "spearman")
grp    <- split(x, k)                           # stops, grouped by interchanges
g_n    <- sapply(grp, length)
g_mean <- sapply(grp, mean)
g_var  <- sapply(grp, var_pop)
comb_k <- combined_stats(g_n, g_mean, g_var)    # combined mean and variance of the four groups: equals mu and sigma^2
eta2   <- sum(g_n * (g_mean - mu)^2) / o2$SST   # between-group share of the variation
res_sd <- tapply(o2$resid, k, sd)               # residual SD in each group (homoscedasticity)
o2t    <- ols_by_hand(k, sqrt(x))               # remedial measure: square-root transform
o2l    <- ols_by_hand(k, log(x))
slopes2 <- c(vertical = o2$b_yx, perpendicular = tls_slope(k, x), horizontal = o2$Syy / o2$Sxy)

reg_s <- sapply(c("Random 30", "Random 300"), function(nm) {  # the regression in 2,000 samples of n pairs
  r <- res[[nm]]
  c(n = as.numeric(sub("Random ", "", nm)), mean_b1 = mean(r[, "b1"]), sd_b1 = sd(r[, "b1"]),
    mean_se_b1 = mean(r[, "se_b1"]), mean_b0 = mean(r[, "b0"]), mean_mse = mean(r[, "mse"]),
    mean_mse_n = mean(r[, "mse_n"]))
})

xr     <- P$stops_at_min_interchanges           # stops on a fewest-interchange route
aff    <- xr > x                                # no single route achieves both minima
extra_stops <- (xr - x)[aff]
ex_i   <- which.max(xr - x)
route  <- c(n_aff = sum(aff), share = mean(aff), mean_extra = mean(extra_stops),
            median_extra = median(extra_stops), max_extra = max(extra_stops),
            r_route = cor(xr, k), mean_route = mean(xr))
same_but_change <- sum(P$same_line & k > 0)     # share a coloured line but need a change of train

cut1   <- median(x)
direct <- factor(ifelse(k == 0, "direct", "needs a change"), levels = c("direct", "needs a change"))
long   <- factor(ifelse(x > cut1, "long", "short"), levels = c("short", "long"))
ct2    <- table(direct, long)                   # cross tabulation
cs2    <- chisq_table(ct2)
OR1    <- odds_ratio(ct2)
p_long <- ct2[, "long"] / rowSums(ct2)
odds_chg <- ct2["needs a change", "long"] / ct2["needs a change", "short"]
odds_dir <- ct2["direct", "long"] / ct2["direct", "short"]
cuts   <- quantile(x, c(0.5, 0.75, 0.9), names = FALSE)
OR_s   <- sapply(cuts, function(ct_) odds_ratio(table(direct, factor(x > ct_, levels = c(FALSE, TRUE)))))
len3   <- cut(x, c(-Inf, qs[1], qs[3], Inf), labels = c("short", "medium", "long"))
ct43   <- table(interchanges = factor(k), length = len3)          # a 4 x 3 table
cs43   <- chisq_table(ct43)
p_long_k <- tapply(x > cut1, k, mean)           # share of long trips for 0, 1, 2, 3 interchanges


## ---- D. Question 3: centrality ----
ord     <- order(stn$centrality)
z_cent  <- (stn$centrality - mean(stn$centrality)) / sd_pop(stn$centrality)
o3      <- ols_by_hand(stn$to_interchange, stn$centrality)
fit3    <- lm(centrality ~ to_interchange, data = stn)
r3_s    <- cor(stn$centrality, stn$to_interchange, method = "spearman")
r_cm    <- cor(stn$centrality, stn$mean_ic, method = "spearman")
slopes3 <- c(vertical = o3$b_yx, perpendicular = tls_slope(stn$to_interchange, stn$centrality),
             horizontal = o3$Syy / o3$Sxy)
low     <- stn$to_interchange <= 3
high    <- stn$to_interchange >= 6
r3_low  <- cor(stn$centrality[low], stn$to_interchange[low])
r3_high <- cor(stn$centrality[high], stn$to_interchange[high])
res3_sd <- c(low = sd(o3$resid[low]), high = sd(o3$resid[high]))
far     <- which.max(stn$centrality)
far_ic  <- which(stn$interchange)[which.max(stn$centrality[stn$interchange])]
pick    <- order(o3$resid)[c(1:3, length(o3$resid) - 2:0)]

type3 <- factor(ifelse(stn$interchange, "interchange", "not an interchange"), levels = c("not an interchange", "interchange"))
cen3  <- factor(ifelse(stn$central, "central", "other"), levels = c("other", "central"))
tab3  <- table(type3, cen3)
cs3   <- chisq_table(tab3)
OR3   <- odds_ratio(tab3)
p_cen <- tab3[, "central"] / rowSums(tab3)
set.seed(20261103)                              # randomization reference: shuffle which stations are interchanges
B_perm   <- 10000
chi_perm <- replicate(B_perm, chisq_table(table(factor(sample(stn$interchange), c(FALSE, TRUE)),
                                                factor(stn$central, c(FALSE, TRUE))))$chi2)
p_perm   <- (1 + sum(chi_perm >= cs3$chi2 - 1e-9)) / (B_perm + 1)


## ---- E. Question 4: real network against 1,000 random networks ----
count_by_line <- function(Z) tapply(Z$station, Z$group, function(s) length(unique(s)))
run_experiment <- function(L, reps = 1000, seed = 1, nbins = 60) {
  set.seed(seed)
  ic0 <- interchange_report(L); ns0 <- length(unique(L$station))
  cbl0 <- count_by_line(L); meet0 <- sort(ic0$lines)
  extra <- function(Dx, icx) {                  # summaries kept for every network, real or random
    d_ic    <- apply(Dx[, icx, drop = FALSE], 1, min)         # stops to the nearest interchange
    cen     <- rowSums(Dx) / (ncol(Dx) - 1)                   # mean stops to the other stations
    central <- cen <= quantile(cen, 0.25, names = FALSE)
    t2 <- table(factor(colnames(Dx) %in% icx, c(FALSE, TRUE)), factor(central, c(FALSE, TRUE)))
    c(mean_to_interchange = mean(d_ic), max_to_interchange = max(d_ic),
      or_central = if (all(t2 > 0)) odds_ratio(t2) else NA_real_,
      slope_centrality = cov(d_ic, cen) / var(d_ic), r2_centrality = cor(d_ic, cen)^2,
      mean_between_ic = mean(Dx[icx, icx][upper.tri(Dx[icx, icx])]))   # mean stops between two interchange stations
  }
  net0 <- make_network(L); D0 <- stops_matrix(net0); I0 <- interchange_matrix(net0)
  real   <- c(summarise_network(D0, I0), extra(D0, ic0$station))
  sims   <- matrix(NA_real_, reps, length(real), dimnames = list(NULL, names(real)))
  inv    <- matrix(NA, reps, 4, dimnames = list(NULL, c("stations", "interchanges", "per_line", "meeting")))
  counts <- matrix(0L, reps, nbins)             # number of pairs at each number of stops
  for (r in seq_len(reps)) {                    # the random networks
    Lr    <- random_network(L)
    net_r <- make_network(Lr)
    D_r   <- stops_matrix(net_r)
    if (anyNA(D_r)) stop("A random network was disconnected (replication ", r, ").")
    if (max(D_r) > nbins) stop("Increase nbins.")
    I_r <- interchange_matrix(net_r); icr <- interchange_report(Lr)
    sims[r, ] <- c(summarise_network(D_r, I_r), extra(D_r, icr$station))
    inv[r, ]  <- c(length(unique(Lr$station)) == ns0, nrow(icr) == nrow(ic0),     # what must not change
                   identical(count_by_line(Lr), cbl0), identical(sort(icr$lines), meet0))
    counts[r, ] <- tabulate(D_r[upper.tri(D_r)], nbins = nbins)
    if (r %% 100 == 0) message("random networks done: ", r, " of ", reps)
  }
  list(real = real, sims = sims, inv = inv, counts = counts,
       real_counts = tabulate(D0[upper.tri(D0)], nbins = nbins),
       prop_at_or_below = vapply(names(real), function(m) mean(sims[, m] <= real[m], na.rm = TRUE), numeric(1)),
       seed = seed, reps = reps, version = 3L)
}

n_reps <- 1000
seed   <- 20261102
cache  <- "random_networks_1000.rds"            # the 1,000 random networks take several minutes, so they are saved
rnd    <- if (file.exists(cache)) readRDS(cache) else NULL
same   <- !is.null(rnd) && identical(rnd$version, 3L) && rnd$reps == n_reps && rnd$seed == seed &&
  identical(as.integer(rnd$real_counts), as.integer(tabulate(D[upper.tri(D)], nbins = length(rnd$real_counts))))
if (!same) {
  rnd <- run_experiment(L, reps = n_reps, seed = seed)
  saveRDS(rnd, cache)
}
real    <- rnd$real                             # the real network
sims    <- rnd$sims                             # one row per random network
prop    <- rnd$prop_at_or_below                 # share of random networks at or below the real value
nrep    <- rnd$reps
mc_tail <- vapply(names(real), function(m) (1 + sum(sims[, m] >= real[m], na.rm = TRUE)) / (nrep + 1), numeric(1))
mc_se   <- sqrt(prop * (1 - prop) / nrep)       # Monte Carlo standard error of a share (SE of a proportion)
inv_all <- all(rnd$inv)                         # the fixed quantities, checked in every one of the networks
## How much of the gap do longer arms explain? Across the 1,000 random networks, regress the mean stops on the
## mean stops to the nearest interchange, then apply the line to the real network's arms.
arm      <- ols_by_hand(sims[, "mean_to_interchange"], sims[, "mean_stops"])
arm_gap  <- unname(real["mean_stops"] - mean(sims[, "mean_stops"]))                              # the gap to explain
arm_part <- arm$b1 * unname(real["mean_to_interchange"] - mean(sims[, "mean_to_interchange"]))  # part due to longer arms
arm_up   <- unname(real["mean_stops"] - (arm$b0 + arm$b1 * real["mean_to_interchange"]))        # distance above the line
arm_z    <- arm_up / sd(arm$resid)                                                              # in residual SDs
## Breakdown of trip length: outer arms (stops to the nearest interchange) and core spacing (stops between two
## interchange stations). Regress the mean stops on both across the 1,000 random networks; because the fitted plane
## passes through the means, the predicted gap splits exactly into an arm part and a core part.
fit_q4    <- lm(mean_stops ~ mean_to_interchange + mean_between_ic, data = as.data.frame(sims))
b_q4      <- coef(fit_q4)
r2_q4     <- summary(fit_q4)$r.squared
gap_arm   <- unname(b_q4["mean_to_interchange"] * (real["mean_to_interchange"] - mean(sims[, "mean_to_interchange"])))
gap_core  <- unname(b_q4["mean_between_ic"] * (real["mean_between_ic"] - mean(sims[, "mean_between_ic"])))
pred_q4   <- unname(mean(sims[, "mean_stops"]) + gap_arm + gap_core)       # equals predict(fit_q4) at the real network
gap_left  <- unname(real["mean_stops"] - pred_q4)                          # part of the gap neither explains
n_core_ge <- sum(sims[, "mean_between_ic"] >= real[["mean_between_ic"]])   # random networks with spacing as wide
cum      <- t(apply(rnd$counts, 1, cumsum)) / sum(rnd$real_counts)   # c.d.f. of every random network
real_cum <- cumsum(rnd$real_counts) / sum(rnd$real_counts)           # c.d.f. of the real network
d_ic0  <- apply(D[, ic$station], 1, min)        # which stations are farthest from every interchange?
far_st <- names(d_ic0)[d_ic0 == max(d_ic0)]
far_ln <- unique(stn$lines[match(far_st, stn$station)])


## ---- F. Non-sampling error in the census: do our modelling decisions matter? ----
scenario <- function(Ls, label) {
  nt <- make_network(Ls); Dd <- stops_matrix(nt); Ii <- interchange_matrix(nt)
  d <- Dd[upper.tri(Dd)]; i <- Ii[upper.tri(Ii)]
  data.frame(Scenario = label, Stations = nt$n, Interchanges = nrow(interchange_report(Ls)),
             MeanStops = mean(d), MedianStops = median(d), Longest = max(d), MeanInterchanges = mean(i),
             stringsAsFactors = FALSE)
}
L_mc <- L; L_mc$station[L_mc$station == "Madhuban Chowk (Magenta)"] <- "Madhuban Chowk"   # one interchange station
L_ae <- L[L$group != "Airport-Express", ]                                                # no Airport Express
L_both <- L_ae; L_both$station[L_both$station == "Madhuban Chowk (Magenta)"] <- "Madhuban Chowk"
sens <- rbind(scenario(L,      "As modelled"),
              scenario(L_mc,   "Madhuban Chowk as one station"),
              scenario(L_ae,   "Without the Airport Express"),
              scenario(L_both, "Both changes"))

## Q4 without the Airport Express, like for like: the 1,000 random networks are drawn again from the network
## without it (same seed), so the real network is compared with random networks that also lack it. Saved, like
## the main experiment, because it takes several minutes.
D_ae     <- stops_matrix(make_network(L_ae))
cache_ae <- "random_networks_1000_no_airport_express.rds"
rnd_ae   <- if (file.exists(cache_ae)) readRDS(cache_ae) else NULL
same_ae  <- !is.null(rnd_ae) && identical(rnd_ae$version, 3L) && rnd_ae$reps == n_reps && rnd_ae$seed == seed &&
  identical(as.integer(rnd_ae$real_counts), as.integer(tabulate(D_ae[upper.tri(D_ae)], nbins = length(rnd_ae$real_counts))))
if (!same_ae) {
  rnd_ae <- run_experiment(L_ae, reps = n_reps, seed = seed)
  saveRDS(rnd_ae, cache_ae)
}
ae_real <- rnd_ae$real[["mean_stops"]]                       # real network without the Airport Express
ae_rand <- mean(rnd_ae$sims[, "mean_stops"])                 # average random network without it
ae_sd   <- sd(rnd_ae$sims[, "mean_stops"])
ae_gap  <- ae_real - ae_rand
ae_ge   <- sum(rnd_ae$sims[, "mean_stops"] >= ae_real)       # random networks with as many stops as the real one
ae_tail <- (1 + ae_ge) / (rnd_ae$reps + 1)
