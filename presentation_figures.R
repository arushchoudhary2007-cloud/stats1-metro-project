## =====================================================================
##  presentation_figures.R      (base R only, no packages needed)
##
##  The figures of the presentation slides: the report's Figures 2, 3, 5
##  and 6 and the supplement's figure of sample means, redrawn with larger
##  text for slides. Every number comes from analysis.R.
##  Run from the project folder:  Rscript presentation_figures.R
##  Writes PNG files (300 dpi) into the folder presentation_figures/.
## =====================================================================
source("analysis.R")

out <- "presentation_figures"
dir.create(out, showWarnings = FALSE)
blue <- "#1f4e8c"; orange <- "#d9822b"; red <- "firebrick"
fmt  <- function(v, d = 2) formatC(v, format = "f", digits = d)

slide_png <- function(name, w, h) {             # w and h in inches, as placed on the slide
  png(file.path(out, name), width = w, height = h, units = "in", res = 300, type = "cairo",
      family = "Carlito", pointsize = 15)
}
base_par <- function() par(cex = 1, mgp = c(1.6, 0.5, 0), tcl = -0.3, las = 1,
                           cex.main = 1.05, cex.axis = 0.85, cex.lab = 0.95, font.main = 2)

## ---- Q1: shape of the minimum stops (report Figure 2) ----
slide_png("q1_shape.png", 8.2, 5.3)
layout(matrix(c(1, 2, 3, 4, 4, 4), 2, 3, byrow = TRUE), heights = c(2.15, 1))
base_par(); par(mai = c(0.75, 0.75, 0.4, 0.12))
plot(Fx, verticals = TRUE, do.points = FALSE, col = blue, lwd = 2.5,
     main = "Empirical c.d.f.", xlab = "minimum stops", ylab = "share of pairs")
abline(h = c(0.25, 0.5, 0.75), lty = 3, col = "grey60")
par(mai = c(0.75, 0.92, 0.4, 0.05))             # room for the four-digit counts
hist(x, breaks = seq(0.5, max(x) + 0.5, 1), col = "grey82", border = "white", main = "Histogram",
     xlab = "minimum stops", ylab = "", ylim = c(0, 1900))
title(ylab = "pairs", line = 2.55)
par(mai = c(0.75, 0.75, 0.4, 0.12))
abline(v = c(mo, median(x), mu), col = c(orange, "darkgreen", red), lwd = 2.5, lty = c(1, 1, 2))
legend("topright", c(paste("mode", mo), paste("median", median(x)), paste("mean", fmt(mu))),
       col = c(orange, "darkgreen", red), lwd = 2.5, lty = c(1, 1, 2), bty = "n", cex = 0.72, seg.len = 1.4)
qqnorm(x, pch = 16, cex = 0.35, col = blue, main = "Normal Q-Q plot",
       xlab = "Normal quantiles", ylab = "minimum stops")
qqline(x, col = red, lwd = 2)
par(mai = c(0.72, 0.3, 0.4, 0.3))
boxplot(x, horizontal = TRUE, col = "grey82", outpch = 16, outcex = 0.6, outcol = blue,
        xlab = "minimum stops", main = paste0("Box plot: ", paste(c(min(x), qs, max(x)), collapse = ", "),
                                              "; outliers beyond ", fence_hi, " stops (dashed)"))
abline(v = fence_hi, lty = 2, col = red, lwd = 2)
dev.off()

## ---- Q1: 2,000 sample means of each kind (supplement figure) ----
slide_png("q1_samples.png", 7.4, 4.7)
base_par(); par(mai = c(0.75, 2.05, 0.45, 0.25))
bx <- list(`Convenience, n = 300` = res[["Convenience 300"]][, "mean"],
           `Random, n = 300` = res[["Random 300"]][, "mean"], `Random, n = 100` = res[["Random 100"]][, "mean"],
           `Random, n = 30` = res[["Random 30"]][, "mean"])
plot.new(); plot.window(xlim = c(11.5, 21.5), ylim = c(0.5, 4.5), xaxs = "i")
rect(mu - 1, 0.5, mu + 1, 4.5, col = adjustcolor(red, alpha.f = 0.08), border = NA)   # within 1 stop of the census mean
boxplot(bx, horizontal = TRUE, add = TRUE, at = 1:4, axes = FALSE,
        col = c("#f6c999", "#c9d6ea", "#c9d6ea", "#c9d6ea"), border = c("#a35a12", blue, blue, blue), pch = 16, cex = 0.45)
axis(1, at = 12:21); axis(2, at = 1:4, labels = names(bx), tick = FALSE); box()
title(xlab = "sample mean (stops)")
abline(v = mu, lty = 2, lwd = 2, col = red)
text(mu, 4.78, paste0("census mean ", fmt(mu), "  (shaded: within 1 stop)"), col = red, cex = 0.8, xpd = NA)
text(mean(x[in_yb]), 1.62, fmt(mean(x[in_yb])), col = "#a35a12", cex = 0.85, font = 2)
dev.off()

## ---- Q2: stops against interchanges (report Figure 3) ----
slide_png("q2_stops_interchanges.png", 12.3, 3.55)
par(mfrow = c(1, 3)); base_par(); par(mai = c(0.75, 0.78, 0.4, 0.15))
set.seed(2)
plot(jitter(k, amount = 0.3), jitter(x, amount = 0.4), pch = 16, cex = 0.3,
     col = adjustcolor(blue, alpha.f = 0.12), xaxt = "n",
     xlab = "minimum interchanges", ylab = "minimum stops", main = "Stops against interchanges")
axis(1, at = 0:3)
abline(fit2, col = red, lwd = 2.5)
points(0:3, g_mean, pch = 18, cex = 2.2, col = "#f2a900")
legend("topleft", c("least squares line", "group mean"), col = c(red, "#f2a900"), lwd = c(2.5, NA),
       pch = c(NA, 18), pt.cex = 1.8, bty = "n", cex = 0.8)
boxplot(x ~ k, col = "grey82", pch = 16, cex = 0.35,
        xlab = "minimum interchanges", ylab = "minimum stops", main = "Box plots by interchanges")
bp <- barplot(100 * p_long_k, names.arg = 0:3, col = blue, border = "white", ylim = c(0, 100), yaxt = "n",
              xlab = "minimum interchanges", ylab = "% of long trips", main = "Share of long trips (> 16 stops)")
abline(h = 100 * mean(long == "long"), lty = 2, col = "grey35", lwd = 1.5); axis(2, at = c(0, 25, 50, 75, 100))
text(bp[1], 100 * mean(long == "long"), paste0("all pairs ", round(100 * mean(long == "long")), "%"), pos = 3, cex = 0.8, col = "grey25")
text(bp, 100 * p_long_k, labels = paste0(round(100 * p_long_k), "%"), pos = 1, offset = 0.5, cex = 0.9, col = "white", font = 2)
dev.off()

## ---- Q3: centrality (report Figure 5) ----
slide_png("q3_centrality.png", 8.3, 4.65)
par(mfrow = c(1, 2)); base_par(); par(mai = c(0.75, 0.75, 0.4, 0.12))
plot(stn$to_interchange, stn$centrality, pch = 16, cex = 0.85,
     col = ifelse(stn$interchange, orange, blue),
     xlab = "stops to the nearest interchange", ylab = "mean stops to the others",
     main = "What explains centrality?", ylim = c(min(stn$centrality) - 0.8, max(stn$centrality) + 1.5))
abline(fit3, col = red, lwd = 2.5)
abline(a = o3$ybar - slopes3[2] * o3$xbar, b = slopes3[2], col = "grey30", lty = 2, lwd = 1.5)
abline(a = o3$ybar - slopes3[3] * o3$xbar, b = slopes3[3], col = "grey30", lty = 3, lwd = 1.5)
legend("bottomright", legend = c("interchange", "other station", "least squares"), pch = c(16, 16, NA),
       lwd = c(NA, NA, 2.5), col = c(orange, blue, red), bty = "n", cex = 0.8)
line_st <- unique(data.frame(station = L$station, colour = colour))
cent_by_line <- split(stn$centrality[match(line_st$station, stn$station)], line_st$colour)
line_order   <- names(sort(sapply(cent_by_line, median)))
cent_by_line <- cent_by_line[rev(line_order)]
names(cent_by_line) <- sub("Airport-Express", "Airport Exp.", names(cent_by_line))
pal <- c(Red = "#d62728", Yellow = "#f2c200", Blue = "#4a90d9", Green = "#3a9d3a", Violet = "#8e5bbd",
         Pink = "#f58fb3", Magenta = "#c2185b", Grey = "grey65", `Airport Exp.` = "#f28c28")
par(mai = c(0.75, 1.15, 0.4, 0.12))
boxplot(cent_by_line, horizontal = TRUE, col = pal[names(cent_by_line)], pch = 16, cex = 0.6,
        xlab = "mean stops to the others", main = "Stations on each line")
abline(v = mean(stn$centrality), lty = 2, col = "grey40")
dev.off()

## ---- Q4: real network against 1,000 random networks (report Figure 6) ----
slide_png("q4_random_networks.png", 7.9, 5.3)
par(mfrow = c(2, 2)); base_par(); par(mai = c(0.75, 0.75, 0.4, 0.12))
hist_real <- function(mm, main, xlab) {
  v <- sims[, mm]; rng <- range(c(v, real[[mm]]))
  hist(v, breaks = 25, col = "grey82", border = "white", xlim = rng + c(-0.05, 0.05) * diff(rng),
       main = main, xlab = xlab, ylab = "random networks")
  abline(v = real[[mm]], col = red, lwd = 3)
}
hist_real("mean_stops", "Mean stops", "mean stops (red: real network)")
xs <- seq_along(real_cum)
plot(xs, real_cum, type = "n", xlim = c(0, 55), ylim = c(0, 1),
     xlab = "minimum stops", ylab = "share of pairs", main = "c.d.f. of stops")
for (r in seq_len(min(50, nrow(cum)))) lines(c(0, xs), c(0, cum[r, ]), type = "s", col = adjustcolor(orange, alpha.f = 0.35))
lines(c(0, xs), c(0, real_cum), type = "s", col = blue, lwd = 3)
legend("bottomright", legend = c("real network", "random networks"), col = c(blue, orange), lwd = c(3, 1.5), bty = "n", cex = 0.8)
qqnorm(sims[, "mean_stops"], pch = 16, cex = 0.4, col = blue, main = "Q-Q plot: random mean stops",
       xlab = "Normal quantiles", ylab = "mean stops")
qqline(sims[, "mean_stops"], col = red, lwd = 2)
hist_real("mean_changes", "Mean interchanges", "mean interchanges (red: real network)")
dev.off()

cat("Slide figures written to", normalizePath(out), "\n")
