## =====================================================================
##  metro_network.R      (base R only, no packages needed)
##
##  Minimum number of stops and minimum number of interchanges between
##  every pair of stations, and a comparison of the real network with
##  randomly rearranged networks (the "hypothetical population").
## =====================================================================

## ---- 1. Input format ------------------------------------------------
## The network is a table with ONE ROW PER STATION PER SERVICE, where a
## "line" below means a TRAIN SERVICE: a path that one train runs end to
## end without changing. A coloured metro line with a branch is therefore
## entered as several services that share the stations they have in common
## (e.g. "Blue-A": trunk + arm 1, "Blue-B": trunk + arm 2).
##
##   line     name of the service
##   pos      order of the station along that service: 1, 2, 3, ...
##   station  station name. A station served by several services (a trunk
##            station, or an interchange) must be spelt IDENTICALLY each time.
##   group    optional: the coloured line the service belongs to ("Blue").
##            Services in the same group share track, so stations they
##            share are NOT interchanges. If missing, group = line.
##   fixed    optional, TRUE/FALSE. TRUE = keep this station where it is in
##            the random networks. Mark ALL rows of such a station.
##   circular optional, TRUE/FALSE. TRUE for a service that runs round a loop
##            (the last station joins back to the first, e.g. the Pink Line).
##            List each station of the loop once.


## ---- 2. Reading and checking ---------------------------------------
read_lines_csv <- function(path) {
  L <- read.csv(path, stringsAsFactors = FALSE, strip.white = TRUE)
  L$station <- gsub("\\s+", " ", trimws(L$station))
  L$line    <- trimws(L$line)
  if (!"group" %in% names(L)) L$group <- L$line
  L$group <- trimws(L$group)
  if (!"fixed" %in% names(L)) L$fixed <- FALSE
  L$fixed <- as.logical(L$fixed)
  L$fixed[is.na(L$fixed)] <- FALSE
  if (!"circular" %in% names(L)) L$circular <- FALSE
  L$circular <- as.logical(L$circular)
  L$circular[is.na(L$circular)] <- FALSE
  L <- L[order(L$line, L$pos), ]
  rownames(L) <- NULL
  check_lines(L)
  L
}

check_lines <- function(L) {
  dup <- tapply(L$station, L$line, function(s) any(duplicated(s)))
  if (any(dup))
    stop("A station is repeated within a line: ",
         paste(names(dup)[dup], collapse = ", "))
  cat("Lines:", length(unique(L$line)),
      "| distinct stations:", length(unique(L$station)),
      "| interchange stations:", nrow(interchange_report(L)), "\n")
  invisible(TRUE)
}

## Stations shared by 2+ coloured lines = interchanges (compare with the real map).
interchange_report <- function(L) {
  if (!"group" %in% names(L)) L$group <- L$line
  gs <- tapply(L$group, L$station, function(x) paste(sort(unique(x)), collapse = " + "))
  gs <- gs[grepl("+", gs, fixed = TRUE)]
  data.frame(station = names(gs), lines = as.character(gs), row.names = NULL)
}

## Pairs of names that are almost identical (catches spelling differences,
## e.g. "Kashmere Gate" vs "Kashmiri Gate" that would hide an interchange).
possible_duplicates <- function(L, max_dist = 2) {
  u <- unique(L$station)
  d <- adist(tolower(u))
  d[lower.tri(d, diag = TRUE)] <- Inf
  w <- which(d <= max_dist, arr.ind = TRUE)
  if (nrow(w) == 0) return(data.frame())
  data.frame(a = u[w[, 1]], b = u[w[, 2]], distance = d[w])
}


## ---- 3. Building the network ---------------------------------------
make_network <- function(L) {
  L <- L[order(L$line, L$pos), ]
  stations <- unique(L$station)
  n        <- length(stations)
  sid      <- match(L$station, stations)   # station index of each row
  lines    <- unique(L$line)
  lid      <- match(L$line, lines)         # line index of each row

  circ <- if ("circular" %in% names(L)) tapply(L$circular, L$line, any)[lines] else rep(FALSE, length(lines))

  ## (a) station graph: neighbours = consecutive stations on a line
  from <- integer(0); to <- integer(0)
  for (l in seq_along(lines)) {
    idx <- which(lid == l)
    if (length(idx) > 1) {
      from <- c(from, sid[idx[-length(idx)]])
      to   <- c(to,   sid[idx[-1]])
      if (circ[l]) { from <- c(from, sid[idx[length(idx)]]); to <- c(to, sid[idx[1]]) }  # close the loop
    }
  }
  e   <- unique(cbind(pmin(from, to), pmax(from, to)))
  adj <- vector("list", n)
  for (i in seq_len(n)) adj[[i]] <- c(e[e[, 1] == i, 2], e[e[, 2] == i, 1])

  ## (b) "state" graph for interchanges: one state per (station, line)
  S <- length(sid)
  along <- vector("list", S)               # neighbours on the same line
  for (l in seq_along(lines)) {
    idx <- which(lid == l); m <- length(idx)
    if (m > 1) for (j in seq_len(m))
      along[[idx[j]]] <- c(if (j > 1) idx[j - 1], if (j < m) idx[j + 1])
    if (m > 2 && circ[l]) {                       # loop: last and first are neighbours
      along[[idx[1]]] <- c(along[[idx[1]]], idx[m])
      along[[idx[m]]] <- c(along[[idx[m]]], idx[1])
    }
  }
  stat_states <- lapply(seq_len(n), function(s) which(sid == s))
  paid <- vector("list", S)                 # same station, other service
  for (k in seq_len(S)) paid[[k]] <- setdiff(stat_states[[sid[k]]], k)
  list(stations = stations, n = n, adj = adj, lines = lines, S = S,
       st_station = sid, st_line = lid, along = along,
       paid = paid, stat_states = stat_states)
}


## ---- 4. Distances ---------------------------------------------------
## Minimum number of stops: breadth-first search from every station.
stops_matrix <- function(net) {
  n <- net$n
  D <- matrix(NA_integer_, n, n, dimnames = list(net$stations, net$stations))
  for (s in seq_len(n)) {
    d <- rep(NA_integer_, n); d[s] <- 0L
    queue <- integer(n); head <- 1L; tail <- 1L; queue[1] <- s
    while (head <= tail) {
      u <- queue[head]; head <- head + 1L
      for (v in net$adj[[u]]) if (is.na(d[v])) {
        d[v] <- d[u] + 1L; tail <- tail + 1L; queue[tail] <- v
      }
    }
    D[s, ] <- d
  }
  D
}

## Minimum number of interchanges: level-by-level search. Riding along a
## service costs 0; changing from one service to another at a station
## costs 1. Boarding any service at the origin is free.
interchange_matrix <- function(net) {
  n <- net$n; S <- net$S
  I <- matrix(NA_integer_, n, n, dimnames = list(net$stations, net$stations))
  for (s in seq_len(n)) {
    dist <- rep(NA_integer_, S)
    frontier <- net$stat_states[[s]]
    dist[frontier] <- 0L
    level <- 0L
    while (length(frontier)) {
      queue <- frontier                      # spread at cost 0 within this level
      while (length(queue)) {
        k <- queue[1]; queue <- queue[-1]
        for (j in net$along[[k]]) if (is.na(dist[j])) {
          dist[j] <- level; queue <- c(queue, j)
        }
      }
      nxt <- integer(0)                      # then one change of service
      for (k in which(dist == level)) for (j in net$paid[[k]]) if (is.na(dist[j])) {
        dist[j] <- level + 1L; nxt <- c(nxt, j)
      }
      frontier <- nxt; level <- level + 1L
    }
    I[s, ] <- vapply(net$stat_states, function(ix) {
      x <- dist[ix]
      if (all(is.na(x))) NA_integer_ else min(x, na.rm = TRUE)
    }, integer(1))
  }
  I
}


## ---- 5. Summaries ---------------------------------------------------
## Each unordered pair of stations counted once.
summarise_network <- function(D, I = NULL) {
  d <- D[upper.tri(D)]
  out <- c(mean_stops = mean(d), median_stops = median(d), max_stops = max(d))
  if (!is.null(I)) {
    i <- I[upper.tri(I)]
    out <- c(out, mean_changes = mean(i), share_2plus_changes = mean(i >= 2))
  }
  out
}

## Mean number of stops from each station to all the others (smaller = more central).
station_centrality <- function(D) sort(rowSums(D) / (ncol(D) - 1))


## ---- 6. Random networks (hypothetical population) ------------------
resample <- function(x, ...) x[sample.int(length(x), ...)]   # safe when length(x) == 1

## Same coloured lines, same number of stations on each, same interchanges
## (the same groups of coloured lines meet at each interchange station), but
## the position of every interchange along each of its coloured lines is
## re-drawn at random. Implemented by permuting station LABELS within each
## coloured line, so branches and shared trunks stay consistent.
random_network <- function(L) {
  if (!"group" %in% names(L)) L$group <- L$line
  fixed_st <- unique(L$station[L$fixed])
  ngroups  <- tapply(L$group, L$station, function(x) length(unique(x)))
  movable  <- setdiff(names(ngroups)[ngroups >= 2], fixed_st)
  R <- L
  for (g in unique(L$group)) {
    rows  <- which(L$group == g)
    slots <- setdiff(unique(L$station[rows]), fixed_st)  # labels allowed to move
    mv    <- slots[slots %in% movable]                   # interchange stations
    ord   <- slots[!(slots %in% movable)]                # ordinary stations
    new_label <- character(length(slots)); names(new_label) <- slots
    new_mv_slots <- resample(slots, length(mv))          # where interchanges go now
    new_label[new_mv_slots] <- resample(mv)
    rest <- setdiff(slots, new_mv_slots)
    if (length(rest)) new_label[rest] <- resample(ord)
    chg <- rows[L$station[rows] %in% slots]
    R$station[chg] <- new_label[L$station[chg]]
  }
  R
}

compare_random <- function(L, reps = 200,
                           with_changes = FALSE, seed = 1) {
  set.seed(seed)
  net <- make_network(L)
  D <- stops_matrix(net)
  if (anyNA(D)) stop("The real network is not connected - check station names.")
  real <- summarise_network(D, if (with_changes) interchange_matrix(net))
  sims <- matrix(NA_real_, reps, length(real), dimnames = list(NULL, names(real)))
  for (r in seq_len(reps)) {
    net_r <- make_network(random_network(L))
    D_r <- stops_matrix(net_r)
    if (anyNA(D_r)) stop("A random network was disconnected (replication ", r, ").")
    sims[r, ] <- summarise_network(D_r, if (with_changes) interchange_matrix(net_r))
  }
  ## proportion of random networks at or below the real value
  prop_below <- vapply(names(real), function(m) mean(sims[, m] <= real[m]), numeric(1))
  list(real = real, sims = sims, prop_at_or_below = prop_below)
}

plot_comparison <- function(res, measure = "mean_stops") {
  x <- res$sims[, measure]
  hist(x, breaks = 20, col = "grey85", border = "white",
       main = "Randomly rearranged networks", xlab = measure)
  abline(v = res$real[measure], lwd = 2, col = "red")
  legend("topright", legend = "real network", col = "red", lwd = 2, bty = "n")
}


## ---- 7. Optional: first draft of the line table from GTFS files ----
## Reads stops.txt, routes.txt, trips.txt, stop_times.txt from 'dir' and,
## for each route, uses its trip with the most stops as the station order.
## Routes often come once per direction, so CHECK the output: keep one
## entry per service (use route_ids = c(...)), and edit the 'group' column
## so that services on the same coloured line share a group name.
gtfs_to_lines <- function(dir, route_ids = NULL) {
  rd <- function(f) read.csv(file.path(dir, f), stringsAsFactors = FALSE, strip.white = TRUE)
  stops <- rd("stops.txt"); routes <- rd("routes.txt")
  trips <- rd("trips.txt"); st <- rd("stop_times.txt")
  tl <- table(st$trip_id)
  trips$n <- as.integer(tl[as.character(trips$trip_id)])
  trips <- trips[order(trips$route_id, -trips$n), ]
  best <- trips[!duplicated(trips$route_id), ]
  if (!is.null(route_ids)) best <- best[best$route_id %in% route_ids, ]
  out <- lapply(seq_len(nrow(best)), function(r) {
    x <- st[st$trip_id == best$trip_id[r], ]
    x <- x[order(x$stop_sequence), ]
    nm <- stops$stop_name[match(x$stop_id, stops$stop_id)]
    nm <- gsub("\\s+", " ", trimws(gsub("metro station", "", nm, ignore.case = TRUE)))
    lab <- routes$route_short_name[match(best$route_id[r], routes$route_id)]
    if (is.null(lab) || is.na(lab) || lab == "") lab <- as.character(best$route_id[r])
    data.frame(line = lab, pos = seq_along(nm), station = nm,
               group = lab, fixed = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
