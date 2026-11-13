## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----pivotal-fn---------------------------------------------------------------
random_pivotal_sample <- function(pik, n = NULL, prn = NULL, ...) {
  tol <- 1e-06
  active <- which(pik > tol & pik < 1 - tol)

  while (length(active) >= 2L) {
    ij <- sample(active, 2L)
    i <- ij[1L]
    j <- ij[2L]
    total <- pik[i] + pik[j]

    if (total < 1) {
      if (runif(1) < pik[i] / total)
        pik[c(i, j)] <- c(total, 0)
      else pik[c(i, j)] <- c(0, total)
    } else {
      if (runif(1) < (1 - pik[j]) / (2 - total)) {
        pik[c(i, j)] <- c(1, total - 1)
      } else {
        pik[c(i, j)] <- c(total - 1, 1)
      }
    }
    active <- which(pik > tol & pik < 1 - tol)
  }
  sort(which(pik > 0.5))
}

## ----pivotal-register---------------------------------------------------------
library(sondage)

register_method(
  "random_pivotal",
  type            = "wor",
  sample_fn       = random_pivotal_sample,
  joint_fn        = he_jip,
  fixed_size      = TRUE,
  variance_family = "pps_brewer",
  probabilities   = "exact"
)

## ----pivotal-use--------------------------------------------------------------
pik <- inclusion_prob(c(2, 3, 4, 5, 6, 7, 8, 9), n = 4)
s <- unequal_prob_wor(pik, method = "random_pivotal")
s

pikl <- joint_inclusion_prob(s)
round(pikl, 4)
round(sampling_cov(s, weighted = TRUE), 4)
joint_inclusion_prob(s, sampled_only = TRUE)

## ----pivotal-verify-----------------------------------------------------------
sim <- unequal_prob_wor(pik, method = "random_pivotal", nrep = 5000)
freq <- tabulate(sim$sample, nbins = length(pik)) / 5000
cbind(target = pik, empirical = freq)

N <- length(pik)
co_occur <- matrix(0, N, N)
for (j in seq_len(5000)) {
  selected <- sim$sample[, j]
  co_occur[selected, selected] <- co_occur[selected, selected] + 1
}
empirical_jip <- co_occur / 5000
he_pikl <- he_jip(pik)

pairs <- data.frame(
  i = c(1, 2, 3, 5),
  j = c(8, 7, 6, 8)
)
pairs$HE <- round(he_pikl[cbind(pairs$i, pairs$j)], 4)
pairs$empirical <- round(empirical_jip[cbind(pairs$i, pairs$j)], 4)
pairs

## ----pivotal-cleanup, include = FALSE-----------------------------------------
unregister_method("random_pivotal")

## ----tille, eval = requireNamespace("sampling", quietly = TRUE)---------------
tille_sample <- function(pik, n = NULL, prn = NULL, ...) {
  which(as.logical(sampling::UPtille(pik)))
}

tille_joint <- function(pik, sample_idx = NULL, ...) {
  pikl <- sampling::UPtillepi2(pik)
  if (!is.null(sample_idx)) {
    pikl <- pikl[sample_idx, sample_idx, drop = FALSE]
  }
  pikl
}

register_method(
  "tille",
  type            = "wor",
  sample_fn       = tille_sample,
  joint_fn        = tille_joint,
  fixed_size      = TRUE,
  variance_family = "pps_brewer",
  probabilities   = "exact"
)

pik <- inclusion_prob(c(2, 3, 4, 5, 6, 7, 8, 9), n = 4)
s <- unequal_prob_wor(pik, method = "tille")
s

# Exact joint inclusion probabilities from UPtillepi2
round(joint_inclusion_prob(s), 4)

# Full variance estimation chain
round(sampling_cov(s), 4)

## ----tille-cleanup, include = FALSE, eval = requireNamespace("sampling", quietly = TRUE)----
unregister_method("tille")

## ----custom-wr----------------------------------------------------------------
custom_multinomial_sample <- function(hits, n = NULL, prn = NULL, ...) {
  sample.int(length(hits), size = n, replace = TRUE, prob = hits)
}

custom_multinomial_joint <- function(hits, sample_idx = NULL, ...) {
  n <- round(sum(hits))
  keep <- if (is.null(sample_idx)) seq_along(hits) else sample_idx
  h <- hits[keep]
  factor <- (n - 1) / n
  joint <- factor * outer(h, h)
  diag(joint) <- h + factor * h^2
  joint
}

register_method(
  "custom_multinomial",
  type            = "wr",
  sample_fn       = custom_multinomial_sample,
  joint_fn        = custom_multinomial_joint,
  variance_family = "wr",
  probabilities   = "exact"
)

hits <- expected_hits(c(2, 3, 5, 10), n = 4)
s_wr <- unequal_prob_wr(hits, method = "custom_multinomial")
s_wr
joint_expected_hits(s_wr)
joint_expected_hits(s_wr, sampled_only = TRUE)

## ----custom-wr-cleanup, include = FALSE---------------------------------------
unregister_method("custom_multinomial")

## ----cube-lp, eval = requireNamespace("sampling", quietly = TRUE)-------------
cube_lp_sample <- function(pik, n = NULL, aux = NULL, ...) {
  X <- cbind(pik, aux)
  which(sampling::samplecube(X, pik, comment = FALSE) == 1)
}

register_method(
  "cube_lp",
  type            = "balanced",
  sample_fn       = cube_lp_sample,
  joint_fn        = he_jip,
  variance_family = "pps_brewer",
  probabilities   = "exact"
)

pik <- inclusion_prob(c(2, 3, 4, 5, 6, 7, 8, 9), n = 4)
x <- matrix(c(10, 20, 15, 25, 30, 35, 40, 45))
s <- balanced_wor(pik, aux = x, method = "cube_lp")
s

# Balancing check: HT estimate of the aux total vs the true total
colSums(x[s$sample, , drop = FALSE] / pik[s$sample]) - colSums(x)

## ----cube-lp-strata, eval = requireNamespace("sampling", quietly = TRUE)------
cube_lp_stratified <- function(pik, n = NULL, aux = NULL, strata = NULL, ...) {
  if (is.null(strata)) {
    return(cube_lp_sample(pik, n = n, aux = aux))
  }
  X <- if (is.null(aux)) matrix(pik, ncol = 1) else cbind(pik, aux)
  which(sampling::balancedstratification(X, strata, pik, comment = FALSE) == 1)
}

register_method(
  "cube_lp_str",
  type            = "balanced",
  sample_fn       = cube_lp_stratified,
  joint_fn        = he_jip,
  variance_family = "pps_brewer",
  supports_strata = TRUE,
  probabilities   = "exact"
)

pik <- rep(0.5, 8)
strata <- rep(1:2, each = 4)
s <- balanced_wor(pik, aux = matrix(as.double(1:8)), strata = strata,
                  method = "cube_lp_str")

# Within-stratum sample sizes are preserved
tabulate(strata[s$sample], nbins = 2)

## ----cube-lp-cleanup, include = FALSE, eval = requireNamespace("sampling", quietly = TRUE)----
unregister_method("cube_lp")
unregister_method("cube_lp_str")

## ----lpm1-fn------------------------------------------------------------------
lpm1_sample <- function(pik, n = NULL, aux = NULL, spread = NULL, ...) {
  if (is.null(spread)) {
    stop("method 'lpm1' requires 'spread'", call. = FALSE)
  }
  d <- as.matrix(dist(spread))
  diag(d) <- Inf
  p <- pik
  eps <- 1e-9
  repeat {
    u <- which(p > eps & p < 1 - eps)
    if (length(u) == 0L) {
      break
    }
    if (length(u) == 1L) {
      p[u] <- as.numeric(runif(1) < p[u])
      break
    }
    i <- u[sample.int(length(u), 1L)]
    v <- u[u != i]
    j <- v[which.min(d[i, v])]
    w <- u[u != j]
    if (w[which.min(d[j, w])] != i) {
      next # not mutual nearest neighbors: redraw i
    }
    s <- p[i] + p[j]
    if (s > 1) {
      if (runif(1) < (1 - p[j]) / (2 - s)) {
        p[i] <- 1
        p[j] <- s - 1
      } else {
        p[j] <- 1
        p[i] <- s - 1
      }
    } else {
      if (runif(1) < p[j] / s) {
        p[j] <- s
        p[i] <- 0
      } else {
        p[i] <- s
        p[j] <- 0
      }
    }
  }
  which(p > 1 - eps)
}

register_method(
  "lpm1",
  type             = "balanced",
  sample_fn        = lpm1_sample,
  variance_family  = "unsupported",
  supports_aux     = FALSE,
  supports_spread  = TRUE,
  probabilities    = "exact"
)

## ----lpm1-run-----------------------------------------------------------------
set.seed(25)
N <- 200
coords <- cbind(runif(N), runif(N))
pik <- rep(0.15, N)

s <- balanced_wor(pik, spread = coords, method = "lpm1")
s

# Spread diagnostic: mean distance to the nearest sampled neighbor
# (larger is better spread)
nn_dist <- function(idx) {
  d <- as.matrix(dist(coords[idx, ]))
  diag(d) <- Inf
  mean(apply(d, 1, min))
}
s2 <- balanced_wor(pik, spread = coords, method = "lpm2")
s3 <- balanced_wor(pik, spread = coords, method = "scps")
c(
  lpm1 = nn_dist(s$sample),
  lpm2 = nn_dist(s2$sample),
  scps = nn_dist(s3$sample),
  srs = nn_dist(sample.int(N, s$n))
)

## ----lpm1-cleanup, include = FALSE--------------------------------------------
unregister_method("lpm1")

## ----mc-joint, eval = FALSE---------------------------------------------------
# mc_joint <- function(pik, sample_idx = NULL, ..., B = 5000) {
#   N <- length(pik)
#   n <- as.integer(round(sum(pik)))
#   co <- matrix(0, N, N)
#   for (b in seq_len(B)) {
#     s <- my_sampler(pik, n = n)
#     co[s, s] <- co[s, s] + 1
#   }
#   # The diagonal of co / B estimates pik from the same draws. When the
#   # design gives pik exactly, diag(pikl) <- pik replaces the estimate.
#   pikl <- co / B
#   if (!is.null(sample_idx)) {
#     pikl <- pikl[sample_idx, sample_idx, drop = FALSE]
#   }
#   pikl
# }

