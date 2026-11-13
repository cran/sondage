## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE------------------------------------------------------
library(sondage)
set.seed(1)

## -----------------------------------------------------------------------------
library(sondage)
s <- equal_prob_wor(N = 50, n = 5)
s

## -----------------------------------------------------------------------------
s$sample
s$pik[1:10]  # all equal to n/N = 0.1

## -----------------------------------------------------------------------------
s_sys <- equal_prob_wor(N = 50, n = 5, method = "systematic")
s_sys

## -----------------------------------------------------------------------------
s_ber <- equal_prob_wor(N = 50, n = 5, method = "bernoulli")
s_ber
length(s_ber$sample)  # may differ from 5

## -----------------------------------------------------------------------------
s_wr <- equal_prob_wr(N = 50, n = 5)
s_wr

## -----------------------------------------------------------------------------
s_wr$hits  # length-N vector of selection counts

## -----------------------------------------------------------------------------
x <- c(10, 20, 5, 40, 25)
n <- 3

# Without replacement: inclusion probabilities (capped at 1)
pik <- inclusion_prob(x, n)
pik
sum(pik)  # equals n

# With replacement: expected hits (can exceed 1)
eh <- expected_hits(x, n)
eh
sum(eh)   # equals n

## -----------------------------------------------------------------------------
s <- unequal_prob_wor(pik, method = "cps")
s
s$sample
inclusion_prob(s)  # extract the stored design-defining target pik

## -----------------------------------------------------------------------------
s_wr <- unequal_prob_wr(eh, method = "chromy")
s_wr
s_wr$hits  # realized selection counts

## -----------------------------------------------------------------------------
eh_big <- expected_hits(c(10, 20, 5, 40, 125), n)
eh_big
unequal_prob_wr(eh_big, method = "systematic_mr")$hits

## -----------------------------------------------------------------------------
pik <- inclusion_prob(c(10, 20, 5, 40, 25), n = 3)
x <- matrix(c(10, 20, 5, 40, 25))  # auxiliary variable to balance on
s_bal <- balanced_wor(pik, aux = x)
s_bal

## -----------------------------------------------------------------------------
N <- 20
pik_strat <- rep(0.4, N)
x_strat <- matrix(as.double(1:N), ncol = 1)
strata <- rep(1:4, each = 5)

s_strat <- balanced_wor(pik_strat, aux = x_strat, strata = strata)
s_strat

# Check per-stratum sample sizes
table(strata[s_strat$sample])

## -----------------------------------------------------------------------------
pik_control <- rep(0.5, 12)
group <- rep(c("a", "b", "c"), each = 4)
B <- sapply(unique(group), function(g) as.double(group == g))
expected_counts <- colSums(B * pik_control)

s_control <- balanced_wor(
  pik_control,
  bounds = list(
    B = B,
    lower = floor(expected_counts),
    upper = ceiling(expected_counts)
  )
)
colSums(B[s_control$sample, , drop = FALSE])

## -----------------------------------------------------------------------------
coordinates <- as.matrix(expand.grid(x = 1:4, y = 1:4))
pik_spatial <- rep(0.25, nrow(coordinates))

s_lpm2 <- balanced_wor(
  pik_spatial,
  spread = coordinates,
  method = "lpm2"
)
s_scps <- balanced_wor(
  pik_spatial,
  spread = coordinates,
  method = "scps"
)
s_lpm2$sample
s_scps$sample

## -----------------------------------------------------------------------------
set.seed(42)
burdened <- as.double(seq_len(16) %in% c(2, 7, 11))
u <- runif(16)
s1 <- balanced_wor(pik_spatial, spread = burdened, method = "scps", prn = u)
s2 <- balanced_wor(pik_spatial, spread = burdened, method = "scps", prn = 1 - u)
intersect(s1$sample, s2$sample)

## -----------------------------------------------------------------------------
pik <- inclusion_prob(c(10, 20, 5, 40, 25), n = 3)
s <- unequal_prob_wor(pik, method = "cps")

pikl <- joint_inclusion_prob(s)
pikl

## -----------------------------------------------------------------------------
eh <- expected_hits(c(10, 20, 5, 40, 25), n = 3)
s_wr <- unequal_prob_wr(eh, method = "chromy")

joint_expected_hits(s_wr)

## -----------------------------------------------------------------------------
pik_large <- inclusion_prob(runif(500), n = 10)
s_large <- unequal_prob_wor(pik_large, method = "brewer")

# Full 500 x 500 matrix
pikl_full <- joint_inclusion_prob(s_large)
dim(pikl_full)

# Only the 10 x 10 submatrix for sampled units
pikl_sub <- joint_inclusion_prob(s_large, sampled_only = TRUE)
dim(pikl_sub)

## -----------------------------------------------------------------------------
delta_sub <- sampling_cov(s_large, sampled_only = TRUE)
dim(delta_sub)

## -----------------------------------------------------------------------------
sampling_cov(s)

## -----------------------------------------------------------------------------
sampling_cov(s, weighted = TRUE)

## -----------------------------------------------------------------------------
sim <- equal_prob_wor(N = 50, n = 5, nrep = 4)
sim  # design object
dim(sim$sample)  # 5 x 4 matrix: each column is one sample

# Generics still work on batch objects
inclusion_prob(sim)

## ----echo = FALSE-------------------------------------------------------------
props <- data.frame(
  Method = c(
    "SRS WOR", "Systematic (EP)", "Bernoulli",
    "SRS WR",
    "CPS", "Brewer", "Sampford", "Systematic PPS", "Poisson",
    "SPS", "Pareto",
    "Chromy", "Multinomial", "Systematic PPS (MR)",
    "Cube", "LPM2", "SCPS"
  ),
  Type = c(
    rep("EP WOR", 3), "EP WR",
    rep("UP WOR", 7),
    rep("UP WR", 3),
    rep("Balanced", 3)
  ),
  `Fixed n` = c(
    "Yes", "Yes", "No",
    "Yes",
    "Yes", "Yes", "Yes", "Yes", "No",
    "Yes", "Yes",
    "Yes", "Yes", "Yes",
    "Yes", "Yes", "Yes"
  ),
  `Exact design marginals?` = c(
    "Yes", "Yes", "Yes",
    "Yes",
    "Yes", "Yes", "Yes", "Yes", "Yes",
    "Target only*", "Target only*",
    "Yes", "Yes", "Yes",
    "Yes", "Yes", "Yes"
  ),
  `Joint probs` = c(
    "Exact", "Exact (zeros)", "Exact",
    "Exact",
    "Exact", "HE approx", "Exact", "Exact (zeros)", "Exact",
    "HE approx", "HE approx",
    "Simulated", "Exact", "Exact (zeros)",
    "HE approx", "Not available", "Not available"
  ),
  check.names = FALSE
)
knitr::kable(props, align = "llccl",
  caption = paste(
    "HE approx = high-entropy approximation (Brewer & Donadio, 2003).",
    "For WOR methods, design marginals are first-order inclusion probabilities pi_k,",
    "and for WR methods, they are expected hits E(N_k).",
    "Target only* = the stored design-defining pik is the target vector.",
    "For sps and pareto the true finite-population first-order probabilities",
    "are only approximately equal to that target, with the discrepancy",
    "vanishing asymptotically.",
    "Simulated = Monte Carlo via nsim parameter."
  )
)

