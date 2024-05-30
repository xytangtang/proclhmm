
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proclhmm: Latent Hidden Markov Model for Process Data

<!-- badges: start -->

<!-- badges: end -->

This package provides functions for simulating from and fitting the
latent hidden Markov models for response process data
([Tang, 2024](https://doi.org/10.1007/s11336-023-09938-1)). It also
includes functions for simulating from and fitting ordinary hidden
Markov models.

## Installation

You can install the development version of `proclhmm` from GitHub with

``` r
devtools::install_github("xytangtang/proclhmm")
```

## Basic Usage

### Simulate Parameters and Data

``` r
library(proclhmm)
N <- 10 # number of actions
K <- 3 # number of hidden states

# generate parameters
set.seed(12345)
paras_true <- sim_lhmm_paras(N, K)

n <- 100 # sample size
# generate data
data0 <- sim_lhmm(n, paras_true, min_len = 4, mean_len = 25)
action_seqs <- data0$seqs # action sequences
```

### Estimating parameters and latent traits

``` r
# generate initial values of parameters
paras_init <- sim_lhmm_paras(N, K)
# model fitting
lhmm_res <- lhmm(action_seqs, K, paras_init)
#> Optimizing obj fun...
#> Computing theta...
# estimated discrimation parameters for state transition probability matrix
lhmm_res$paras_est$para_a
#>            state2       state3
#> state1 -0.6088158 -15.52584332
#> state2 -2.3456255  -2.83225548
#> state3 -4.5802817   0.08433227
# estimated location parameters for state-action probability matrix
lhmm_res$paras_est$para_beta
#>                 4          0          5         9         6         7
#> state1 -0.9987221  4.4206005  2.3740381 -1.472788 -6.503268  4.059780
#> state2  0.5821316 -1.4579551  0.5479029  1.214155  1.044911  0.736849
#> state3 -3.2111616 -0.5614971 -0.3598079 -6.618605 -4.027255 -1.426849
#>                  1          2         8
#> state1  3.24121223  4.1195709 -5.666249
#> state2 -0.32826612  0.2443809 -8.464710
#> state3 -0.05355758 -1.0904959 -9.719579
```

### Most likely hidden state sequences

``` r
# compute state-transition and state-action probability matrix for the first action sequnece
paras_est <- lhmm_res$paras_est
paras_PQ <- compute_PQ_lhmm(lhmm_res$theta_est[1], paras_est$para_a, paras_est$para_b, paras_est$para_alpha, paras_est$para_beta)
P <- paras_PQ$P
Q <- paras_PQ$Q
# compute initial state probability
P1 <- compute_P1_lhmm(paras_est$para_P1)
# find the most likely hidden state sequences for the first action sequence
find_state_seq(action_seqs[[1]], P1, P, Q)
#>  [1] 2 1 2 1 2 1 1 2 1 2 1 1 2 1 1 2 1 2 1 2 1 2 1 1 2 1
```
