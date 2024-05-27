
#' generate HMM parameters
#'
#' generate logit scale parameters from Uniform(-0.5, 0.5)
#'
#' @param N number of distinct actions
#' @param K number of hidden states
#' @param return_prob logical. indicates to return parameters in probability scale (\code TRUE, default) or logit scale.
#'
#' @return a list of three elements.
#'   If \code{return_prob = TRUE}, the element names are \code{P1}, \code{P}, and \code{Q}.
#'   If \code{return_prob = FALSE}, the element names are \code{para_P1}, \code{para_P}, and \code{oara_Q}.
#' @export
sim_hmm_paras <- function(N, K, return_prob = TRUE) {

  para_P1 <- runif(K-1) - 0.5
  para_P <- matrix(runif(K*(K-1))-0.5, K, K-1)
  para_Q <- matrix(runif(K*(N-1))-0.5, K, N-1)

  out <- list(para_P1 = para_P1, para_P = para_P, para_Q = para_Q)

  if (return_prob) out <- compute_paras_hmm(para_P, para_Q, para_P1)

  return(out)
}

#' generate LHMM parameters
#'
#' The parameters are generated from Uniform(-0.5, 0.5)
#'
#' @inheritParams sim_hmm_paras
#'
#' @return a list of five elements, \code{para_a}, \code{para_b}, \code{para_alpha}, \code{para_beta}, and \code{para_P1}.
#'
#' @export
sim_lhmm_paras <- function(N, K) {

  para_P1 <- runif(K-1) - 0.5
  para_a <- matrix(runif(K*(K-1))-0.5, K, K-1)
  para_b <- matrix(runif(K*(K-1))-0.5, K, K-1)
  para_alpha <- matrix(runif(K*(N-1))-0.5, K, N-1)
  para_beta <- matrix(runif(K*(N-1))-0.5, K, N-1)

  return(list(para_a = para_a, para_b = para_b, para_alpha = para_alpha, para_beta = para_beta, para_P1 = para_P1))

}



#' @export
sim_hmm <- function(n, paras, min_len, mean_len, return_hidden_state = TRUE) {

  K <- nrow(paras$P)
  N <- ncol(paras$Q)
  P <- paras$P
  Q <- paras$Q
  P1 <- paras$P1

  seqs <- list()
  state_seqs <- list()
  l_list <- rpois(n, mean_len)
  l_list[l_list < min_len] <- min_len

  for (i in 1:n) {
    l <- l_list[[i]]
    action_seq <- numeric(l)
    state_seq <- numeric(l)
    state_seq[1] <- sample(1:K, 1, prob = P1)
    action_seq[1] <- sample(1:N, 1, prob = Q[state_seq[1],])
    for (j in 2:l) {
      state_seq[j] <- sample(1:K, 1, prob = P[state_seq[j-1], ])
      action_seq[j] <- sample(1:N, 1, prob = Q[state_seq[j], ])
    }
    seqs[[i]] <- action_seq - 1
    state_seqs[[i]] <- state_seq
  }

  out <- seqs
  if (return_state) out <- list(seqs = seqs, state_seqs = state_seqs)

  out
}

#' @export
sim_lhmm <- function(n, paras, min_length, mean_len, return_hidden_state = TRUE) {
  K <- nrow(paras$para_a)
  N <- ncol(paras$para_alpha) + 1

  seqs <- list()
  state_seqs <- list()
  l_list <- rpois(n, mean_len)
  l_list[l_list < min_len] <- min_len

  theta <- rnorm(n)

  for (i in 1:n) {
    paras_mat <- compute_PQ_lhmm(theta[i], paras$para_a, paras$para_b, paras$para_alpha, paras$para_beta)
    P <- paras_mat$P
    Q <- paras_mat$Q
    P1 <- compute_P1_lhmm(paras$para_P1)
    l <- l_list[[i]]
    action_seq <- numeric(l)
    state_seq <- numeric(l)
    state_seq[1] <- sample(1:K, 1, prob = P1)
    action_seq[1] <- sample(1:N, 1, prob = Q[state_seq[1],])
    for (j in 2:l) {
      state_seq[j] <- sample(1:K, 1, prob = P[state_seq[j-1], ])
      action_seq[j] <- sample(1:N, 1, prob = Q[state_seq[j], ])
    }
    seqs[[i]] <- action_seq - 1
    state_seqs[[i]] <- state_seq
  }

  out <- seqs
  if (return_state) out <- list(seqs = seqs, state_seqs = state_seqs)

  out

}
