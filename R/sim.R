
#' generate HMM parameters
#'
#' \code{sim_hmm_paras} generates logit scale parameters of HMM with \code{K} hidden states and
#' \code{N} distinct actions from Uniform(-0.5, 0.5).
#'
#' @param N number of distinct actions
#' @param K number of hidden states
#' @param return_prob logical. indicates to return parameters in probability scale (\code{TRUE}, default) or logit scale.
#'
#' @return a list of three elements.
#'   If \code{return_prob = TRUE}, the element names are \code{P1}, \code{P}, and \code{Q}.
#'   If \code{return_prob = FALSE}, the element names are \code{para_P1}, \code{para_P}, and \code{oara_Q}.
#'
#' @examples
#' # generate probability parameters
#' set.seed(12345)
#' paras1 <- sim_hmm_paras(5, 2)
#' names(paras1)
#'
#' # generate parameters in the logit scale
#' set.seed(12345)
#' paras2 <- sim_hmm_paras(5, 2, return_prob = FALSE)
#' names(paras2)
#'
#' paras1$P1
#' paras2$para_P1
#'
#' # logit scale parameters can be transformed to probability parameters
#' all.equal(compute_paras_hmm(paras2$para_P, paras2$para_Q, paras2$para_P1), paras1)
#'
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
#' \code{sim_hmm_paras} generates the parameters of LHMM with \code{K} hidden states
#' and \code{N} distinct actions from Uniform(-0.5, 0.5).
#'
#' @inheritParams sim_hmm_paras
#'
#' @return a list of five elements, \code{para_a}, \code{para_b}, \code{para_alpha}, \code{para_beta}, and \code{para_P1}.
#'
#' @examples
#' paras <- sim_lhmm_paras(5, 2)
#' paras
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


#' Simulating action sequences using HMM
#'
#' \code{sim_hmm} generate \code{n} action sequences from HMM based on given parameters.
#' The lengths of the generated sequences are simulated from a Poission distribution with
#' mean \code{mean_len} and at least \code{min_len}.
#'
#' @param n number of action sequences to be generated
#' @param paras a list containing specified HMM parameters: state-transition probability matrix (\code{P}),
#' state-action probability matrix (\code{Q}), and initial state probability (\code{P1}).
#' @param min_len minimum length of generated sequences
#' @param mean_len mean length of generated sequences
#' @param return_state logical. Whether generated hidden state sequences should be returned or not.
#' @return \code{sim_hmm} returns a list of \code{n} generated action sequences if \code{return_state = FALSE}.
#' If \code{return_state = TRUE}, it returns a list of two lists, \code{seqs} and \code{state_seqs}. \code{seqs} gives
#' the generated action sequences. \code{state_seqs} gives the corresponding hidden state sequences.
#'
#' @examples
#' paras <- sim_hmm_paras(5,2)
#' sim_data <- sim_hmm(20, paras, 3, 10)
#'
#' @export
sim_hmm <- function(n, paras, min_len, mean_len, return_state = TRUE) {

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

#' Simulating action sequences using LHMM
#'
#' \code{sim_lhmm} generate \code{n} action sequences from LHMM based on given parameters.
#' The lengths of the generated sequences are simulated from a Poission distribution with
#' mean \code{mean_len} and at least \code{min_len}. The latent trait is generated from standard
#' normal.
#'
#' @param n number of action sequences to be generated
#' @param paras a list containing specified LHMM parameters: \code{para_a}, \code{para_b}, \code{para_alpha},
#' \code{para_beta}, and \code{para_P1}.
#' @param min_len minimum length of generated sequences
#' @param mean_len mean length of generated sequences
#' @param return_state logical. Whether generated hidden state sequences should be returned or not.
#' @return If \code{return_state = TRUE}, \code{sim_hmm} returns a list of three elements
#' \tabular{ll}{
#' \code{seqs} \tab a list of \code{n} generated action sequences \cr
#' \tab \cr
#' \code{theta} \tab latent traits as a vector of length \code{n} \cr
#' \tab \cr
#' \code{state_seqs} \tab a list of \code{n} hidden state sequences
#' }
#' If \code{return_state = FALSE}, the returned list only contains \code{seqs} and \code{theta}.
#'
#' @examples
#' paras <- sim_lhmm_paras(5,2)
#' sim_data <- sim_lhmm(20, paras, 4, 10)
#'
#' @export
sim_lhmm <- function(n, paras, min_len, mean_len, return_state = TRUE) {
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

  out <- list(seqs = seqs, theta = theta)
  if (return_state) out <- list(seqs = seqs, theta = theta, state_seqs = state_seqs)

  out

}
