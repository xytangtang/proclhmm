rehmm_obj_fun <- function(paras, Y, N, K, quad_out) {

  para_mat <- inflate_paras_rehmm(paras, N, K)
  para_P1 <- para_mat$P1
  para_a <- para_mat$a
  para_b <- para_mat$b
  para_alpha <- para_mat$alpha
  para_beta <- para_mat$beta

  -compute_llh_rehmm(Y, para_a, para_b, para_alpha, para_beta, para_P1, quad_out$nodes, quad_out$weights)
}

rehmm_gr_fun <- function(paras, Y, N, K, quad_out) {
  para_mat <- inflate_paras_rehmm(paras, N, K)
  para_P1 <- para_mat$P1
  para_a <- para_mat$a
  para_b <- para_mat$b
  para_alpha <- para_mat$alpha
  para_beta <- para_mat$beta

  n <- length(Y)
  para_a_gr <- matrix(0, K, K-1)
  para_b_gr <- matrix(0, K, K-1)
  para_alpha_gr <- matrix(0, K, N-1)
  para_beta_gr <- matrix(0, K, N-1)
  para_P1_gr <- rep(0, K-1)

  for (index_seq in 1:n) {
    seq_gr_res <- seq2llh_gr_rehmm(Y[[index_seq]], para_a, para_b, para_alpha, para_beta, para_P1, quad_out$nodes, quad_out$weights)
    para_a_gr <- para_a_gr + seq_gr_res$para_a_gr
    para_b_gr <- para_b_gr + seq_gr_res$para_b_gr
    para_alpha_gr <- para_alpha_gr + seq_gr_res$para_alpha_gr
    para_beta_gr <- para_beta_gr + seq_gr_res$para_beta_gr
    para_P1_gr <- para_P1_gr + seq_gr_res$para_P1_gr
  }

  -deflate_paras_rehmm(para_a_gr, para_b_gr, para_alpha_gr, para_beta_gr, para_P1_gr)

}

#' Estimate latent traits in LHMM
#'
#' Compute MAP estimates of latent traits given LHMM parameters
#'
#' @param int_seqs a list of \code{n} action sequences where actions are coded as integers 0, ..., N-1
#' @param para_a \code{K} by \code{K-1} matrix. discrimination parameters of state transition probability matrix
#' @param para_b \code{K} by \code{K-1} matrix. location parameters of state transition probability matrix
#' @param para_alpha \code{K} by \code{N-1} matrix. discrimination parameters of state-action (emission) probability matrix
#' @param para_beta \code{K} by \code{N-1} matrix. location parameters of state-action (emission) probability matrix
#' @param para_P1 a vector of length \code{K-1}. parameters of initial state probability vector
#' @param n_pts number of quadrature points
#'
#' @return a vector of length \code{n}. Estimated latent traits.
#' @export
compute_theta <- function(int_seqs, para_a, para_b, para_alpha, para_beta, para_P1, n_pts) {
  n_seq <- length(int_seqs)
  N <- ncol(para_beta) + 1
  K <- nrow(para_a)

  out <- numeric(n_seq)

  quad_out <- statmod::gauss.quad(n_pts, kind="hermite")
  theta_list <- quad_out$nodes
  weight_list <- quad_out$weights

  P1 <- compute_P1_lhmm(para_P1)

  # this loop is computationally intensive
  for (index_seq in 1:n_seq) {
    d_list <- rep(0, n_pts)
    for (index_pt in 1:n_pts) {
      paras <- compute_PQ_lhmm(sqrt(2)*theta_list[index_pt], para_a, para_b, para_alpha, para_beta)
      d_list[index_pt] <- seq2llh(int_seqs[[index_seq]], P1, paras$P, paras$Q)
    }

    d_list <- d_list + log(weight_list)
    max_d <- max(d_list)
    theta_d_log <- max_d + log(sum(exp(d_list - max_d)))

    out[index_seq] <- sum(exp(d_list - theta_d_log) * sqrt(2) * theta_list)
  }

  out
}

#' MMLE of LHMM
#'
#' Maximum marginalized likelihood estimation of LHMM.
#' Marginalization over latent trait is computed numerically using Guassian-Hermite quadratures from \code{\link{statmod}}.
#' Optimization is performed through \code{\link{optim}}.
#'
#'
#' @param action_seqs a list of \code{n} action sequences
#' @param action_set set of possible actions; if \code{NULL}, derived as the set of distinct actions from \code{action_seqs}
#' @param K number of hidden states
#' @param paras a list of elements named \code{para_a}, \code{para_b}, \code{para_alpha}, \code{para_beta}, and \code{para_P1},
#'   providing initial values of model parameters
#' @param n_pts number of quadrature points
#' @param ... additional arguments passed to \code{\link{optim}}
#'
#' @return A list containing the following elements
#' \tabular{ll}{
#'   \code{seqs} \tab action sequences coded in integers \cr
#'   \tab \cr
#'   \code{K} \tab number of hidden states \cr
#'   \tab \cr
#'   \code{N} \tab number of distinct actions \cr
#'   \tab \cr
#'   \code{paras_init} \tab a list containing initial values of parameters \cr
#'   \tab \cr
#'   \code{paras_est} \tab a list containing parameter estimates \cr
#'   \tab \cr
#'   \code{theta_est} \tab a vector of length \code{n}. estimated latent traits \cr
#'   \tab \cr
#'   \code{init_mllh} \tab initial value of the marginalized likelihood function \cr
#'   \tab \cr
#'   \code{opt_mllh} \tab maximized marginalized likelihood function \cr
#'   \tab \cr
#'   \code{opt_res} \tab object returned by \code{\link{optim}} \cr
#'   }
#'
#' @export
#'

lhmm <- function(action_seqs, action_set, K, paras, n_pts, ...) {

  n <- length(action_seqs) # sample size

  # create action set
  if (is.null(action_set)) {
    action_set <- unique(unlist(action_seqs))
  }else {
    # print warning messages if provided non-null action_set does not match action_seqs
  }
  N <- length(action_set)
  action2id <- 1:N
  names(action2id) <- action_set

  # convert action sequences to integer sequences consisting of 0, 1, ..., N-1
  int_seqs <- list(n)
  for (i in 1:n) {
    seqt <- action_seqs[[i]]
    int_seqs[[i]] <- action2id[seqt] - 1
  }

  quad_out <- statmod::gauss.quad(n_pts, kind="hermite")

  # flatten parameters for optimization
  paras_init_rehmm <- deflate_paras_rehmm(paras$para_a, paras$para_b, paras$para_alpha, paras$para_beta, paras$para_P1)
  # objective function value at intial parameter values
  obj_val_init_rehmm <- rehmm_obj_fun(paras_init_rehmm, int_seqs, N, K, quad_out)

  cat("Optimizing obj fun...\n")
  # minimize the marginalized maximum likelihood function
  opt_res_rehmm <- optim(paras_init_rehmm, fn = rehmm_obj_fun, gr = rehmm_gr_fun, Y = int_seqs, N = N, K = K, quad_out = quad_out, method="BFGS", ...)

  # objective function value at optimized parameter values
  obj_val_opt_rehmm <- opt_res_rehmm$value

  # optimized parameter values
  opt_paras_rehmm <- inflate_paras_rehmm(opt_res_rehmm$par, N, K)

  # compute MAE of theta
  cat("Computing theta...\n")
  theta_est <- compute_theta(int_seqs, opt_paras_rehmm$a, opt_paras_rehmm$b, opt_paras_rehmm$alpha, opt_paras_rehmm$beta, opt_paras_rehmm$P1, n_pts)

  out <- list(seqs = int_seqs, K = K, N = N,
              paras_init = inflate_paras_rehmm(paras_init_rehmm, N, K),
              paras_est = opt_paras_rehmm,
              theta_est = theta_est,
              init_mllh = -obj_val_init_rehmm,
              opt_mllh = -obj_val_opt_rehmm,
              opt_res = opt_res_rehmm)

  return(out)
}
