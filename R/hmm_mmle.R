hmm_obj_fun <- function(paras, seqs, n, N, K) {
  para_list <- inflate_paras_hmm(paras, N, K)
  para_list <- compute_paras_hmm(para_list$para_P, para_list$para_Q, para_list$para_P1)
  P <- para_list$P
  Q <- para_list$Q
  P1 <- para_list$P1

  out <- compute_llh(seqs, n, P1, P, Q)

  -out

}

hmm_gr_fun <- function(paras, seqs, n, N, K) {

  para_list <- inflate_paras_hmm(paras, N, K)
  para_list <- compute_paras_hmm(para_list$para_P, para_list$para_Q, para_list$para_P1)
  P <- para_list$P
  Q <- para_list$Q
  P1 <- para_list$P1

  gr_res <- compute_llh_gr(seqs, P1, P, Q)
  para_P1_gr <- P1[-1] * (gr_res$P1_gr[-1] - sum(gr_res$P1_gr * P1))
  para_P_gr <- matrix(0, K, K-1)
  para_Q_gr <- matrix(0, K, N-1)

  for (k in 1:K) {
    para_P_gr[k, ] <- P[k, -1] * (gr_res$P_gr[k, -1] - sum(gr_res$P_gr[k,] * P[k, ]))
    para_Q_gr[k, ] <- Q[k, -1] * (gr_res$Q_gr[k, -1] - sum(gr_res$Q_gr[k,] * Q[k, ]))
  }

  -deflate_paras_hmm(para_P1_gr, para_P_gr, para_Q_gr)
}

#' MMLE of HMM
#'
#' Maximum marginalized likelihood estimation of HMM.
#' Optimization is performed through \code{\link{optim}}.

#'
#' @inheritParams lhmm
#' @param paras a list of elements named \code{para_P1}, \code{para_P}, and \code{para_Q},
#'   providing initial values of model parameters
#'
#' @return a list containing the following elements
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
#'   \code{init_mllh} \tab initial value of the marginalized likelihood function \cr
#'   \tab \cr
#'   \code{opt_mllh} \tab maximized marginalized likelihood function \cr
#'   \tab \cr
#'   \code{opt_res} \tab object returned by \code{\link{optim}} \cr
#'   }
#'
#' @examples
#' # generate data
#' paras_true <- sim_hmm_paras(5, 2)
#' sim_data <- sim_hmm(20, paras_true, 4, 10)
#' # randomly generate initial values of parameters
#' paras_init <- sim_hmm_paras(5, 2, return_prob=FALSE)
#' # fit hmm
#' hmm_res <- hmm(sim_data$seqs, 2, paras_init)
#'
#' @export

hmm <- function(action_seqs, K, paras, ...) {
  n <- length(action_seqs) # sample size

  # create action set
  action_set <- unique(unlist(action_seqs))

  N <- length(action_set)
  action2id <- 1:N
  names(action2id) <- action_set

  # convert action sequences to integer sequences consisting of 0, 1, ..., N-1
  int_seqs <- list(n)
  for (i in 1:n) {
    seqt <- action_seqs[[i]]
    int_seqs[[i]] <- action2id[seqt] - 1
  }

  paras_init_hmm <- compute_paras_hmm(paras$para_P, paras$para_Q, paras$para_P1)
  paras_vec <- deflate_paras_hmm(paras$para_P1, paras$para_P, paras$para_Q)
  obj_val_init_hmm <- hmm_obj_fun(paras_vec, int_seqs, n, N, K)

  cat("optimizing obj function...\n")
  opt_res_hmm <- optim(paras_vec, fn = hmm_obj_fun, gr = hmm_gr_fun, seqs = int_seqs, N = N, K = K, n = n, method="BFGS", ...)

  opt_paras_raw_hmm <- inflate_paras_hmm(opt_res_hmm$par, N, K)
  opt_paras_hmm <- compute_paras_hmm(opt_paras_raw_hmm$para_P, opt_paras_raw_hmm$para_Q, opt_paras_raw_hmm$para_P1)
  obj_val_opt_hmm <- opt_res_hmm$val

  colnames(opt_paras_hmm$P) <- rownames(opt_paras_hmm$P) <- paste0("state", 1:K)
  rownames(opt_paras_hmm$Q) <- paste0("state", 1:K)
  colnames(opt_paras_hmm$Q) <- action_set

  out <- list(seqs = int_seqs, K = K, N = N,
              paras_init = paras_init_hmm,
              paras_est = opt_paras_hmm,
              init_llh = -obj_val_init_hmm,
              opt_llh = -obj_val_opt_hmm,
              opt_res = opt_res_hmm)

  return(out)
}
