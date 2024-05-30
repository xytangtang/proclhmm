#' Compute LHMM probabilities from parameters
#'
#' Compute state-transition and state-action (emission) probability matrices
#' from LHMM parameters
#'
#' @param theta latent trait
#' @inheritParams compute_theta
#'
#' @return A list of two elements
#' \itemize{
#'   \item {\code{P}} {\code{K} by \code{K} state-transition probability matrix}
#'   \item {\code{Q}} {\code{K} by \code{N} state-action probability matrix}
#' }
#' @seealso \code{\link{compute_P1_lhmm}} for initial state probabilities of
#'   LHMM, \code{\link{compute_paras_hmm}} for computing probabilities in HMM.
#' @export
compute_PQ_lhmm <- function(theta, para_a, para_b, para_alpha, para_beta) {
	K <- nrow(para_a)
	N <- ncol(para_alpha) + 1

	a_mat <- cbind(rep(0, K), para_a)
	b_mat <- cbind(rep(0, K), para_b)
	alpha_mat <- cbind(rep(0, K), para_alpha)
	beta_mat <- cbind(rep(0, K), para_beta)

	P <- exp(a_mat * theta + b_mat)
	P <- P / rowSums(P)
	Q <- exp(alpha_mat * theta + beta_mat)
	Q <- Q / rowSums(Q)

	list(P = P, Q = Q)
}

#' Compute LHMM probabilities from parameters
#'
#' Compute initial state probability from LHMM parameters
#'
#' @inheritParams compute_theta
#' @return initial state probability vector of length \code{K}
#' @seealso \code{\link{compute_PQ_lhmm}} for state-transition and state-action
#'   probabilities of LHMM, \code{\link{compute_paras_hmm}} for computing
#'   probabilities in HMM.
#' @export
compute_P1_lhmm <- function(para_P1) {
	P1 <- c(0, para_P1)
	P1 <- exp(P1)
	P1 <- P1 / sum(P1)
	P1
}

# Gauss Hermite quadrature for calculating the expectation
# of a function of multivariate normal random vector
mgauss_hermite_quad <- function(n, mu, sigma) {
	gh <- gauss.quad(n, kind="hermite")
	dm <- length(mu)
	idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
	pts <- matrix(gh$nodes[idx],nrow(idx),dm)
	wts <- apply(matrix(gh$weights[idx],nrow(idx),dm), 1, prod)
	rot <- chol(sigma)
	pts <- t(t(pts %*% rot) + mu)

	return(list(pts = pts, wts = wts))
}

inflate_paras_rehmm <- function(paras, N, K) {
	para_a <- matrix(paras[1:(K*(K-1))], K, K-1)
	para_b <- matrix(paras[K*(K-1) + 1:(K*(K-1))], K, K-1)
	para_alpha <- matrix(paras[2*K*(K-1) + 1:(K*(N-1))], K, N-1)
	para_beta <- matrix(paras[2*K*(K-1) + K*(N-1) + 1:(K*(N-1))], K, N-1)
	para_P1 <- matrix(paras[2*K*(K-1) + 2*K*(N-1) + 1:(K-1)])

	list(para_a = para_a, para_b = para_b, para_alpha = para_alpha, para_beta = para_beta, para_P1 = para_P1)
}

deflate_paras_rehmm <- function(para_a, para_b, para_alpha, para_beta, para_P1) {
	c(as.vector(para_a), as.vector(para_b), as.vector(para_alpha), as.vector(para_beta), as.vector(para_P1))
}


#' Compute probabilities from logit scale parameters in HMM
#'
#' @param para_P \code{K} by \code{K-1} matrix. parameters of state-transition probability matrix
#' @param para_Q \code{K} by \code{N-1} matrix. parameters of state-action (emission) probability matrix
#' @param para_P1 \code{K-1} vector. parameters of initial state probability distribution
#' @return a list of three elements:
#' \itemize{
#'   \item {\code{P}} {\code{K} by \code{K} state-transition probability matrix}
#'   \item {\code{Q}} {\code{K} by \code{N} state-action (emission) probability matrix}
#'   \item {\code{P1}} {initial state probability vector of length \code{K}}
#' }
#' @seealso \code{\link{compute_PQ_lhmm}}, \code{\link{compute_P1_lhmm}} for computing probabilities in LHMM
#' @export
compute_paras_hmm <- function(para_P, para_Q, para_P1) {
	K <- nrow(para_P)
	P <- cbind(rep(0, K), para_P)
	P <- exp(P)
	P <- P / rowSums(P)

	Q <- cbind(rep(0, K), para_Q)
	Q <- exp(Q)
	Q <- Q / rowSums(Q)

	P1 <- c(0, para_P1)
	P1 <- exp(P1)
	P1 <- P1 / sum(P1)

	list(P = P, Q = Q, P1 = P1)

}

deflate_paras_hmm <- function(para_P1, para_P, para_Q) {
	c(as.vector(para_P), as.vector(para_Q), as.vector(para_P1))
}

inflate_paras_hmm <- function(paras, N, K) {
	para_P <- matrix(paras[1:(K*(K-1))], K, K-1)
	para_Q <- matrix(paras[K*(K-1) + 1:(K*(N-1))], K, N-1)
	para_P1 <- paras[K*(K-1) + K*(N-1) + 1:(K-1)]

	list(para_P = para_P, para_Q = para_Q, para_P1 = para_P1)
}



