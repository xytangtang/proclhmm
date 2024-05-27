#' proclhmm: Latent Hidden Markov Models for Response Process Data
#'
#' This package provides functions for simulating from and fitting the latent hidden
#' Markov models for response process data (Tang, 2024).
#' It also includes functions for simulating from and fitting ordinary hidden Markov models.
#'
#' @section Data Simulation Functions:
#' \itemize{
#'   \item \code{\link{sim_hmm_paras}} generates parameters of HMM
#'   \item \code{\link{sim_hmm}} generates actions sequences from HMM.
#'   \item \code{\link{sim_lhmm_paras}} generates parameters of LHMM
#'   \item \code{\link{sim_lhmm}} generates actions sequences from LHMM.
#' }
#'
#' @section Model Fitting Functions:
#' \itemize{
#'   \item \code{\link{hmm}} fits HMM models. Parameters are estimated
#'   through marginalized maximum likelihood estimation.
#'   \item \code{\link{lhmm}} fits LHMM models. Parameters are estimated
#'   through marginalized maximum likelihood estimation.
#'   \item \code{\link{compute_theta}} compute MAP estimates of latent traits in LHMM.
#'   \item \code{\link{find_state_seq}} compute the most likely hidden state sequence.
#' }
#'
#' @section References:
#'
#' Tang, X. (2024) Latent Hidden Markov Models for Response Process Data. Psychometrika 89, 205-240. <https://doi.org/10.1007/s11336-023-09938-1>
#'
#' @docType _PACKAGE
#' @name proclhmm
#' @useDynLib proclhmm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
#>NULL
