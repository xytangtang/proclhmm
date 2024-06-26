% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{compute_PQ_lhmm}
\alias{compute_PQ_lhmm}
\title{Compute LHMM probabilities from parameters}
\usage{
compute_PQ_lhmm(theta, para_a, para_b, para_alpha, para_beta)
}
\arguments{
\item{theta}{latent trait}

\item{para_a}{\code{K} by \code{K-1} matrix. discrimination parameters of state transition probability matrix}

\item{para_b}{\code{K} by \code{K-1} matrix. location parameters of state transition probability matrix}

\item{para_alpha}{\code{K} by \code{N-1} matrix. discrimination parameters of state-action (emission) probability matrix}

\item{para_beta}{\code{K} by \code{N-1} matrix. location parameters of state-action (emission) probability matrix}
}
\value{
A list of two elements
\tabular{ll}{
{\code{P}} \tab {\code{K} by \code{K} state-transition probability matrix} \cr
\tab \cr
{\code{Q}} \tab {\code{K} by \code{N} state-action probability matrix} \cr
\tab \cr
}
}
\description{
Compute state-transition and state-action (emission) probability matrices
from LHMM parameters
}
\examples{
paras <- sim_lhmm_paras(5, 2)
prob_paras <- compute_PQ_lhmm(1.5, paras$para_a, paras$para_b, paras$para_alpha, paras$para_beta)

}
\seealso{
\code{\link{compute_P1_lhmm}} for initial state probabilities of
LHMM, \code{\link{compute_paras_hmm}} for computing probabilities in HMM.
}
