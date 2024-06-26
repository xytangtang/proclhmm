% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\name{sim_hmm_paras}
\alias{sim_hmm_paras}
\title{generate HMM parameters}
\usage{
sim_hmm_paras(N, K, return_prob = TRUE)
}
\arguments{
\item{N}{number of distinct actions}

\item{K}{number of hidden states}

\item{return_prob}{logical. indicates to return parameters in probability scale (\code{TRUE}, default) or logit scale.}
}
\value{
a list of three elements.
If \code{return_prob = TRUE}, the element names are \code{P1}, \code{P}, and \code{Q}.
If \code{return_prob = FALSE}, the element names are \code{para_P1}, \code{para_P}, and \code{oara_Q}.
}
\description{
\code{sim_hmm_paras} generates logit scale parameters of HMM with \code{K} hidden states and
\code{N} distinct actions from Uniform(-0.5, 0.5).
}
\examples{
# generate probability parameters
set.seed(12345)
paras1 <- sim_hmm_paras(5, 2)
names(paras1)

# generate parameters in the logit scale
set.seed(12345)
paras2 <- sim_hmm_paras(5, 2, return_prob = FALSE)
names(paras2)

paras1$P1
paras2$para_P1

# logit scale parameters can be transformed to probability parameters
all.equal(compute_paras_hmm(paras2$para_P, paras2$para_Q, paras2$para_P1), paras1)

}
