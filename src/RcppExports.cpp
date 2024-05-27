// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_fb_matrix
void compute_fb_matrix(IntegerVector Y, int L, int N, int K, NumericVector P1, NumericMatrix P, NumericMatrix Q, NumericMatrix forward, NumericMatrix backward, NumericMatrix gamma, NumericVector scaling);
RcppExport SEXP _proclhmm_compute_fb_matrix(SEXP YSEXP, SEXP LSEXP, SEXP NSEXP, SEXP KSEXP, SEXP P1SEXP, SEXP PSEXP, SEXP QSEXP, SEXP forwardSEXP, SEXP backwardSEXP, SEXP gammaSEXP, SEXP scalingSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type forward(forwardSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type backward(backwardSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scaling(scalingSEXP);
    compute_fb_matrix(Y, L, N, K, P1, P, Q, forward, backward, gamma, scaling);
    return R_NilValue;
END_RCPP
}
// seq2llh
double seq2llh(IntegerVector Y, NumericVector P1, NumericMatrix P, NumericMatrix Q);
RcppExport SEXP _proclhmm_seq2llh(SEXP YSEXP, SEXP P1SEXP, SEXP PSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(seq2llh(Y, P1, P, Q));
    return rcpp_result_gen;
END_RCPP
}
// compute_llh
double compute_llh(List Y, int n, NumericVector P1, NumericMatrix P, NumericMatrix Q);
RcppExport SEXP _proclhmm_compute_llh(SEXP YSEXP, SEXP nSEXP, SEXP P1SEXP, SEXP PSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_llh(Y, n, P1, P, Q));
    return rcpp_result_gen;
END_RCPP
}
// seq2llh_gr
List seq2llh_gr(IntegerVector Y, NumericVector P1, NumericMatrix P, NumericMatrix Q);
RcppExport SEXP _proclhmm_seq2llh_gr(SEXP YSEXP, SEXP P1SEXP, SEXP PSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(seq2llh_gr(Y, P1, P, Q));
    return rcpp_result_gen;
END_RCPP
}
// compute_llh_gr
List compute_llh_gr(List Y, NumericVector P1, NumericMatrix P, NumericMatrix Q);
RcppExport SEXP _proclhmm_compute_llh_gr(SEXP YSEXP, SEXP P1SEXP, SEXP PSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_llh_gr(Y, P1, P, Q));
    return rcpp_result_gen;
END_RCPP
}
// compute_PQ_cpp
List compute_PQ_cpp(double theta, NumericMatrix para_a, NumericMatrix para_b, NumericMatrix para_alpha, NumericMatrix para_beta);
RcppExport SEXP _proclhmm_compute_PQ_cpp(SEXP thetaSEXP, SEXP para_aSEXP, SEXP para_bSEXP, SEXP para_alphaSEXP, SEXP para_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_a(para_aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_b(para_bSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_alpha(para_alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_beta(para_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_PQ_cpp(theta, para_a, para_b, para_alpha, para_beta));
    return rcpp_result_gen;
END_RCPP
}
// compute_P_cpp
NumericMatrix compute_P_cpp(double theta, NumericMatrix para_a, NumericMatrix para_b);
RcppExport SEXP _proclhmm_compute_P_cpp(SEXP thetaSEXP, SEXP para_aSEXP, SEXP para_bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_a(para_aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_b(para_bSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_P_cpp(theta, para_a, para_b));
    return rcpp_result_gen;
END_RCPP
}
// compute_Q_cpp
NumericMatrix compute_Q_cpp(double theta, NumericMatrix para_alpha, NumericMatrix para_beta);
RcppExport SEXP _proclhmm_compute_Q_cpp(SEXP thetaSEXP, SEXP para_alphaSEXP, SEXP para_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_alpha(para_alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_beta(para_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Q_cpp(theta, para_alpha, para_beta));
    return rcpp_result_gen;
END_RCPP
}
// compute_P1_cpp
NumericVector compute_P1_cpp(NumericVector para_P1);
RcppExport SEXP _proclhmm_compute_P1_cpp(SEXP para_P1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type para_P1(para_P1SEXP);
    rcpp_result_gen = Rcpp::wrap(compute_P1_cpp(para_P1));
    return rcpp_result_gen;
END_RCPP
}
// seq2llh_re_cpp
double seq2llh_re_cpp(IntegerVector Y, NumericMatrix para_a, NumericMatrix para_b, NumericMatrix para_alpha, NumericMatrix para_beta, NumericVector para_P1, NumericVector quad_nodes, NumericVector quad_weights);
RcppExport SEXP _proclhmm_seq2llh_re_cpp(SEXP YSEXP, SEXP para_aSEXP, SEXP para_bSEXP, SEXP para_alphaSEXP, SEXP para_betaSEXP, SEXP para_P1SEXP, SEXP quad_nodesSEXP, SEXP quad_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_a(para_aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_b(para_bSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_alpha(para_alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_beta(para_betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para_P1(para_P1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type quad_nodes(quad_nodesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type quad_weights(quad_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(seq2llh_re_cpp(Y, para_a, para_b, para_alpha, para_beta, para_P1, quad_nodes, quad_weights));
    return rcpp_result_gen;
END_RCPP
}
// compute_llh_rehmm
double compute_llh_rehmm(List Y, NumericMatrix para_a, NumericMatrix para_b, NumericMatrix para_alpha, NumericMatrix para_beta, NumericVector para_P1, NumericVector quad_nodes, NumericVector quad_weights);
RcppExport SEXP _proclhmm_compute_llh_rehmm(SEXP YSEXP, SEXP para_aSEXP, SEXP para_bSEXP, SEXP para_alphaSEXP, SEXP para_betaSEXP, SEXP para_P1SEXP, SEXP quad_nodesSEXP, SEXP quad_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_a(para_aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_b(para_bSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_alpha(para_alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_beta(para_betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para_P1(para_P1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type quad_nodes(quad_nodesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type quad_weights(quad_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_llh_rehmm(Y, para_a, para_b, para_alpha, para_beta, para_P1, quad_nodes, quad_weights));
    return rcpp_result_gen;
END_RCPP
}
// seq2llh_gr_rehmm
List seq2llh_gr_rehmm(IntegerVector Y, NumericMatrix para_a, NumericMatrix para_b, NumericMatrix para_alpha, NumericMatrix para_beta, NumericVector para_P1, NumericVector quad_nodes, NumericVector quad_weights);
RcppExport SEXP _proclhmm_seq2llh_gr_rehmm(SEXP YSEXP, SEXP para_aSEXP, SEXP para_bSEXP, SEXP para_alphaSEXP, SEXP para_betaSEXP, SEXP para_P1SEXP, SEXP quad_nodesSEXP, SEXP quad_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_a(para_aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_b(para_bSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_alpha(para_alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type para_beta(para_betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type para_P1(para_P1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type quad_nodes(quad_nodesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type quad_weights(quad_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(seq2llh_gr_rehmm(Y, para_a, para_b, para_alpha, para_beta, para_P1, quad_nodes, quad_weights));
    return rcpp_result_gen;
END_RCPP
}
// HMM_mle
void HMM_mle(List Y, int n, int N, int K, NumericVector P1, NumericMatrix P, NumericMatrix Q, double tot, int maxit);
RcppExport SEXP _proclhmm_HMM_mle(SEXP YSEXP, SEXP nSEXP, SEXP NSEXP, SEXP KSEXP, SEXP P1SEXP, SEXP PSEXP, SEXP QSEXP, SEXP totSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type tot(totSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    HMM_mle(Y, n, N, K, P1, P, Q, tot, maxit);
    return R_NilValue;
END_RCPP
}
// find_state_seq
IntegerVector find_state_seq(IntegerVector seq, NumericVector P1, NumericMatrix P, NumericMatrix Q);
RcppExport SEXP _proclhmm_find_state_seq(SEXP seqSEXP, SEXP P1SEXP, SEXP PSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(find_state_seq(seq, P1, P, Q));
    return rcpp_result_gen;
END_RCPP
}
// seq2cllh
double seq2cllh(IntegerVector Y, NumericVector P1, NumericMatrix P, NumericMatrix Q);
RcppExport SEXP _proclhmm_seq2cllh(SEXP YSEXP, SEXP P1SEXP, SEXP PSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(seq2cllh(Y, P1, P, Q));
    return rcpp_result_gen;
END_RCPP
}
// compute_cllh
double compute_cllh(List Y, int n, NumericVector P1, NumericMatrix P, NumericMatrix Q);
RcppExport SEXP _proclhmm_compute_cllh(SEXP YSEXP, SEXP nSEXP, SEXP P1SEXP, SEXP PSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_cllh(Y, n, P1, P, Q));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_proclhmm_compute_fb_matrix", (DL_FUNC) &_proclhmm_compute_fb_matrix, 11},
    {"_proclhmm_seq2llh", (DL_FUNC) &_proclhmm_seq2llh, 4},
    {"_proclhmm_compute_llh", (DL_FUNC) &_proclhmm_compute_llh, 5},
    {"_proclhmm_seq2llh_gr", (DL_FUNC) &_proclhmm_seq2llh_gr, 4},
    {"_proclhmm_compute_llh_gr", (DL_FUNC) &_proclhmm_compute_llh_gr, 4},
    {"_proclhmm_compute_PQ_cpp", (DL_FUNC) &_proclhmm_compute_PQ_cpp, 5},
    {"_proclhmm_compute_P_cpp", (DL_FUNC) &_proclhmm_compute_P_cpp, 3},
    {"_proclhmm_compute_Q_cpp", (DL_FUNC) &_proclhmm_compute_Q_cpp, 3},
    {"_proclhmm_compute_P1_cpp", (DL_FUNC) &_proclhmm_compute_P1_cpp, 1},
    {"_proclhmm_seq2llh_re_cpp", (DL_FUNC) &_proclhmm_seq2llh_re_cpp, 8},
    {"_proclhmm_compute_llh_rehmm", (DL_FUNC) &_proclhmm_compute_llh_rehmm, 8},
    {"_proclhmm_seq2llh_gr_rehmm", (DL_FUNC) &_proclhmm_seq2llh_gr_rehmm, 8},
    {"_proclhmm_HMM_mle", (DL_FUNC) &_proclhmm_HMM_mle, 9},
    {"_proclhmm_find_state_seq", (DL_FUNC) &_proclhmm_find_state_seq, 4},
    {"_proclhmm_seq2cllh", (DL_FUNC) &_proclhmm_seq2cllh, 4},
    {"_proclhmm_compute_cllh", (DL_FUNC) &_proclhmm_compute_cllh, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_proclhmm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}