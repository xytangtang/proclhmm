#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// function for calculating forward and backward matrix for a single sequence
//[[Rcpp::export]]
void compute_fb_matrix(IntegerVector Y, int L, int N, int K,
                       NumericVector P1, NumericMatrix P, NumericMatrix Q,
                       NumericMatrix forward, NumericMatrix backward, NumericMatrix gamma, NumericVector scaling) {
  int k = 0, j = 0, l = 0;
  // forward prob
  for (k = 0; k < K; k++) {
    forward(k, 0) = P1[k] * Q(k, Y[0]);
    scaling[0] += forward(k, 0);
  }
  for (k = 0; k < K; k++) {
    forward(k, 0) /= scaling[0];
  }
  for (j = 1; j < L; j++) {
    for (k = 0; k < K; k++) {
      for (l = 0; l < K; l++) {
        forward(k, j) += forward(l, j-1) * P(l, k) * Q(k, Y[j]);
      }
      scaling[j] += forward(k, j);
    }
    for (k = 0; k < K; k++) {
      forward(k, j) /= scaling[j];
    }
  }


  // backward prob
  for (k = 0; k < K; k++) {
    backward(k, L - 1) = 1.0;
  }
  for (j = 1; j < L; j++) {
    for (k = 0; k < K; k++) {
      for (l = 0; l < K; l++) {
        backward(k, L - j - 1) += backward(l, L - j) * P(k, l) * Q(l, Y[L - j]);
      }
    }
    for (k = 0; k < K; k++) {
      backward(k, L - j - 1) /= scaling[L - j];
    }
  }

  for (int j = 0; j < L; j++) {
    for (int k = 0; k < K; k++) {
      gamma(k, j) = forward(k, j) * backward(k, j);
    }
  }

}

// [[Rcpp::export]]
double seq2llh(IntegerVector Y, NumericVector P1, NumericMatrix P, NumericMatrix Q) {

  double llh = 0.0;
  int N = Q.ncol();
  int K = P1.size();
  int L = Y.size();

  NumericMatrix forward(K, L);
  NumericMatrix backward(K, L);
  NumericMatrix gamma(K, L);
  NumericVector scaling(L);

  compute_fb_matrix(Y, L, N, K, P1, P, Q, forward, backward, gamma, scaling);

  for (int j = 0; j < L; j++) {
    llh += log(scaling[j]);
  }

  return llh;

}

//[[Rcpp::export]]
double compute_llh(List Y, int n, NumericVector P1, NumericMatrix P, NumericMatrix Q) {

  double llh = 0;
  IntegerVector seq;
  for (int i = 0; i < n; i++) {
    seq = Y[i];
    llh += seq2llh(seq, P1, P, Q);
  }

  return llh;
}



// [[Rcpp::export]]
List seq2llh_gr(IntegerVector Y, NumericVector P1, NumericMatrix P, NumericMatrix Q) {

  int N = Q.ncol();
  int K = P1.size();
  int L = Y.size();
  int k, l, t;

  NumericMatrix forward(K, L);
  NumericMatrix backward(K, L);
  NumericMatrix gamma(K, L);
  NumericVector scaling(L);

  compute_fb_matrix(Y, L, N, K, P1, P, Q, forward, backward, gamma, scaling);

  NumericVector P1_gr(K);
  NumericMatrix P_gr(K, K);
  NumericMatrix Q_gr(K, N);

  for (k = 0; k < K; k++) {
    P1_gr[k] = Q(k, Y[0]) * backward(k, 0) / scaling[0];
  }
  for (k = 0; k < K; k++) {
    for (l = 0; l < K; l++) {
      for (t = 0; t < L-1; t++) {
        P_gr(k, l) += (forward(k, t) * backward(l, t+1) * Q(l, Y[t+1]) / scaling(t+1));
      }
    }
  }
  for (k = 0; k < K; k++) {
    for (t = 0; t < L; t++) {
      Q_gr(k, Y[t]) += (forward(k, t) * backward(k, t) / Q(k, Y[t]));
    }
  }

  List out;
  out["P1_gr"] = P1_gr;
  out["P_gr"] = P_gr;
  out["Q_gr"] = Q_gr;

  return out;

}



// [[Rcpp::export]]
List compute_llh_gr(List Y, NumericVector P1, NumericMatrix P, NumericMatrix Q) {

  int N = Q.ncol();
  int K = P1.size();
  int n = Y.size();
  int i, k, l, j;
  IntegerVector seq;
  NumericVector P1_gr(K);
  NumericMatrix P_gr(K, K);
  NumericMatrix Q_gr(K, N);

  for (i = 0; i < n; i++) {
    seq = Y[i];
    List seq_gr = seq2llh_gr(seq, P1, P, Q);

    NumericVector seq_P1_gr = seq_gr["P1_gr"];
    NumericMatrix seq_P_gr = seq_gr["P_gr"];
    NumericMatrix seq_Q_gr = seq_gr["Q_gr"];

    for (k = 0; k < K; k++) {
      P1_gr[k] += seq_P1_gr[k];
      for (l = 0; l < K; l++) {
        P_gr(k, l) += seq_P_gr(k, l);
      }
      for (j = 0; j < N; j++) {
        Q_gr(k, j) += seq_Q_gr(k, j);
      }
    }
  }

  List out;
  out["P1_gr"] = P1_gr;
  out["P_gr"] = P_gr;
  out["Q_gr"] = Q_gr;

  return out;
}


// [[Rcpp::export]]
List compute_PQ_cpp(double theta,
                    NumericMatrix para_a,
                    NumericMatrix para_b,
                    NumericMatrix para_alpha,
                    NumericMatrix para_beta) {
  int K = para_a.nrow();
  int N = para_alpha.ncol() + 1;

  NumericMatrix P(K, K);
  NumericMatrix Q(K, N);
  NumericVector P_rowsum(K);
  NumericVector Q_rowsum(N);

  for (int i = 0; i < K; i++) {
    P(i, 0) = 1.0;
    P_rowsum[i] = P(i, 0);
    for (int j = 1; j < K; j++) {
      P(i,j) = exp(para_a(i, j-1) * theta + para_b(i, j-1));
      P_rowsum[i] += P(i, j);
    }

    for (int j = 0; j < K; j++) {
      P(i, j) = P(i, j) / P_rowsum[i];
    }
  }

  for (int i = 0; i < K; i++) {
    Q(i, 0) = 1.0;
    Q_rowsum[i] = Q(i, 0);
    for (int j = 1; j < N; j++) {
      Q(i,j) = exp(para_alpha(i, j-1) * theta + para_beta(i, j-1));
      Q_rowsum[i] += Q(i, j);
    }

    for (int j = 0; j < N; j++) {
      Q(i, j) = Q(i, j) / Q_rowsum[i];
    }
  }

  List out;
  out["P"] = P;
  out["Q"] = Q;

  return out;

}


// [[Rcpp::export]]
NumericMatrix compute_P_cpp(double theta,
                            NumericMatrix para_a,
                            NumericMatrix para_b) {
  int K = para_a.nrow();

  NumericMatrix P(K, K);
  NumericVector P_rowsum(K);

  for (int i = 0; i < K; i++) {
    P(i, 0) = 1.0;
    P_rowsum[i] = P(i, 0);
    for (int j = 1; j < K; j++) {
      P(i,j) = exp(para_a(i, j-1) * theta + para_b(i, j-1));
      P_rowsum[i] += P(i, j);
    }

    for (int j = 0; j < K; j++) {
      P(i, j) = P(i, j) / P_rowsum[i];
    }
  }

  return P;

}

// [[Rcpp::export]]
NumericMatrix compute_Q_cpp(double theta,
                            NumericMatrix para_alpha,
                            NumericMatrix para_beta) {
  int K = para_alpha.nrow();
  int N = para_alpha.ncol() + 1;

  NumericMatrix Q(K, N);
  NumericVector Q_rowsum(N);

  for (int i = 0; i < K; i++) {
    Q(i, 0) = 1.0;
    Q_rowsum[i] = Q(i, 0);
    for (int j = 1; j < N; j++) {
      Q(i,j) = exp(para_alpha(i, j-1) * theta + para_beta(i, j-1));
      Q_rowsum[i] += Q(i, j);
    }

    for (int j = 0; j < N; j++) {
      Q(i, j) = Q(i, j) / Q_rowsum[i];
    }
  }

  return Q;

}


// [[Rcpp::export]]
NumericVector compute_P1_cpp(NumericVector para_P1) {

  int K = para_P1.size() + 1;
  NumericVector P1(K);
  double P1_sum = 0.0;
  P1[0] = 1.0;
  P1_sum = P1[0];
  for (int i = 1; i < K; i++) {
    P1[i] = exp(para_P1[i-1]);
    P1_sum += P1[i];
  }
  for (int i = 0; i < K; i++) {
    P1[i] = P1[i] / P1_sum;
  }
  return P1;
}

// [[Rcpp::export]]
double seq2llh_re_cpp(IntegerVector Y,
                      NumericMatrix para_a,
                      NumericMatrix para_b,
                      NumericMatrix para_alpha,
                      NumericMatrix para_beta,
                      NumericVector para_P1,
                      NumericVector quad_nodes,
                      NumericVector quad_weights) {

  int n_pts = quad_nodes.size();
  NumericVector llh_list(n_pts);

  NumericVector P1 = compute_P1_cpp(para_P1);

  for (int index_pt = 0; index_pt < n_pts; index_pt++) {
    List paras = compute_PQ_cpp(sqrt(2.0)*quad_nodes[index_pt], para_a, para_b, para_alpha, para_beta);
    llh_list[index_pt] = seq2llh(Y, P1, paras["P"], paras["Q"]) + log(quad_weights[index_pt]);
  }

  double max_llh = max(llh_list);

  return max_llh + log(sum(exp(llh_list - max_llh)));

}

// [[Rcpp::export]]
double compute_llh_rehmm(List Y,
                         NumericMatrix para_a,
                         NumericMatrix para_b,
                         NumericMatrix para_alpha,
                         NumericMatrix para_beta,
                         NumericVector para_P1,
                         NumericVector quad_nodes,
                         NumericVector quad_weights) {
  int N = para_a.ncol() + 1;
  int K = para_a.nrow();
  int n_pts = quad_nodes.size();
  int n = Y.size();
  double out = 0.0;

  for (int i = 0; i < n; i++) {

    IntegerVector seq = Y[i];
    NumericVector llh_list(n_pts);

    NumericVector P1(K);
    NumericMatrix P(K, K);
    NumericMatrix Q(K, N);

    P1 = compute_P1_cpp(para_P1);

    for (int index_pt = 0; index_pt < n_pts; index_pt++) {
      List paras = compute_PQ_cpp(sqrt(2.0)*quad_nodes[index_pt], para_a, para_b, para_alpha, para_beta);
      llh_list[index_pt] = seq2llh(seq, P1, paras["P"], paras["Q"]) + log(quad_weights[index_pt]);
    }

    double max_llh = max(llh_list);

    out += (max_llh + log(sum(exp(llh_list - max_llh))));
  }

  return out;

}


// [[Rcpp::export]]
List seq2llh_gr_rehmm(IntegerVector Y,
                      NumericMatrix para_a,
                      NumericMatrix para_b,
                      NumericMatrix para_alpha,
                      NumericMatrix para_beta,
                      NumericVector para_P1,
                      NumericVector quad_nodes,
                      NumericVector quad_weights) {

  int N = para_beta.ncol() + 1;
  int K = para_a.nrow();
  int n_pts = quad_nodes.size();

  NumericVector llh_list(n_pts);
  NumericVector llh_c_list(n_pts);

  NumericVector P1 = compute_P1_cpp(para_P1);

  for (int index_pt = 0; index_pt < n_pts; index_pt++) {
    List paras = compute_PQ_cpp(sqrt(2.0)*quad_nodes[index_pt], para_a, para_b, para_alpha, para_beta);

    llh_c_list[index_pt] = seq2llh(Y, P1, paras["P"], paras["Q"]);
    llh_list[index_pt] = llh_c_list[index_pt] + log(quad_weights[index_pt]);
  }
  double llh = max(llh_list);
  llh = llh + log(sum(exp(llh_list - llh)));

  NumericVector para_P1_gr(K-1);
  NumericMatrix para_a_gr(K, K-1);
  NumericMatrix para_b_gr(K, K-1);
  NumericMatrix para_alpha_gr(K, N-1);
  NumericMatrix para_beta_gr(K, N-1);

  for (int index_pt = 0; index_pt < n_pts; index_pt++) {
    List paras = compute_PQ_cpp(sqrt(2.0)*quad_nodes[index_pt], para_a, para_b, para_alpha, para_beta);
    NumericMatrix P = paras["P"];
    NumericMatrix Q = paras["Q"];

    List seq_gr = seq2llh_gr(Y, P1, P, Q);
    NumericVector seq_P1_gr = seq_gr["P1_gr"];
    NumericMatrix seq_P_gr = seq_gr["P_gr"];
    NumericMatrix seq_Q_gr = seq_gr["Q_gr"];

    double sum_P1 = 0.0;
    for (int k = 0; k < K; k++) sum_P1 += (seq_P1_gr[k] * P1[k]);
    for (int k = 0; k < K - 1; k++) {
      para_P1_gr[k] += (exp(llh_c_list[index_pt] - llh) * P1[k+1] * (seq_P1_gr[k+1] - sum_P1) * quad_weights[index_pt]);
    }

    for (int k = 0; k < K; k++) {
      double sum_P = 0.0;
      for (int l = 0; l < K; l++) sum_P += (seq_P_gr(k,l) * P(k,l));
      for (int l = 0; l < K-1; l++) {
        para_a_gr(k,l) += (exp(llh_c_list[index_pt] - llh) * sqrt(2.0)*quad_nodes[index_pt]*P(k,l+1) * (seq_P_gr(k,l+1) - sum_P) * quad_weights[index_pt]);
        para_b_gr(k,l) += (exp(llh_c_list[index_pt] - llh) * P(k,l+1) * (seq_P_gr(k,l+1) - sum_P) * quad_weights[index_pt]);
      }

      double sum_Q = 0.0;
      for (int j = 0; j < N; j++) sum_Q += (seq_Q_gr(k,j) * Q(k, j));
      for (int j = 0; j < N-1; j++) {
        para_alpha_gr(k,j) += (exp(llh_c_list[index_pt] - llh) * sqrt(2.0)*quad_nodes[index_pt]*Q(k,j+1) * (seq_Q_gr(k,j+1) - sum_Q) * quad_weights[index_pt]);
        para_beta_gr(k,j) += (exp(llh_c_list[index_pt] - llh) * Q(k,j+1) * (seq_Q_gr(k,j+1) - sum_Q) * quad_weights[index_pt]);
      }
    }
  }

  List out;
  out["para_a_gr"] = para_a_gr;
  out["para_b_gr"] = para_b_gr;
  out["para_alpha_gr"] = para_alpha_gr;
  out["para_beta_gr"] = para_beta_gr;
  out["para_P1_gr"] = para_P1_gr;

  return out;

}

// forward-backward algorithm for hidden markov model
//[[Rcpp::export]]
void HMM_mle(List Y, int n, int N, int K,
             NumericVector P1, NumericMatrix P, NumericMatrix Q,
             double tot, int maxit) {

  int i, j, k, l;
  int count = 0;

  double llh = 0.0;
  double llh_new = compute_llh(Y, n, P1, P, Q); // calculate loglikelihood likelihood

  while ((fabs(llh - llh_new) > tot) & (count < maxit)) {
    double temp = 0.0;
    count += 1;
    if (count % 100 == 0) Rcout << "Iteration: " << count << "\n";
    NumericVector P1_new(K);
    NumericMatrix P_new(K, K);
    NumericMatrix P_den(K);
    NumericMatrix Q_new(K, N);
    NumericMatrix Q_den(K);

    for (i=0; i < n; i++) {
      IntegerVector seq = Y[i];
      int L = seq.size();
      NumericMatrix forward_matrix(K, L);
      NumericMatrix backward_matrix(K, L);
      NumericMatrix gamma_matrix(K, L);
      NumericVector scaling(L);
      compute_fb_matrix(seq, L, N, K, P1, P, Q, forward_matrix, backward_matrix, gamma_matrix, scaling);

      // update num and den of probabilities

      for (k = 0; k < K; k++) {
        P1_new[k] += gamma_matrix(k, 0);

        for (j = 0; j < L - 1; j++) {

          for (int l = 0; l < K; l++){
            temp = forward_matrix(k, j) * P(k, l) * Q(l, seq[j + 1]) * backward_matrix(l, j + 1) / scaling[j + 1];
            P_new(k, l) += temp;
            P_den[k] += temp;
          }
          Q_new(k, seq[j]) += gamma_matrix(k, j);
          Q_den[k] += gamma_matrix(k, j);
        }

        Q_new(k, seq[L - 1]) += gamma_matrix(k, L - 1);
        Q_den[k] += gamma_matrix(k, L - 1);
      }

    }

    // update probabilities
    for (k = 0; k < K; k++) {
      P1[k] = P1_new[k] / n;

      for (l = 0; l < K; l++) {
        P(k, l) = P_new(k, l) / P_den[k];
      }

      for (j = 0; j < N; j++) {
        Q(k, j) = Q_new(k, j) / Q_den[k];
      }
    }

    // calculate loglikelihood for checking convergence
    llh = llh_new;
    llh_new = compute_llh(Y, n, P1, P, Q);
    if (count % 100 == 0) Rcout << "llh diff: " << fabs(llh - llh_new) << "\n";

  }
}

//' Viterbi algorithm for HMM
//'
//' Find the most likely hidden state sequence of an observed sequence under HMM
//'
//' @param seq An action sequence coded in integers
//' @param P1 initial state probability vector of length \code{K}
//' @param P \code{K} by \code{K} state transition probability matrix
//' @param Q \code{K} by \code{N} state-action (emission) probability matrix
//'
//' @return a hidden state sequence coded in integers
//'
//' @export
//[[Rcpp::export]]
IntegerVector find_state_seq(IntegerVector seq, NumericVector P1, NumericMatrix P, NumericMatrix Q) {
  int L = seq.size();
  int K = P1.size();
  int k = 0, j = 0, l = 0;
  IntegerVector best_state_seq(L);

  NumericMatrix v_mat(K, L); // keep track of best prob
  NumericVector scaling(L); // keep track of scaling factor
  IntegerMatrix p_mat(K, L); // keep track of best state

  NumericVector v_vec(K);

  // initialization
  for (k = 0; k < K; k++) {
    p_mat(k, 0) = 0;
    v_mat(k, 0) = P1[k] * Q(k, seq[0]);
    scaling[0] += v_mat(k, 0);
  }
  for (k = 0; k < K; k++) {
    v_mat(k, 0) /= scaling[0];
  }

  if (L > 1) {
    for (j = 1; j < L; j++) {
      for (k = 0; k < K; k++) {
        double max_value = 0.0;
        int max_id = 0;
        for (l = 0; l < K; l++) {
          v_vec[l] = v_mat(l, j - 1) * P(l, k) * Q(k, seq[j]);
          if (v_vec[l] > max_value) {
            max_value = v_vec[l];
            max_id = l;
          }
        }
        v_mat(k, j) = max_value;
        p_mat(k, j) = max_id;
        scaling[j] += v_mat(k, j);
      }
      for (k = 0; k < K; k++) {
        v_mat(k, j) /= scaling[j];
      }
    }
  }

  double max_value = 0.0;
  int max_id = 0;
  for (k = 0; k < K; k++) {
    if (v_mat(k, L - 1) > max_value) {
      max_value = v_mat(k, L - 1);
      max_id = k;
    }
  }

  best_state_seq[L - 1] = max_id;
  if (L > 1) {
    for (j = L - 1; j > 0; j--) {
      best_state_seq[j - 1] = p_mat(best_state_seq[j], j);
    }
  }

  return best_state_seq;
}

// [[Rcpp::export]]
double seq2cllh(IntegerVector Y, NumericVector P1, NumericMatrix P, NumericMatrix Q) {

  int L = Y.size();
  IntegerVector Y_state(L);
  Y_state = find_state_seq(Y, P1, P, Q);

  double cllh = 0.0;

  int state_curr;
  int state_prev;
  state_curr = Y_state[0];
  cllh = log(P1[state_curr]) + log(Q(state_curr, Y[0]));

  for (int i = 1; i < L; i++) {
    state_prev = state_curr;
    state_curr = Y_state[i];

    cllh += (log(P(state_prev, state_curr)) + log(Q(state_curr, Y[i])));
  }

  return cllh;

}

//[[Rcpp::export]]
double compute_cllh(List Y, int n, NumericVector P1, NumericMatrix P, NumericMatrix Q) {

  double cllh = 0.0;
  IntegerVector seq;
  for (int i = 0; i < n; i++) {
    seq = Y[i];
    cllh += seq2cllh(seq, P1, P, Q);
  }

  return cllh;
}

