#include <R.h>
#include <TMB.hpp>
#include <omp.h>

// config.trace.parallel     = true;
// config.trace.optimize     = true;
// config.trace.atomic       = true;
// config.optimize.instantly = false; // Always optimize just after tape creation?
// config.optimize.parallel  = false; // Allow optimize in parallel (memory consuming)?
// config.tape.parallel      = false; // Enable parallel tape creation?

#include <ktools.hpp>

// prepare Q
template <class Type>
Eigen::SparseMatrix<Type> prepare_Q(matrix<Type> R, Type tau) {
  R = tau * R.array();
  R.diagonal().array() += 1e-8;
  return tmbutils::asSparseMatrix(R);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);
  // #ifdef _OPENMP
  //   int threadnum = omp_get_thread_num();
  //   Rprintf("OpenMP supported, thread #%d.\n", threadnum);
  // #endif
  // Model data
  DATA_VECTOR(afs);
  DATA_VECTOR(afs_u);
  DATA_VECTOR(afs_l);
  DATA_IVECTOR(age);
  DATA_IVECTOR(event);
  DATA_IVECTOR(yob);

  DATA_IVECTOR(cc_id);
  DATA_IVECTOR(cc_yob_id);
  DATA_VECTOR(svw);

  // priors
  DATA_VECTOR(sd_beta);
  DATA_VECTOR(sd_yob);
  DATA_VECTOR(sd_age);

  DATA_VECTOR(sd_cc); // spatial
  DATA_VECTOR(sd_cc_yob); // interaction

  DATA_VECTOR(palpha);
  // DATA_VECTOR(alpha_prec);

  DATA_VECTOR(p_a);
  // DATA_VECTOR(a_prec);

  DATA_MATRIX(R_age);
  DATA_MATRIX(R_yob);

  DATA_MATRIX(R_cc);
  DATA_INTEGER(R_cc_rank);

  DATA_MATRIX(R_cc_yob); // R_cc X R_yob
  DATA_INTEGER(R_cc_yob_rank);

  // Data model - log-logistic parameters
  PARAMETER(intercept);

  Type prior = 0.0; // this kills parallel
  prior -= dnorm(intercept, sd_beta(0), sd_beta(1), true);

  // PARAMETER(log_alpha);
  // prior -= dlgamma(log_alpha, palpha(0), palpha(1), true);

  PARAMETER_VECTOR(log_alpha_vec); // do Q will be nicer
  prior -= dlgamma(log_alpha_vec, palpha(0), palpha(1), true).sum();
  vector<Type> alpha_vec = exp(log_alpha_vec);
  // PARAMETER(log_alpha_vec_prec);
  // Type alpha_vec_sd = pow(exp(log_alpha_vec_prec), -0.5);
  // prior -= ktools::pc_prec(exp(log_alpha_vec_prec), alpha_prec(0), alpha_prec(1));
  // prior -= dnorm(alpha_vec, Type(0), alpha_vec_sd).sum();
  // prior -= ktools::soft_zero_sum(alpha_vec);
  // for (int i = 0; i < alpha_vec.size(); ++i)
    // alpha_vec(i) += exp(log_alpha);

  // PARAMETER(log_a);
  // prior -= dlgamma(log_a, p_a(0), p_a(1), true);

  PARAMETER_VECTOR(log_a_vec);
  prior -= dlgamma(log_a_vec, p_a(0), p_a(1), true).sum();
  vector<Type> a_vec = exp(log_a_vec);
  // PARAMETER(log_a_vec_prec);
  // prior -= ktools::pc_prec(exp(log_a_vec_prec), a_prec(0), a_prec(1));
  // Type a_vec_sd = pow(exp(log_a_vec_prec), -0.5);
  // prior -= dnorm(a_vec, Type(0), a_vec_sd).sum();
  // prior -= ktools::soft_zero_sum(a_vec);
  // for (int i = 0; i < a_vec.size(); ++i)
  //   a_vec(i) += exp(log_a);

  // yob rw2
  PARAMETER_VECTOR (yob_rw2);
  PARAMETER        (log_yob_rw2_e);
  Type yob_rw2_e = exp(log_yob_rw2_e);
  prior -= ktools::pc_prec(yob_rw2_e, sd_yob(0), sd_yob(1));
  prior += ktools::rw(yob_rw2, R_yob, yob_rw2_e);

  // age rw2
  PARAMETER_VECTOR (age_rw2);
  PARAMETER        (log_age_rw2_e);
  Type age_rw2_e = exp(log_age_rw2_e);
  prior -= ktools::pc_prec(age_rw2_e, sd_age(0), sd_age(1));
  prior += ktools::rw(age_rw2, R_age, age_rw2_e);

  // countries spatial
  PARAMETER_VECTOR  (cc_vec);
  PARAMETER         (log_cc_e);
  Type cc_e = exp(log_cc_e);
  prior -= ktools::pc_prec(cc_e, sd_cc(0), sd_cc(1));
  prior -= ktools::soft_zero_sum(cc_vec);
  prior += density::GMRF(prepare_Q(R_cc, cc_e))(cc_vec);
  prior += (R_cc_rank - cc_vec.size()) * log(sqrt(2*M_PI)); // ktools::GMRF would be nice
  
  // countries x yob interaction
  PARAMETER_VECTOR  (cc_yob);
  PARAMETER         (log_cc_yob_e);
  Type cc_yob_e = exp(log_cc_yob_e);
  prior -= ktools::pc_prec(cc_yob_e, sd_cc_yob(0), sd_cc_yob(1));
  for (int j = 0; j < cc_vec.size(); ++j) {
    vector<Type> v_j = cc_yob.segment(j * yob_rw2.size(), yob_rw2.size());
    prior -= ktools::soft_zero_sum(v_j);
  }
  prior += density::GMRF(prepare_Q(R_cc_yob, cc_yob_e))(cc_yob);
  prior += (R_cc_yob_rank - cc_yob.size()) * log(sqrt(2*M_PI)); // ktools::GMRF would be nice

  // Data likelihood
  // vector<Type> pwdens(afs.size()); // for LOO, WAIC
  for (int i = 0; i < afs.size(); i++) {
    Type eta = intercept + yob_rw2(yob(i)) + age_rw2(age(i)) + cc_vec(cc_id(i)) + cc_yob(cc_yob_id(i));
    Type lambda = exp(eta);
    if (event(i)) {
      dll -= log(
        svw(i) * (
          ktools::St_llogisI(afs_l(i), alpha_vec(cc_id(i)), lambda, a_vec(cc_id(i))) -
          ktools::St_llogisI(afs_u(i), alpha_vec(cc_id(i)), lambda, a_vec(cc_id(i)))
        ));
    } else {
      dll -= log(svw(i) * ktools::St_llogisI(afs(i), alpha_vec(cc_id(i)), lambda, a_vec(cc_id(i))));
    }
  }
  dll += prior;
  // dll += prior - sum(log(pwdens));
  REPORT(prior);;
  return dll;
}