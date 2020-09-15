#include <TMB.hpp>
#include <ktools.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);

  // data
  DATA_VECTOR(afs);
  DATA_VECTOR(afs_u);
  DATA_VECTOR(afs_l);
  DATA_IVECTOR(age);
  DATA_IVECTOR(event);
  DATA_IVECTOR(yob);

  DATA_IVECTOR(cc_id);
  DATA_IVECTOR(ccxyob_id);
  DATA_IVECTOR(ccxage_id);
  DATA_VECTOR(svw);

  // priors
  DATA_VECTOR(sd_beta);
  DATA_VECTOR(sd_yob);
  DATA_VECTOR(sd_age);

  DATA_VECTOR(sd_cc); // spatial
  DATA_VECTOR(sd_ccxyob); // interaction
  DATA_VECTOR(sd_ccxage); // interaction

  DATA_VECTOR(palpha);
  DATA_VECTOR(p_a);

  DATA_MATRIX(R_age);
  DATA_MATRIX(R_yob);

  DATA_MATRIX(R_cc);
  DATA_INTEGER(R_cc_rank);

  DATA_MATRIX(R_ccxyob); // R_cc X R_yob
  DATA_INTEGER(R_ccxyob_rank);
  DATA_MATRIX(R_ccxage); // R_cc X R_age
  DATA_INTEGER(R_ccxage_rank);

  // Data model - log-logistic parameters
  PARAMETER(intercept);

  Type prior = 0.0; // this kills parallel
  prior -= dnorm(intercept, sd_beta(0), sd_beta(1), true);

  PARAMETER_VECTOR(log_alpha_vec); // do Q will be nicer
  prior -= dlgamma(log_alpha_vec, palpha(0), palpha(1), true).sum();
  vector<Type> alpha_vec = exp(log_alpha_vec);

  PARAMETER_VECTOR(log_a_vec);
  prior -= dlgamma(log_a_vec, p_a(0), p_a(1), true).sum();
  vector<Type> a_vec = exp(log_a_vec);

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
  prior += density::GMRF(ktools::prepare_Q(R_cc, cc_e))(cc_vec);
  prior += (R_cc_rank - cc_vec.size()) * log(sqrt(2*M_PI)); // ktools::GMRF would be nice
  
  // countries x yob interaction
  PARAMETER_VECTOR  (ccxyob);
  PARAMETER         (log_ccxyob_e);
  Type ccxyob_e = exp(log_ccxyob_e);
  prior -= ktools::pc_prec(ccxyob_e, sd_ccxyob(0), sd_ccxyob(1));
  for (int j = 0; j < cc_vec.size(); ++j) {
    vector<Type> v_j = ccxyob.segment(j * yob_rw2.size(), yob_rw2.size());
    prior -= ktools::soft_zero_sum(v_j);
  }
  prior += density::GMRF(ktools::prepare_Q(R_ccxyob, ccxyob_e))(ccxyob);
  prior += (R_ccxyob_rank - ccxyob.size()) * log(sqrt(2*M_PI)); // ktools::GMRF would be nice

  // countries x age interaction
  PARAMETER_VECTOR  (ccxage);
  PARAMETER         (log_ccxage_e);
  Type ccxage_e = exp(log_ccxage_e);
  prior -= ktools::pc_prec(ccxage_e, sd_ccxage(0), sd_ccxage(1));
  for (int j = 0; j < cc_vec.size(); ++j) {
    vector<Type> v_j = ccxage.segment(j * age_rw2.size(), age_rw2.size());
    prior -= ktools::soft_zero_sum(v_j);
  }
  prior += density::GMRF(ktools::prepare_Q(R_ccxage, ccxage_e))(ccxage);
  prior += (R_ccxage_rank - ccxage.size()) * log(sqrt(2*M_PI)); // ktools::GMRF would be nice

  // Data likelihood
  for (int i = 0; i < afs.size(); i++) {
    Type eta = intercept + 
      yob_rw2(yob(i)) + 
      age_rw2(age(i)) + 
      cc_vec(cc_id(i)) + 
      ccxyob(ccxyob_id(i)) +
      ccxage(ccxage_id(i));
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
  REPORT(prior);;
  return dll;
}