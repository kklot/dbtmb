#include <TMB.hpp>
#include "ktools.hpp"

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
  
  DATA_SCALAR(age_order);
  DATA_SCALAR(yob_order);

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

  Type prior = 0.0;
  prior -= dnorm(intercept, sd_beta(0), sd_beta(1), true);

  // Shape
  PARAMETER_VECTOR(log_alpha_vec); 
  prior -= dnorm(log_alpha_vec, palpha(0), palpha(1), true).sum();
  vector<Type> alpha_vec = exp(log_alpha_vec);

  // Skewness *
  PARAMETER_VECTOR(a_vec_star);
  prior -= dnorm(a_vec_star, p_a(0), p_a(1), true).sum();
  // Skewness real
  vector<Type> a_vec = exp(a_vec_star - 1.1*log_alpha_vec);

  // yob rw2
  PARAMETER_VECTOR (yob_rw2);
  PARAMETER        (log_yob_rw2_e);
  Type yob_rw2_e = exp(log_yob_rw2_e);
  prior -= ktools::pc_prec(yob_rw2_e, sd_yob(0), sd_yob(1));
  prior += ktools::rw(yob_rw2, R_yob, yob_rw2_e, yob_order);

  // age rw2
  PARAMETER_VECTOR (age_rw2);
  PARAMETER        (log_age_rw2_e);
  Type age_rw2_e = exp(log_age_rw2_e);
  prior -= ktools::pc_prec(age_rw2_e, sd_age(0), sd_age(1));
  prior += ktools::rw(age_rw2, R_age, age_rw2_e, age_order);

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
  prior -= ktools::constraint2D(ccxyob.data(), yob_rw2.size(), cc_vec.size(), true, true, false, false);
  prior += density::GMRF(ktools::prepare_Q(R_ccxyob, ccxyob_e))(ccxyob);
  prior += (R_ccxyob_rank - ccxyob.size()) * log(sqrt(2*M_PI)); // ktools::GMRF would be nice

  // countries x age interaction
  PARAMETER_VECTOR  (ccxage);
  PARAMETER         (log_ccxage_e);
  Type ccxage_e = exp(log_ccxage_e);
  prior -= ktools::pc_prec(ccxage_e, sd_ccxage(0), sd_ccxage(1));
  prior -= ktools::constraint2D(ccxage.data(), age_rw2.size(), cc_vec.size());
  prior += density::GMRF(ktools::prepare_Q(R_ccxage, ccxage_e))(ccxage);
  prior += (R_ccxage_rank - ccxage.size()) * log(sqrt(2*M_PI)); // ktools::GMRF would be nice

  // Data likelihood
  for (int i = 0; i < afs.size(); i++) {
    Type eta = intercept + yob_rw2(yob(i)) + age_rw2(age(i)) + cc_vec(cc_id(i)) + ccxyob(ccxyob_id(i)) + ccxage(ccxage_id(i));
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
  // Reporting
  int nC = cc_vec.size(), nT = yob_rw2.size(), nA = age_rw2.size();
  vector<Type> rdims(3), median(nC * nT * nA), lambdas(nC * nT * nA);
  rdims << nA, nT, nC;
  for (int cc = 0; cc < nC; ++cc)
    for (int yb = 0; yb < nT; ++yb)
      for (int ag = 0; ag < nA; ++ag) {
        Type ate = intercept + yob_rw2(yb) + age_rw2(ag) + cc_vec(cc) + ccxyob(cc*nT+yb) + ccxage(cc*nA+ag);
        median[cc*nT*nA + yb*nA + ag] = 1/exp(ate) * pow(-1 + pow(0.5, -1/a_vec(cc)), -1/alpha_vec(cc));
        lambdas[cc*nT*nA + yb*nA + ag] = exp(ate);
      }
  REPORT(rdims);
  REPORT(median);
  REPORT(lambdas);
  // just a way to extract parameters without indexing headache
  REPORT(alpha_vec); REPORT(a_vec); REPORT(yob_rw2); REPORT(age_rw2); REPORT(cc_vec); 
  REPORT(ccxyob); REPORT(ccxage);
  ADREPORT(alpha_vec);
  ADREPORT(a_vec);
  ADREPORT(age_rw2_e);
  ADREPORT(yob_rw2_e);
  ADREPORT(cc_e);
  ADREPORT(ccxyob_e);
  ADREPORT(ccxage_e);
  return dll;
}
