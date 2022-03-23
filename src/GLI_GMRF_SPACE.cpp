#include <TMB.hpp>
#include "ktools.hpp"

// https://en.wikipedia.org/wiki/Partial_autocorrelation_function
template <class Type>
vector<Type> to_phi(vector<Type> thetas)
{ // order 2
  vector<Type> psi(2), phi(2);
  psi[0] = 2. * exp(thetas[0]) / (1. + exp(thetas[0])) - 1.;
  psi[1] = 2. * exp(thetas[1]) / (1. + exp(thetas[1])) - 1.;
  phi[1] = psi[1];
  phi[0] = psi[0] * (1.0 - phi[1]);
  // https://github.com/kaskr/adcomp/issues/360#issuecomment-1073667612
  if (phi[1] == -1) // this should not happen sample from MVN and transformation
    phi[1] += FLT_EPSILON;
  if (phi[1] == 1 - abs(phi[0])) // this might happen
    phi[1] -= DBL_EPSILON;
  return phi;
}

template <class Type>
Eigen::SparseMatrix<Type> make_ok(matrix<Type> R, Type tau)
{
  R = tau * R.array();
  return tmbutils::asSparseMatrix(R);
}

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
  DATA_VECTOR(singlesvy);

  // priors
  DATA_VECTOR(sd_beta);
  DATA_VECTOR(palpha);
  DATA_VECTOR(p_a);
  DATA_VECTOR(pc_main);
  DATA_VECTOR(pc_interx);

  DATA_MATRIX(R_age);
  DATA_MATRIX(R_cc);
  DATA_MATRIX(R_ccxyob); // R_cc X R_yob
  DATA_INTEGER(R_ccxyob_rank);
  DATA_MATRIX(R_ccxage); // R_cc X R_age

  // Data model - log-logistic parameters
  PARAMETER(intercept);
  dll -= dnorm(intercept, sd_beta(0), sd_beta(1), true);

  // Shape
  PARAMETER_VECTOR(log_alpha_vec); 
  dll -= dnorm(log_alpha_vec, palpha(0), palpha(1), true).sum();
  vector<Type> alpha_vec = exp(log_alpha_vec);

  // Skewness *
  PARAMETER_VECTOR(a_vec_star);
  dll -= dnorm(a_vec_star, p_a(0), p_a(1), true).sum();
  // Skewness real
  vector<Type> a_vec = exp(a_vec_star - 1.1 * log_alpha_vec);

  // age rw2
  PARAMETER_VECTOR (age_rw2);
  PARAMETER        (log_age_rw2_e);
  Type age_rw2_e = exp(log_age_rw2_e);
  dll -= ktools::pc_logprec(log_age_rw2_e, pc_main[0], pc_main[1]);
  dll -= ktools::soft_zero_sum(age_rw2);
  Eigen::SparseMatrix<Type> Qage = make_ok(R_age, age_rw2_e);
  dll += density::GMRF(Qage)(age_rw2);

  // yob ARk
  // - hyper
  PARAMETER_VECTOR(pacf_vec); // * theta * //
  matrix<Type> Sigma(2,2);
  Sigma.fill(0);
  Sigma(0, 0) = Type(1); // TODO: fix upstream in TMB examples
  Sigma(1, 1) = Type(1);
  dll += density::MVNORM(Sigma)(pacf_vec);
  vector<Type> phi_yob = to_phi(pacf_vec);
  // - main params
  PARAMETER_VECTOR(yob_rw2);
  dll += density::ARk(phi_yob)(yob_rw2);
  dll -= ktools::soft_zero_sum(yob_rw2);

  // countries spatial
  PARAMETER_VECTOR  (cc_vec);
  PARAMETER         (log_cc_e);
  Type cc_e = exp(log_cc_e);
  dll -= ktools::pc_logprec(log_cc_e, pc_main[0], pc_main[1]);
  dll -= ktools::soft_zero_sum(cc_vec);
  Eigen::SparseMatrix<Type> Qcc = make_ok(R_cc, cc_e);
  dll += density::GMRF(Qcc)(cc_vec);
  
  // countries x yob interaction
  PARAMETER_VECTOR  (ccxyob);
  PARAMETER         (log_ccxyob_e);
  Type ccxyob_e = exp(log_ccxyob_e);
  prior -= ktools::pc_prec(ccxyob_e, sd_ccxyob(0), sd_ccxyob(1));
  prior -= ktools::constraint2D(ccxyob.data(), yob_rw2.size(), cc_vec.size());
  prior += density::GMRF(ktools::prepare_Q(R_ccxyob, ccxyob_e))(ccxyob);
  prior += (R_ccxyob_rank - ccxyob.size()) * log(sqrt(2*M_PI)); // ktools::GMRF would be nice

  // countries x age interaction
  PARAMETER_VECTOR  (ccxage);
  PARAMETER         (log_ccxage_e);
  Type ccxage_e = exp(log_ccxage_e);
  dll -= ktools::pc_logprec(log_ccxage_e, pc_interx[0], pc_interx[1]);
  dll -= ktools::constraint2D_singleton(ccxage.data(), singlesvy, age_rw2.size(), cc_vec.size(), true, true, true, false);
  Eigen::SparseMatrix<Type> Qccxage = make_ok(R_ccxage, ccxage_e);
  dll += density::GMRF(Qccxage)(ccxage);

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
  REPORT(intercept);
  REPORT(alpha_vec); REPORT(a_vec); 
  REPORT(yob_rw2); REPORT(age_rw2); REPORT(cc_vec); 
  REPORT(ccxyob); REPORT(ccxage); 
  REPORT(rdims);
  REPORT(lambdas); REPORT(median);
  return dll;
}
