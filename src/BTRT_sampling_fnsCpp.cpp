#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include <RcppGSL.h>

// This is a script with functions for performing the MCMC for the BTRTucker model

// Rcpp::NumericVector BTRT_draw_gamCpp(Rcpp::NumericMatrix eta,
//                                      Rcpp::NumericMatrix Sig0,
//                                      Rcpp::NumericVector mu_gam,
//                                      Rcpp::NumericVector y_til,
//                                double sig_y2) {
//   Rcpp::NumericMatrix Sig0_inv = Sig0.;
//   Rcpp::NumericMatrix Sig_gam_inv = eta.transpose() * eta / sig_y2;
//   Sig_gam_inv += Sig0_inv;
//   Rcpp::NumericMatrix Sig_gam = Sig_gam_inv.solve();
//   Rcpp::NumericVector mu_gam_p = y_til.transpose() * eta;
//   mu_gam_p /= sig_y2;
//   mu_gam_p += mu_gam.transpose() * Sig0_inv;
//   mu_gam_p *= Sig_gam;
//   Rcpp::NumericVector Z = rnorm(eta.cols());
//   Rcpp::NumericVector gam = Sig_gam.Cholesky() * Z;
//   gam += mu_gam_p.transpose();
//   return gam;
// }
