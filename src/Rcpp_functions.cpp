#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;


//' Ornstein Uhelenbeck Likelihood calculation
//'
//' @export
// [[Rcpp::export]]
List logl_cxx(NumericVector ThetaI, NumericVector YY, NumericMatrix XX){
  #include <cmath>  	// to use sqrt and exp
  using Eigen::Map;   	// to map input variable to an existing array of data
  using Eigen::MatrixXd;  // to use MatrixXd
  using Eigen::VectorXd;  // to use VectorXd
  using Eigen::ArrayXd;   // to use ArrayXd
  using Eigen::LLT;   	// to do the LLT decomposition
  using Eigen::Lower;
  const Map<VectorXd> Theta(Rcpp::as<Map<VectorXd> > (ThetaI));   //Map vector ThetaI to matrixXd Theta
  const Map<VectorXd> Y(Rcpp::as<Map<VectorXd> > (YY));       	//Map vector YY to vectorXd Y
  const Map<MatrixXd> X(Rcpp::as<Map<MatrixXd> > (XX));       	//Map matrix QQ to matrixXd X

  int N= Y.size();  //number of points
  double s_f2, l, s_n2;     //hyperparameters / function amplitude, characteristic length, noise amplitude, s_c
  ArrayXd dF = ArrayXd::Constant( Theta.size(),0);

  s_f2 = exp( 2.0* Theta(0) ); //exponentiate to make sure they are positive
  l   = exp(      Theta(1) );
  s_n2 = exp( 2.0* Theta(2) );
  double rr = .0;

  MatrixXd K = MatrixXd::Zero(N,N) ;        //Covariance K
  MatrixXd K_x_x = MatrixXd::Zero(N,N) ;    //Covariance K helper
  for (int i=0; i<N; i++){
    for (int j=i; j<N; j++){
      rr = (X(i,j));
      K(i,j) = s_f2 * exp(-( rr )/l) ;
    }
  }

  K_x_x = K.transpose() + K;                    // Built full matrix
  K = K_x_x + MatrixXd::Identity(N,N) * ( s_n2 - s_f2 +.000001);     // Add Noise //took out  -s_c
  LLT<MatrixXd> llt_K(K); //compute the Cholesky decomposition of the matrix
  double logLikelihood = 0.;

  if (llt_K.info()==Eigen::Success) {       // if the Cholesky decomposition exists.
    VectorXd alpha = llt_K.solve(Y);                   //Get alpha = inv(K)Y
    double A1 = 0.5*(Y.transpose()*alpha).value();         //Fit term
    double A2 = llt_K.matrixLLT().diagonal().array().log().sum();  //Complexity term
    logLikelihood = ( A1 + A2 + (N*.5)*log(2*M_PI));       //Final log-likelihood

  }
  else {
    return (List::create(Named("F") =1234567890. ));
  }
  return (List::create(Named("F") = logLikelihood));
}

