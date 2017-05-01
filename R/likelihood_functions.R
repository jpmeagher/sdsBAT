#' Phylogenetic Likelihood
#'
#' Calculates the Phylogenetic Ornsetin-Uhelenbeck Process negative log
#' likelihood, given the process hyperparameters, a phylogenetic tree, and some
#' data.
#'
#' @param dat Observations from the tips of the phylogenetic tree
#' @inheritParams phylo_ou_cov
#'
#' @return A positive value which is the negative log likelihood for the data.
pou_logl_slow <- function(ln_hyperparameters, phylogenetic_tree, observations){
  K <- pou_covariance(ln_hyperparameters, phylogenetic_tree)

  ntips <- length(phylogenetic_tree$tip.label)

  C <- chol(K)
  K_inv <- chol2inv(C)
  y <- matrix(observations, nrow = ntips, ncol = 1)

  A.1 <- -.5*((t(y) %*% K_inv) %*% y)
  A.2 <- -.5*log(prod(diag(C))^2)
  A.3 <- -(ntips/2) * log(2*pi)
  lik <- -(A.1 + A.2  + A.3)

  lik
}


pou_logl_cxxf <- inline::cxxfunction(signature(ThetaI = "vector",YY = "vector", XX= "matrix"),
  '
  #include <cmath>  	// to use sqrt and exp
  using Eigen::Map;   	// to map input variable to an existing array of data
  using Eigen::MatrixXd;  // to use MatrixXd
  using Eigen::VectorXd;  // to use VectorXd
  using Eigen::ArrayXd;   // to use ArrayXd
  using Eigen::LLT;   	// to do the LLT decomposition
  using Eigen::Lower; 	// to get the lower triangular view
  const Map<VectorXd> Theta(Rcpp::as<Map<VectorXd> > (ThetaI));   //Map vector ThetaI to matrixXd Theta
  const Map<VectorXd> Y(Rcpp::as<Map<VectorXd> > (YY));       	//Map vector YY to vectorXd Y
  const Map<MatrixXd> X(Rcpp::as<Map<MatrixXd> > (XX));       	//Map matrix QQ to matrixXd X
  //const Map<VectorXd> B(Rcpp::as<Map<VectorXd> > (BB));       	//Map vector BB to vectorXd B

  using namespace std;

  int N= Y.size();  //number of points
  double s_f2, l, s_n2;     //hyperparameters / function amplitude, characteristic length, noise amplitude, s_c
  ArrayXd dF = ArrayXd::Constant( Theta.size(),0);

  s_f2 = exp( 2*Theta(0) );  //exponentiate to make sure they are positive
  l   = exp( Theta(1) );
  s_n2 = exp( 2*Theta(2) );
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
  return (List::create(Named("F") = logLikelihood));',
  plugin = "RcppEigen")



pou_logl_cxx <- function(ln_hyperparameters, phylogenetic_tree, observations){
  distances <- ape::cophenetic.phylo(phylogenetic_tree)

  storage.mode(distances) <- "double";        #Save stuff as doubles
  storage.mode(observations) <- "double";
  storage.mode(ln_hyperparameters) <- "double";

  return ( pou_logl_cxxf(ln_hyperparameters, observations, distances));   #LogLikelihood
}

#' Fast Phylogenetic Likelihood
#'
#' Using the inline and RcppEigen packages performs a fast calculation of the
#' Phylogenetic Ornsetin-Uhelenbeck Process negative log likelihood, given the
#' process hyperparameters, a phylogenetic tree, and some data. The C++ code was
#' taken from code accompanying Hajipantelis et al.
#'
#' @inheritParams phylo_ou_likelihood
#' @inheritParams phylo_ou_cov
#'
#' @return A positive value which is the negative log likelihood for the data.
#'
#' @source P. Z. Hadjipantelis, N. S. Jones, J. Moriarty, D. A. Springate, and
#'   C. G. Knight, Function-valued traits in evolution, Journal of The Royal
#'   Society Interface. 10(82), 20121032 (2013). \url{https://github.com/fpgpr}
pou_logl_fast <- function(ln_hyperparameters, phylogenetic_tree, observations){
  if(!requireNamespace("inline", quietly = TRUE)){
    stop("inline needed for this function to work. Please install it, or use pou_logl_slow()",
      call. = FALSE)
  }
  if(!requireNamespace("RcppEigen", quietly = TRUE)){
    stop("RcppEigen needed for this function to work. Please install it, or use pou_logl_slow()",
      call. = FALSE)
  }
  distances <- ape::cophenetic.phylo(phylogenetic_tree)

  storage.mode(distances) <- "double";        #Save stuff as doubles
  storage.mode(observations) <- "double";
  storage.mode(ln_hyperparameters) <- "double";

  return ( pou_logl_cxxf(ln_hyperparameters, observations, distances)$F);   #LogLikelihood
}


