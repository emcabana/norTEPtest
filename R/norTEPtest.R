#'@title Normality tests based on TEP's
#'
#' @description Performs tests of normality consistent
#' and focused to detect skewness (type='skew', the default)
#' or kurtosis (type='kurt')
#'
#' @details The tests statistics are quadratic forms in the vector of evaluations
#' of the normalized Hermite polynomials at the scaled sample points,
#' introduced by the authors in
#' \emph{Tests of normality based on Transformed Empirical Processes},
#' Methodology and Computing in Applied Probability, Kluwer Academic
#' Publishers,vol.5, pp 309-335, 2003.
#' @param X Vector containing the elements of the sample
#' @param type Character specifying the alternative of focusing:
#'   \itemize{
#'     \item `"skew"` — Test focused to detect skewness
#'     \item `"kurt"` — Test focused to detect changes in kurtosis
#'   }
#' @param l Integer, Number of Hermite polynomials evaluated at the sample points -1
#' @param nrep Integer, Number of Monte Carlo evaluations to estimate the p-value
#' (defalts to 1000)
#'
#' @return A list with
#'   \itemize{
#'     \item `"pval"` — The p-value of the sample estimated by MC
#'     \item `"type"` — The alternative of focusing
#'     \item `"est"` — The value of the test statistic 
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- rnorm(100)
#' TEPnor(X)
#' }
#'
#' @export
TEPnor <- function(X,type='skew',l=10,nrep=1000){
	Qdat <- TEPest(X,type=type,l=l)
	Qrep <- sapply(1:nrep,function(i){Xsim=rnorm(length(X))
		TEPest(Xsim,type=type,l=l)})
		pval <- mean(Qrep>rep(Qdat,length(Qrep)))
		return(list(pval=pval,type=type,est=Qdat))
}
#'
#' @title Normalized Hermite polynomials
#' @description Computes the normalized Hermite polynomials
#' of degrees 0 to l at the entries of the vector x_vec
#' @param x_vec Numeric vector
#' @param l Integer
#' @return A matrix with rows corresponding to degrees 0..l
#' @export
eval_hvec_stable <- function(x_vec, l) {
  n  <- length(x_vec)
  Hs <- matrix(0, l + 1, n)
  Hs[1, ] <- 1
  if (l >= 1) Hs[2, ] <- x_vec
  for (k in seq_len(l - 1)) {
    Hs[k + 2, ] <- (x_vec * Hs[k + 1, ] - sqrt(k) * Hs[k, ]) / sqrt(k + 1)
  }
  Hs
}
# Funcion 2: matriz C por Monte Carlo paralelo
#' @title Monte Carlo approximation of C matrix
#' @description Computes the covariance matrix C by parallel Monte Carlo simulation
#' @param l Integer, maximum Hermite degree
#' @param N_per_core Number of samples per core
#' @param block_size Size of each block
#' @param n_cores Number of cores used
#' @param seed Random seed
#' @return Symmetric matrix C of dimension (l+1) x (l+1)
#' @export
C_matrix_MC_par <- function(l, N_per_core = 1e8, block_size = 1e5,
                             n_cores = detectCores()-1, seed = 42) {
  n_blocks <- N_per_core / block_size

  C_list <- mclapply(seq_len(n_cores), function(b) {
    set.seed(seed + b)
    C_loc <- matrix(0, l+1, l+1)
    for (i in seq_len(n_blocks)) {
      s  <- rnorm(block_size); t <- rnorm(block_size)
      u  <- pnorm(s);          v <- pnorm(t)
      g  <- (abs(u - v) - 0.5)^2 + 0.25
      Hs <- eval_hvec_stable(s, l)
      Ht <- eval_hvec_stable(t, l)
      C_loc <- C_loc + (Hs * rep(g, each = l+1)) %*% t(Ht)
    }
    C_loc
  }, mc.cores = n_cores)

  N_total <- N_per_core * n_cores
  C <- Reduce("+", C_list) / N_total
  C <- (C + t(C)) / 2

  # Anular i+j impar
  mask <- outer(0:l, 0:l, function(i,j) (i+j) %% 2 == 0)
  C <- C * mask

  C[1,] <- 0; C[,1] <- 0; C[1,1] <- 1/3
  C
}
#' @title TEP-based CvM statistic
#' @description Computes the TEP-based Cramér–von Mises statistic
#' @param X Numeric vector, the sample
#' @param type Character specifying the alternative of focusing
#' @param l Integer, number of Hermite polynomials evaluated at the sample points -1
#' @return Numeric value of the test statistic
#' @export
TEPest=function(X,type='skew',l=10){
	Y=X-mean(X);Y=Y/sqrt(mean(Y^2))
	herm <- rowSums(eval_hvec_stable(Y,l+3)[4:(l+4),])/sqrt(length(X))
	if(l<=10)Cmat=Cdat[1:(l+1),1:(l+1)] else Cmat<-C_matrix_MC_par(l)
	if(type=='kurt'){Cmat=Cmat[c(1:3,5,4,6:(l+1)),c(1:3,5,4,6:(l+1))]
		herm=c(herm[2],herm[1],herm[3:length(herm)])}
	Q=t(herm)%*%Cmat%*%herm	
}

