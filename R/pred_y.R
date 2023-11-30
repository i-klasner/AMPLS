#' Find Predicted Fitted Values
#' 
#' Generate a list of Y matrices for each number of principle components
#' @param testX The matrix of predictor variables
#' @param testY The matrix of response variables
#' @param coeffs The list of coefficient matrices to generate fitted values
#' @return A list of fitted values.
#' @examples
#' data(nutrimouse);
#' simpls_B <- simpls(nutrimouse$genotype, nutrimouse$lipid, cov(nutrimouse$genotype))$coefficients;
#' fitted <- pred_y(nutrimouse$genotype, nutrimouse$lipid, simpls_B);
#' @export
pred_y <- function(testX,testY,coeffs){
  # 30x120, 30x21, 120x21x29
  n <- nrow(testX) # 30
  p <- ncol(testX) # 120
  q <- ncol(testY) # 21
  a <- dim(coeffs)[3] # 29
  
  Y <- array(rep(0, n*q*a),c(n,q,a)) # 30 x 21 x 29
  for (i in 1:q) { # iterating through response variables
    B <- unlist(coeffs[,i,])  # matrix of the coefficients, 120x29
    
    mean_X <- matrix(colMeans(testX),nrow=1)
    rep_mean_Y <- matrix(colMeans(testY)[i],ncol=a)
    
    b0 <- rep_mean_Y - mean_X%*%B  ## estimates for intercept
    B0 <- matrix(c(b0),nrow=n,ncol=a,byrow=T)  ## a matrix of coefficients for the intercept
    
    fitted_val <-  B0 + testX%*%B
    Y[,i,] <- fitted_val
  }
  
  return(Y);
}