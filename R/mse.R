#' Get Root Mean Squared Error of Matrix
#' 
#' Generate RMSEs for a given true and estimated matrix or vector
#' @param est.val The matrix or vector of estimated values
#' @param true.val The matrix or vector of true values
#' @return The RMSE value
#' @export
mse <- function(est.val, true.val) {
  
  if(class(est.val)[1]!="numeric") {
    
    p <- ncol(est.val)
    Y.mat <- matrix(rep(true.val,p),ncol=p,byrow=F)
    
    MSE <- sqrt(apply((est.val - Y.mat)^2,2,mean))
    
  } else {
    
    diff <- est.val - true.val
    MSE <- sqrt(mean(diff^2))
    
  }
  return(MSE)
} 