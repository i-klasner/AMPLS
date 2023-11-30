# R Package Practice

#' Proposed SIMPLS Method
#'
#' Find matrix of coefficients and fitted values from given data.
#' @param X The matrix of predictor variables
#' @param Y The matrix of response variables
#' @param dd The covariance matrix
#' @return A list of the following: the coefficient matrix from the maximum number of principle components
#' a scores matrix, leverages, the x variance and y variance matrices, the number of principle 
#' components used, a list of coefficient matrices by number of principle components used, the P
#' matrix, R matrix, V matrix, and a list of fitted values by number of principle components used.
#' @examples
#' data(nutrimouse);
#' simpls_B <- simpls(nutrimouse$genotype, nutrimouse$lipid, cov(nutrimouse$genotype))$coefficients;
#' fitted <- pred_y(nutrimouse$genotype, nutrimouse$lipid, simpls_B);
#' @export
##### SIMPLS Code #####
simpls <- function(X,Y,dd){
  n <- nrow(X) # sample size
  p <- ncol(X) # number of predictor variables
  resp <- ncol(Y) # number of response variables
  Y0 <- Y-rep(colMeans(Y), each=n) # center Y
  if(n > p) {
    a <- p
  } else {
    a <- n-1
  }
  coeffs <- array(rep(0, p*resp*a),c(p,resp,a)) # coefficients for number of PC
  
  cov_dd <- dd
  S_ww <- t(X)%*%X
  S_wy <- t(X)%*%Y0
  
  W <- X
  Omega.u <- ginv(cov_dd)
  tol <- 0.001
  M.iter <- 1000
  for(iter in 1:M.iter){ 
    if(iter==1){
      Omega.x <- ginv(S_ww/n)
    }else{
      Omega.x <- ginv(prev.Sigma.x)
    }
    Omega.x <- (Omega.x+t(Omega.x))/2 # ensure Omega.x is symmetric
    Lambda <- (Omega.x + Omega.u) 
    inv.Lambda <- solve(Lambda) 
    
    for(i in 1:n){
      w.vec <- matrix(W[i,],p,1)
      tmp.eval1 <- inv.Lambda%*%Omega.u%*%w.vec
      tmp.eval2 <- tmp.eval1%*%t(tmp.eval1)
      
      if(i==1){
        B.mat <- tmp.eval2
      }else{
        B.mat <- B.mat + tmp.eval2
      }
      
      if(i==n){
        Sigma.x <- (inv.Lambda + B.mat/n)
      }
    }
    
    if((iter > 1) && (max(abs(prev.Sigma.x - Sigma.x)) < tol)){
      break;
    } else{
      prev.Sigma.x <- Sigma.x
    }
  }  
  
  ### M-Step ###
  
  D_w <- diag(sqrt(pmax(0,eigen(S_ww)$values)),p,p)
  diag(D_w)[-c(which(diag(D_w) > 10^(-3)))]=0
  Q_w <- eigen(S_ww)$vector
  S_wh <- Q_w%*%D_w%*%t(Q_w)
  D_x <- diag(sqrt(eigen(n*Sigma.x)$values),p,p) 
  diag(D_x)[-c(which(diag(Re(D_x)) > 10^(-3)))]=0
  Q_x <- eigen(n*Sigma.x)$vector
  S_xh <- Q_x%*%D_x%*%t(Q_x)
  inv.D_w <- diag(0,p,p)
  diag(inv.D_w)[which(diag(D_w)>10^(-3))] <- 1/diag(D_w)[which(diag(D_w)>10^(-3))]
  inv.D_x <- diag(0,p,p)
  diag(inv.D_x)[which(diag(Re(D_x))>10^(-3))] <- 1/diag(D_x)[which(diag(Re(D_x))>10^(-3))]
  inv.S_wh <- Q_w%*%inv.D_w%*%t(Q_w)
  inv.S_xh <- Q_x%*%inv.D_x%*%t(Q_x)
  X <- W%*%inv.S_wh%*%S_xh  # calibration of X
  Y0 <- W%*%inv.S_wh%*%inv.S_xh%*%S_wy
  
  S <- t(X)%*%Y0
  
  for (f in 1:a) {
    eigen <- eigen(t(S)%*%S)
    comp <- which.max(eigen$values) # find the index of the largest eigenvalue in the vector of values
    q <- eigen$vectors[,comp] # Y factor weights
    r <- S%*%q # X factor weights
    tt <- X%*%r # X factor scores
    tt <- tt-mean(tt) # center scores
    ttnorm <- norm(tt,type='2') # norm of t
    tt <- tt/ttnorm # normalize scores
    r <- r/ttnorm[[1]] # normalize weights with norm of t
    p <- t(X)%*%tt
    q <- t(Y0)%*%tt
    u <- Y0%*%q
    v <- p # orthogonal loadings
    if(f > 1) {
      v <- v - V%*%(t(V)%*%p)
      u <- u - TT%*%(t(TT)%*%u)
    }
    v <- v/norm(v,type='2') # normalize orthogonal loadings
    S <- S - v%*%(t(v)%*%S)
    if (f == 1) { #if first iteration, create new vectors
      R <- r
      TT <- tt
      P <- p
      Q <- q
      U <- u
      V <- v
    } else { #otherwise, add columns to matrices
      R <- cbind(R,r)
      TT <- cbind(TT,tt)
      P <- cbind(P,p)
      Q <- cbind(Q,q)
      U <- cbind(U,u)
      V <- cbind(V,v)
    }
    #find coefficients
    B <- R%*%t(Q) # regression coefficients
    coeffs[,,f] <- B
  }
  h <- diag(TT%*%t(TT)) + length(X) # leverages
  varX <- diag(t(P)%*%P)/(length(X)-1) # variance explained for X
  varY <- diag(t(Q)%*%Q)/(length(X)-1) # variance explained for Y
  
  Y <- pred_y(X, Y, coeffs) # get fitted values for training set
  
  return(list(this_coeff=B,
              scores=TT,
              leverages=h,
              xvariance=varX,
              yvaraince=varY,
              numpc=f,
              coefficients=coeffs,
              P=P, 
              R=R,
              V=V,
              fitted.values=Y))
}