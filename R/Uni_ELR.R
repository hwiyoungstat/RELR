#' Uni_ELR
#'
#' @param gird Evaluation point
#' @param X Covariate
#' @param Response  shape (Matrix N*K)
#' @param H bandwidths (column)
#'
#'
#' @return h optimal bandwidth
#' @return RMSE with observed
#' @return RMSE_True with true value
#' @export
#'
Uni_ELR <- function(grid,X,Response,H){

  K <- dim(Response)[2]
  N <- length(X) # Number of data
  N.eval <- length(grid) # Number of Evaluation points
  NH <- length(H)


  J <- function(x){return(x%*%t(Conj(x)))} # VW embedding
  Y <- sapply(1:N, function(i) J(Response[i,]),simplify="array") # Embeded shape


  # Define preshape function
  pre <- function(x){
    trans <- x-mean(x)
    norm <- Re(trans%*%Conj(trans))
    return(trans/sqrt(norm))
  }




  # K-fold Cross Valitaion (5)

  RMSE <- rep(0,NH)

  F <- 5 # Number of replication

  for(f in 1:F){

  CV.idx <- sample(N,round(4*N/5),replace=FALSE)
  VA.idx <- seq(1:N)[-CV.idx]
  CV.N <- length(CV.idx)
  VA.N <- N-CV.N
  CV.X <- X[CV.idx]
  VA.X <- X[-CV.idx]
  CV.Y <- Y[,,CV.idx]
  VA.Y <- Y[,,-CV.idx]
  Est <- array(0,c(K,VA.N,NH))
  Validate.Response <- Response[VA.idx,]

  for(b in 1:NH){
    h <- H[b]

  # Calculate the kernel Weight

  W <- sapply(1:VA.N, function(i) dnorm((CV.X-VA.X[i])/h)/h,simplify="matrix") # Calculate Weight (Row = CV for model, COl=VAlidataion(evaluation point))
  W.Sum <- colSums(W)
  Weight <- sapply(1:VA.N, function(i) W[,i]/W.Sum[i],simplify="matrix")


  tr <- sapply(1:VA.N, function(j) sapply(1:CV.N, function(i) Weight[i,j]*CV.Y[,,i],simplify="array"),simplify="array")
  Est_J <- sapply(1:VA.N, function(i) apply(tr[,,,i],c(1,2), sum), simplify="array") # Estimated embeded matrix
  Est.temp <- sapply(1:VA.N, function(i) eigen(Est_J[,,i])$vectors[,1],simplify="array")
  Est[,,b] <- sapply(1:VA.N, function(i) pre(Est.temp[,i]),simplify="matrix") # preshape array of estimated embeded matrix



  Temp <- sapply(1:VA.N, function(i) procdist(Est[,i,b], Validate.Response[i,],type="full",reflect=FALSE)^2)
  RMSE[b] <- RMSE[b] + sum(Temp)

  }

  }

  # Fitting

  h <- H[which.min(RMSE)]


  # Calculate the kernel Weight

  W <- sapply(1:N.eval, function(i) dnorm((X-grid[i])/h)/h,simplify="matrix") # Calculate Weight (Row = CV for model, COl=VAlidataion(evaluation point))
  W.Sum <- colSums(W)
  Weight <- sapply(1:N.eval, function(i) W[,i]/W.Sum[i],simplify="matrix")


  tr <- sapply(1:N.eval, function(j) sapply(1:N, function(i) Weight[i,j]*Y[,,i],simplify="array"),simplify="array")
  Est_J <- sapply(1:N.eval, function(i) apply(tr[,,,i],c(1,2), sum), simplify="array") # Estimated embeded matrix
  Est.temp <- sapply(1:N.eval, function(j) eigen(Est_J[,,j])$vectors[,1],simplify="array")
  Est <- sapply(1:N.eval, function(i) pre(Est.temp[,i]),simplify="matrix") # preshape array of estimated embeded matrix

  # Calculate the RMSE with Observed values

  Temp <- sapply(1:N.eval, function(i) procdist(Est[,i], Response[i,],type="full",reflect=FALSE)^2)
  RMSE <- sqrt(mean(Temp))


  return(list(Est,h,RMSE))


}




