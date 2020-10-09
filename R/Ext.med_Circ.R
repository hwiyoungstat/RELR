#' Ext.med_Circ
#'
#' @param X sample data on unit circle
#'
#' @return geometric median on unit circle
#' @export

Ext.med_Circ <- function(X){
  N <- dim(X)[1]
  K <- dim(X)[2]

  Old <- colMeans(X)/sqrt(sum(colMeans(X)^2)) # Starting Extrinsic mean
  alpha <- 1

  Ite <- 100
  Tol <- 1e-7

  # Distance function
  Distance <- function(A,B){
    return(sqrt(sum((A-B)^2)))
  }

  # G function
  G <- function(A,Old,dist){
    Normalized <- sapply(1:dim(A)[1], function(i) (Old-A[i,])/dist[i] , simplify="matrix")
    Num <- rowSums(Normalized)
    Deno <- sum(1/dist)

    G <- Num/Deno
    return(G)
  }

  # Loop

  for(j in 1:Ite){

    D <- sapply(1:N, function(i) Distance(X[i,], Old))
    Gk <- G(X,Old,D)
    New <- Old-alpha*Gk
    New <- New/sqrt(sum(New^2))

    if (Distance(New,Old)^2 < Tol) break
    Old <- New

    print(j)
  }

  return(New)
}

