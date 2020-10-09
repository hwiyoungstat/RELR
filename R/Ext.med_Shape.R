#' Ext.med_Shape
#'
#' @param X sample data on shape space
#'
#' @return geometric median of planar shape
#' @export

Ext.med_Shape <- function(X){
  K <- dim(X)[1]
  N <- dim(X)[3]

  Old <- matrix(0,K,K)
  alpha <- 1

  Ite <- 100
  Tol <- 1e-5

  # Distance function
  Distance <- function(A,B){
    return(Re(sqrt(tr((A-B)%*%t(Conj(A-B))))))
  }

  # G function
  G <- function(A,Old,dist){
    Normalized <- sapply(1:dim(A)[3], function(i) (Old-A[,,i])/dist[i] , simplify="array")
    Num <- apply(Normalized, c(1,2), sum)
    Deno <- sum(1/dist)

    G <- Num/Deno
    return(G)
  }

  # Loop

  for(j in 1:Ite){

  D <- sapply(1:N, function(i) Distance(X[,,i], Old))
  Gk <- G(X,Old,D)
  New <- Old-alpha*Gk

  if (Distance(New,Old)^2 < Tol) break
  Old <- New

  print(j)
  }

  return(New)
}







