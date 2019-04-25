sample_from_P <- function(P) {
  n = ncol(P)
  A = Matrix(0, n, n)
  A[upper.tri(A)] <- 1*(runif(sum(upper.tri(P))) < P[upper.tri(P)])
  A = A + t(A)
  return(A)
}
