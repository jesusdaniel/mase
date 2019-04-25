#' Omnibus embedding for only two adjacency matrices
#' @param A1 adjacency matrix
#' @param A2 adjacency matrix
#' @param d number of embedding dimensions of the omnibus matrix. If NULL, dimension is chosen automatically
#' @return A list containing the estimated latent positions of each matrix
OMNI_embedding <- function(A1, A2, d = NULL) {
  n <- ncol(A1)
  M <- rbind(cbind(A1, (A1+ A2)/2),
             cbind((A1+ A2)/2, A2))
  Z <- ase(M, d)
  X1 <- Z[1:n, ]
  X2 <- Z[(n+1):(2*n), ]
  return(list(X1 = X1, X2 = X2))
}


#' Omnibus embedding for a list of adjacency matrices
#' @param Adjacencies list of matrices
#' @param d number of embedding dimensions of the omnibus matrix. If NULL, dimension is chosen automatically
#' @return A list containing the estimated generalised latent positions of each matrix
OMNI_matrix <- function(Adjacencies, d = NULL) {
  m = length(Adjacencies)
  n = ncol(Adjacencies[[1]])
  Omni = compute_OMNIBUS_mat(Adjacencies)
  Z = g.ase(Omni, d)
  Xhats <- lapply(1:m, function(i) list(X=Z$X[((i-1)*n + 1):(i*n), ],
                                        D = Z$D))
  return(Xhats)
}





bootstrap_omni <- function(X.gase, B = 1000) {
  d <- ncol(X.gase$X)
  n <- nrow(X.gase$X)
  P1 <- P_from_g.ase(X.gase)
  P1 = P1 * (P1>0)
  P1 = P1*(P1 <1) + 1* (P1>1)
  Ab1 <- 0*P1
  Ab2 <- 0*P1
  bootstrap <- sapply(1:B, function(i) {
    Ab1[upper.tri(Ab1)] <- 1*(runif(sum(upper.tri(P1))) < P1[upper.tri(P1)])
    Ab2[upper.tri(Ab1)] <- 1*(runif(sum(upper.tri(P1))) < P1[upper.tri(P1)])
    Ab1 = (Ab1 + t(Ab1))
    Ab2 = (Ab2 + t(Ab2))
    omni <- OMNI_embedding(Ab1, Ab2, d)
    return(norm(omni$X1 - omni$X2, "F"))
  })
  return(bootstrap)
}

bootstrap_omni_from_P <- function(P, d,  B = 1000) {
  n <- nrow(P)
  P1 <- P
  P1 = P1 * (P1>0)
  P1 = P1*(P1 <1) + 1* (P1>1)
  Ab1 <- 0*P1
  Ab2 <- 0*P1
  bootstrap <- sapply(1:B, function(i) {
    Ab1[upper.tri(Ab1)] <- 1*(runif(sum(upper.tri(P1))) < P1[upper.tri(P1)])
    Ab2[upper.tri(Ab1)] <- 1*(runif(sum(upper.tri(P1))) < P1[upper.tri(P1)])
    Ab1 = (Ab1 + t(Ab1))
    Ab2 = (Ab2 + t(Ab2))
    omni <- OMNI_embedding(Ab1, Ab2, d)
    return(norm(omni$X1 - omni$X2, "F"))
  })
  return(bootstrap)
}


compute_OMNIBUS_mat <- function(Adjacencies) {
  m = length(Adjacencies)
  n = ncol(Adjacencies[[1]])
  Omni1 <- lapply(Adjacencies, function(A) (Reduce(function(U,j) cbind(U, A), c(list(A), lapply(1:(m-1), function(i) i) ))))
  Omni2 <- lapply(Adjacencies, function(A) (Reduce(function(U,j) rbind(U, A), c(list(A), lapply(1:(m-1), function(i) i) ))))
  Omni = (Reduce(rbind, Omni1) + Reduce(cbind, Omni2))/2
  return(Omni)
}