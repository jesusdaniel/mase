# Author: Jesus Arroyo & Shangsi Wang
#' Function to perform joint embedding of graphs into d dimensions.
#' @param Adj_list a list of adjacency matrices with the same size n x n
#' @param d number of joint embedding dimensions
#' @param K number of embedding dimensions of each graph
#' @param par whether to run in parallel (TRUE/FALSE)
#' @param numpar number of clusters for parallel implmentation
#' @return A list containing the estimated subspace V and a list R with the individual parameters for each graph
joint_latent_positions <- function(Adj_list, d, K=d, ASE = FALSE,
                                   par = FALSE, numpar = 12) {
  if(par) {
    require(parallel)
    cl <- makeCluster(numpar)
    if(ASE) {
      clusterExport(cl = cl, varlist = list("ase"))
      latpos.list <- parLapply(cl = cl, Adj_list, ase, d = K)
    }else{
      clusterExport(cl = cl, varlist = list("eig"))
      latpos.list <- parLapply(cl = cl, Adj_list, eig, d = K)
    }
    stopCluster(cl)
  } else {
    if(ASE) {
      latpos.list <- lapply(Adj_list, ase, d = K)
    }else{
      latpos.list <- lapply(Adj_list, eig, d = K)
    }
  }
  V_all  <- Reduce(cbind, latpos.list)
  require(rARPACK)
  jointsvd <- svd(V_all, d)
  R <- project_networks(Adj_list, jointsvd$u)
  return(list(V = jointsvd$u, sigma = jointsvd$d, R = R))
}


ase <- function(A, d) {
  require(rARPACK) 
  eig <- eigs_sym(as.matrix(A), k = d)
  X <- eig$vectors %*% diag(sqrt(abs(eig$values)), nrow = d)
  return(X)
}

eig <- function(A, d) {
  require(rARPACK) 
  eig <- eigs_sym(as.matrix(A), k = d)$vectors
  return(eig)
}


graph_distance <- function(A, B) {
  return(sum(abs(A-B)))
}

cluster_positions <- function(V, K) {
  require(mclust)
  clus <- Mclust(V, G = K)
  return(clus$classification)
}


plot.3d <- function(V, col = NULL, title = NULL) {
  require(plotly)
  if(is.null(col)) {
    pl <- plot_ly(x = V[, 1], y = V[,2], z = V[,3], main = "hi")
  }else {
    pl <- plot_ly(x = V[, 1], y = V[,2], z = V[,3], color = col)  
  }
  pl %>%
    layout(title = title)
}


project_networks <- function(Adj_list, V) {
  require(Matrix)
  lapply(Adj_list, function(A) crossprod(crossprod(A, V), V))
}

dissimilarity.graphs <- function(Adj_list) {
  D <- sapply(Adj_list, function(A) 
    sapply(Adj_list, function(B) graph_distance(A, B)))
}

fit_joint_sbms <- function(Adj_list, community_memberships) {
  K <- length(unique(community_memberships))
  N <- ncol(Adj_list[[1]])
  Z <- Matrix(0, nrow = N, ncol = K)
  for(k in 1:K) {
    Z[which(community_memberships==k), k] <- 1
  }
  c_size <- colSums(Z)
  if(length(c_size)> 1) {
    cell_sizes <- tcrossprod(c_size) - diag(c_size)
  }else{
    cell_sizes <- c_size^2 - c_size
  }
  Bs <- lapply(Adj_list, function(A)  
    crossprod(crossprod(A, Z), Z) / cell_sizes)
  return(list(Z = Z, Bs = Bs))
}

subspace_distance <- function(V1, V2) {
  U1 = svd(V1)$u
  U2 = svd(V2)$u
  sqrt(ncol(U1) - sum(svd(crossprod(U1, U2))$d))
}
