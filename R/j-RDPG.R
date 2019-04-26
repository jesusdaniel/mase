#' Function to perform multiple adjacency spectral embedding
#' @param Adj_list a list of adjacency matrices with the same size n x n
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param K maximum number of embedding dimensions of each graph
#' @param par whether to run in parallel (TRUE/FALSE)
#' @param numpar number of clusters for parallel implmentation
#' @return A list containing the estimated subspace V and a list R with the individual parameters for each graph
joint_latent_positions <- function(Adj_list, d = NULL, K=d, ASE = FALSE,
                                   par = FALSE, numpar = 12) {
  if(par) {
    require(parallel)
    cl <- makeCluster(numpar)
    if(ASE) {
      clusterExport(cl = cl, varlist = list("ase", "getElbows"))
      latpos.list <- parLapply(cl = cl, Adj_list, ase, d = K)
    }else{
      clusterExport(cl = cl, varlist = list("eig_embedding", "getElbows"))
      latpos.list <- parLapply(cl = cl, Adj_list, eig_embedding, d = K)
    }
    stopCluster(cl)
  } else {
    if(ASE) {
      latpos.list <- lapply(Adj_list, ase, d = K)
    }else{
      latpos.list <- lapply(Adj_list, eig_embedding, d = K)
    }
  }
  V_all  <- Reduce(cbind, latpos.list)
  require(rARPACK)
  jointsvd <- svd(V_all)
  if(is.null(d)) {
    hist(sapply(latpos.list, ncol))
    d = getElbows(jointsvd$d, plot = TRUE)[1]
  }
  V = jointsvd$u[, 1:d, drop = FALSE]
  R <- project_networks(Adj_list, V)
  return(list(V = V, sigma = jointsvd$d, R = R))
}



#' Adjacency spectral embedding
#' @param A adjacency matrix
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param d.max maximum number of embedding dimensions to choose (if d is given, this number is ignored)
#' @param diag.augment whether to do diagonal augmentation (TRUE/FALSE)
#' @return A matrix with n rows and d columns containing the estimated latent positions
ase <- function(A, d = NULL, d.max = ncol(A), diag.augment = TRUE) {
  require(rARPACK)
  # Diagonal augmentation
  if(diag.augment & sum(abs(diag(A))) == 0) {
    deg = colSums(A)
    n = ncol(A)
    diag(A) = deg / (n-1)
  }
  if(is.null(d)) {
    eig <- eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x =  abs(eig$values), decreasing = TRUE)
    #d = getElbow_GMM(vals)
    d = getElbows(vals, plot = F)[1]
    selected.eigs <- which(abs(eig$values) >= vals[d])
    X <- eig$vectors[,selected.eigs] %*% diag(sqrt(abs(eig$values[selected.eigs])), nrow = d)
    return(X)
  } else {
    eig <- eigs(as(A, "dgeMatrix"), k = d)
    X <- eig$vectors %*% diag(sqrt(abs(eig$values)), nrow = d)
    return(X) 
  }
}

#' Generalised adjacency spectral embedding
#' @param A adjacency matrix
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param d.max maximum number of embedding dimensions to choose (if d is given, this number is ignored)
#' @param diag.augment whether to do diagonal augmentation (TRUE/FALSE)
#' @return A list containing a matrix with n rows and d columns representing the estimated latent positions, and the estimated
#' indefinite orthogonal projection matrix
g.ase <- function(A, d = NULL, d.max = ncol(A), diag.augment = T) {
  require(rARPACK) 
  if(is.null(d)) {
    if(diag.augment & sum(abs(diag(A))) == 0) {
      deg = colSums(A)
      n = ncol(A)
      diag(A) = deg / (n-1)
    }
    eigv <- eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x =  abs(eigv$values), decreasing = TRUE)
    #d = getElbow_GMM(vals)
    d = getElbows(vals, plot = F)[1]
    selected.eigs <- which(abs(eigv$values) >= vals[d])
    X <- eigv$vectors[,selected.eigs] %*% diag(sqrt(abs(eigv$values[selected.eigs])), nrow = d)
    D <- sign(eigv$values[selected.eigs])
  } else{
    eig <- eigs(as(A, "dgeMatrix"), k =  d)
    X <- eig$vectors %*% diag(sqrt(abs(eig$values)), nrow = d)
    D <- sign(eig$values) 
  }
  return(list(X=X, D=D))
}

#' Obtain expected adjacency matrix from generalised RDPG
#' @param g_ase output of g.ase
#' @return A matrix with expected adjacency 
P_from_g.ase <- function(g_ase) {
  g_ase$X %*% diag(g_ase$D) %*% t(g_ase$X)
}

#' Unscaled adjacency spectral embedding (top eigenvectors)
#' @param A adjacency matrix
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param d.max maximum number of embedding dimensions to choose (if d is given, this number is ignored)
#' @param diag.augment whether to do diagonal augmentation (TRUE/FALSE)
#' @return A list containing a matrix with n rows and d columns representing the estimated latent positions, and the estimated
#' indefinite orthogonal projection matrix
eig_embedding <- function(A, d = NULL, d.max = ncol(A), diag.augment = FALSE) {
  require(rARPACK)
  n <- ncol(A)
  if(diag.augment & sum(abs(diag(A))) == 0) {
    deg = colSums(A)
    n = ncol(A)
    diag(A) = deg / (n-1)
  }
  if(is.null(d)) {
    eig <- eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x =  abs(eig$values), decreasing = TRUE)#[1:sqrt(n)]
    #d = getElbow_GMM(vals)
    d = getElbows(vals, plot = F)[1]
    selected.eigs <- which(abs(eig$values) >= vals[d])
    eig <- eig$vectors[,selected.eigs]
  }else {
    eig <- eigs(as(A, "dgeMatrix"), k = d)$vectors 
  }
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

distance.Rs <- function(Rs) {
  D <- sapply(Rs, function(A) 
    sapply(Rs, function(B) norm(A-B, "F")))
  return(D)
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
