# Bootstrap code
test_equality_jrdpg <- function(V, R1, R2, B = 1000, P0 = NULL) {
  null_1 <- bootstrap_jrdpg2(V, R1, P0 = P0, B = B)
  null_2 <- bootstrap_jrdpg2(V, R2, P0 = P0, B = B)
  test <- norm(R1 - R2, type = "F")
  pval <- max(p.values(test, null_1), p.values(test, null_2))
  return(list(p.value = pval, null.dist.1 = null_1, null.dist.2 = null_2, 
              test.statistic = test))
}

frobenius <- function(A, B) norm(A-B, "F")

test_equality_jrdpg.ase <- function(A1, A2, d, B = 1000, P0 = NULL) {
  jr1 = joint_latent_positions(list(A1), d)
  jr2 = joint_latent_positions(list(A2), d)
  null_1 <- bootstrap_jrdpg1(jr1$V, jr1$R[[1]], P0 = P0, B = B)
  null_2 <- bootstrap_jrdpg1(jr2$V, jr2$R[[1]], P0 = P0, B = B)
  jr = joint_latent_positions(list(A1, A2), d)
  test <- norm(jr$R[[1]] - jr$R[[2]], type = "F")
  pval <- max(p.values(test, null_1), p.values(test, null_2))
  return(list(p.value = pval, null.dist.1 = null_1, null.dist.2 = null_2, 
              test.statistic = test))
}



test_equality_jrdpg.montecarlo <- function(A1, A2, V, R1, R2, B = 1000) {
  d <- ncol(V)
  null_1 <- bootstrap_jrdpg(V, R1, B = B)
  null_2 <- bootstrap_jrdpg(V, R2, B = B)
  jrdpg <- joint_latent_positions(list(A1, A2), d)
  test <- norm(jrdpg$R[[1]] - jrdpg$R[[2]], type = "F")
  pval <- max(p.values(test, null_1), p.values(test, null_2))
  return(list(p.value = pval, null.dist.1 = null_1, null.dist.2 = null_2, 
              test.statistic = test))
}

test_equality_omni <- function(A1, A2, d, B = 1000) {
  X1.gase = g.ase(A1, d)
  X2.gase = g.ase(A2, d)
  null_1 <- bootstrap_omni(X1.gase, B = B)
  null_2 <- bootstrap_omni(X2.gase, B = B)
  omni <- OMNI_embedding(A1, A2, d)
  test <- norm(omni$X1 - omni$X2, type = "F")
  pval <- max(p.values(test, null_1), p.values(test, null_2))
  return(list(p.value = pval, null.dist.1 = null_1, null.dist.2 = null_2, 
              test.statistic = test))
}

test_equality_omni.montecarlo <- function(A1, A2, X1.gase, X2.gase, B = 1000) {
  d <- ncol(X1.gase$X)
  null_mc1 <- bootstrap_omni(X1.gase, B = B)
  null_mc2 <- bootstrap_omni(X2.gase, B = B)
  omni <- OMNI_embedding(A1, A2, d)
  test <- norm(omni$X1 - omni$X2, type = "F")
  pval <- max(p.values(test, null_mc1), p.values(test, null_mc2))
  return(list(p.value = pval, null.dist.1 = null_mc1, null.dist.2 = null_mc2, 
              test.statistic = test))
}

bootstrap_jrdpg_from_P <- function(P, d, K = d, ASE = F, P0 = NULL, B = 1000) {
  n = nrow(P)
  
  if(is.null(P0)) P0 = Matrix(0, n, n)
  P1 <- P + P0
  #P1 = P1 * (P1>0)
  #P1 = P1*(P1 <1) + 1* (P1>1)
  Ab1 <- 0*P1
  Ab2 <- 0*P1
  bootstrap <- sapply(1:B, function(i) {
    Ab1[upper.tri(Ab1)] <- 1*(runif(sum(upper.tri(P1))) < P1[upper.tri(P1)])
    Ab2[upper.tri(Ab1)] <- 1*(runif(sum(upper.tri(P1))) < P1[upper.tri(P1)])
    Ab1 = (Ab1 + t(Ab1))
    Ab2 = (Ab2 + t(Ab2))
    jrdpg <- joint_latent_positions(list(Ab1-P0, Ab2-P0), d = d, K = K, ASE = ASE)
    return(norm(jrdpg$R[[1]]- jrdpg$R[[2]], "F"))
  })
  return(bootstrap)
}




bootstrap_jrdpg1 <- function(V, R, P0 = NULL, B = 1000) {
  n = nrow(V)
  d = ncol(V)
  
  if(is.null(P0)) P0 = Matrix(0, n, n)
  P1 <- V %*% R %*% t(V) + P0
  #P1 = P1 * (P1>0)
  #P1 = P1*(P1 <1) + 1* (P1>1)
  Ab1 <- 0*P1
  Ab2 <- 0*P1
  bootstrap <- sapply(1:B, function(i) {
    Ab1[upper.tri(Ab1)] <- 1*(runif(sum(upper.tri(P1))) < P1[upper.tri(P1)])
    Ab2[upper.tri(Ab1)] <- 1*(runif(sum(upper.tri(P1))) < P1[upper.tri(P1)])
    Ab1 = (Ab1 + t(Ab1))
    Ab2 = (Ab2 + t(Ab2))
    jrdpg <- joint_latent_positions(list(Ab1-P0, Ab2-P0), d = d, K = d, ASE = F)
    return(norm(jrdpg$R[[1]]- jrdpg$R[[2]], "F"))
  })
  return(bootstrap)
}

bootstrap_jrdpg2 <- function(V, R, P0 = NULL, B = 1000) {
  n = nrow(V)
  d = ncol(V)
  if(is.null(P0)) P0 = Matrix(0, n, n)
  P1 <- V %*% R %*% t(V) + P0
  #P1 = P1 * (P1>0)
  #P1 = P1*(P1 <1) + 1* (P1>1)
  Ab1 <- 0*P1
  Ab2 <- 0*P1
  bootstrap <- sapply(1:B, function(i) {
    Ab1[upper.tri(Ab1)] <- 1*(runif(sum(upper.tri(P1))) < P1[upper.tri(P1)])
    Ab2[upper.tri(Ab1)] <- 1*(runif(sum(upper.tri(P1))) < P1[upper.tri(P1)])
    Ab1 = (Ab1 + t(Ab1))
    Ab2 = (Ab2 + t(Ab2))
    Rb1 <- crossprod(V, crossprod(Ab1-P0, V))
    Rb2 <- crossprod(V, crossprod(Ab2-P0, V))
    return(norm(Rb1 - Rb2, "F"))
  })
  return(bootstrap)
}

bootstrap_jrdpg3 <- function(V, R, P0 = NULL, B = 1000) {
  n = nrow(V)
  if(is.null(P0)) P0 = Matrix(0, n, n)
  P1 <- V %*% R %*% t(V) + P0
  Ab1 <- 0*P1
  bootstrap <- sapply(1:B, function(i) {
    Ab1[upper.tri(Ab1)] <- 1*(runif(sum(upper.tri(P1))) < P1[upper.tri(P1)])
    Ab1 = (Ab1 + t(Ab1))
    Rb1 <- crossprod(V, crossprod(Ab1-P0, V))
    Rdiff <- Rb1 - R
    return(norm(Rdiff, "F"))
  })
  return(bootstrap)
}
p.values <- function(test.statistics, null.distribution) {
  return(sapply(test.statistics, function(test) 
    (sum(test <= null.distribution) + 0.5)/ length(null.distribution)))
}
