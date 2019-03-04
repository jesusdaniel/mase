### Comparing 4 classes
library(igraph)
library(graphclass) #github jesusdaniel/graphclass (only needed for the plots)
library(ggplot2)
library(igraph)
library(multiRDPG)
source("j-RDPG.R")
#download from neurodata
source("../../Joint-embedding/R/joint_embed_ja.R")


m = 40
n = 200
d=2
alpha = 0.5


# plot eigenvectors on 2 dimensions
cols <- c(rep("orange", n/2),rep("red", n/2))
xlims = c(-1, 1)
xlims = c(-0.2, 0.2)
ylims = c(-1, 1)
ylims = c(-0.2, 0.2)
plot_eig <- function(V) {
  df <- data.frame(x = V[,1], y = V[,2], col=cols)
  print(ggplot(df, aes(x=x, y=y)) + 
          #geom_hline(yintercept = 0)+
          #geom_vline(xintercept = 0)+
          geom_point(aes(col="black",fill=cols),shape=21,size=3) +
          xlim(xlims) + ylim(ylims) + xlab("") +
          ylab("") +
          scale_color_manual(values=c("black"))+
          scale_fill_manual(values=c("orange", "red"))+
          theme_bw() + theme(legend.position="none")) 
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #      panel.grid.minor = element_blank(), 
  #      axis.line = element_blank())
  #        element_line(colour = "black"))
}





Bs = lapply(1:m, function(i) {
  B0 = matrix(0.15,  d, d)
  b1 <- if(i%%4==1) {
    c(0.1, 0.0, 0.1)
  } else{ if(i%%4==2) {
    -c(0.1, 0, 0.1)
  } else{ if(i%%4==3) {
    c(0.1, 0.0, 0.0)
  }else {
    c(0.0, 0.0, 0.1)
  }}}
  B = matrix(c(b1[1:2], b1[2], b1[3]), nrow = 2)
  B0 + alpha * B
})
labels <- factor(rep(1:4, m/4))

Z <- matrix(0, n, 2)
Z[1:(n/2), 1] = Z[(n/2 + 1):n,2] = 1

Adj_list <- lapply(Bs, function(B) {
  get.adjacency(sample_sbm(n, B, block.sizes = rep(n/2, 2),
                           directed = FALSE))
})

graphclass:::plot_adjmatrix(Adj_list[[1]])
graphclass:::plot_adjmatrix(Adj_list[[2]])
graphclass:::plot_adjmatrix(Adj_list[[3]])
graphclass:::plot_adjmatrix(Adj_list[[4]])


# fitting the three methodds
jrdpg <- joint_latent_positions(Adj_list, d = 2, K = d, ASE = F)

Adj_list_m <- lapply(Adj_list, as.matrix)
mrdpg <- multiRDPG(Adj_list_m, d = 3)


jeg <- multidembed(Adj_list, d = 3, Innitialize = 1)

# 2d plots of the latent positions
plot_eig(jrdpg$V)
plot_eig(jeg$h[,1:2])
plot_eig(mrdpg$U[,1:2])


# get network coordinates
mrdpg_L = t(sapply(mrdpg$Lambda, diag))
jrdpg_L = t(sapply(jrdpg$R, function(R) R[upper.tri(R, diag=TRUE)]))


# 3d plots of the network coordinates
plot.3d((mrdpg_L), col = labels)
plot.3d((jrdpg_L), col = labels)
plot.3d((jeg$lambda), col = labels)


##### Measure accuracies
subspace_distance(jrdpg$V, Z)
subspace_distance(mrdpg$U[, 1:2], Z)
subspace_distance(jeg$h[, 1:2], Z)

##classification accuracy

# create test set
Adj_list_test <- lapply(Bs, function(B) {
  get.adjacency(sample_sbm(n, B, block.sizes = rep(n/2, 2),
                           directed = FALSE))
})

Bhat_test <- project_networks(Adj_list_test, jrdpg$V)
labels_test = labels


# use knn
Xtrain <- t(sapply(jrdpg$R, as.vector))
Xtest <- t(sapply(Bhat_test, as.vector))

library(class)
# classification accuracy
sum(knn(train = Xtrain, test = Xtest, cl = labels, k = 1) == labels) /  length(labels)


# mrdpg
# project mrdpg test set
Xtrain.mrdpg <- mrdpg_L
Xtest.mrdpg <- t(sapply(Adj_list_test, function(A) diag((t(mrdpg$U) %*% crossprod(A, mrdpg$U) ))))
sum(knn(train = Xtrain.mrdpg, test = Xtest.mrdpg, cl = labels, k = 1) == labels) /  length(labels)



# jeg
Xtrain.jeg <- jeg$lambda
Xtest.jeg <- multiembed_test(A = Adj_list_test, h = jeg$h)
sum(knn(train = Xtrain.jeg, test = Xtest.jeg, cl = labels, k = 1) == labels) /  length(labels)
