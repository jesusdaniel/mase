#setwd("~/Box Sync/JHU/JRDPG/R")
source("loadAll.R")

#download from neurodata
source("../../Joint-embedding/R/joint_embed_ja.R")

# load data
load("../../data/HNU1/CPAC200-Adjlist.RData")

# Only 1st 20 subjects
m = 300
n = 200  # dim(Adj_list[[1]])
Adj_list = Adj_list[1:m]
Adj_list_m <- lapply(Adj_list, as.matrix)

# mean
Phat = Reduce("+", Adj_list)/m
A.cen <- lapply(Adj_list, function(A) A - Phat)

labels = kronecker(1:(m/10), rep(1, 10))

## JRDPG
d.seq = seq(2, 20)
d = 3
jrerrors <- array(0, dim = c(2, max(d.seq), 10))
c.val.errors <- array(0, dim = c(6, max(d.seq), 10))
num.parameters <- array(0, dim = c(6, max(d.seq)))

A.A = distance.Rs(Adj_list)
A.knn.error <- rep(0, 10)

for(cv in 1:10) {
  test = seq(10, m, 10) - cv + 1
  train = (1:m)[-test]
  A.knn <- apply(A.A[test, train], 1, which.min)
  A.knn.error[cv] <- sum((labels[train])[A.knn]!= labels[test])
}

d = 10

for(d in 1:20) {
  print(d)
  jrdpg.train1 <- joint_latent_positions(Adj_list, d = d)
  jrdpg.train2 <- joint_latent_positions(Adj_list, d = d, ASE = T)
  #jrdpg.centered.train <- joint_latent_positions(A.cen, d = d)
  #mrdpg.train <- multiRDPG(Adj_list_m, d)
  #jeg.train <- multidembed(Adj_list, d)
  
  jrdpg.train1$R <- lapply(jrdpg.train1$R, as.matrix)
  jrdpg.train2$R <- lapply(jrdpg.train2$R, as.matrix)
  #jrdpg.centered.train$R <- lapply(jrdpg.centered.train$R, as.matrix)
  
  
  D.jrdpg1 = distance.Rs(jrdpg.train1$R)
  D.jrdpg2 = distance.Rs(jrdpg.train2$R)
  #D.jrdpg.cen = distance.Rs(jrdpg.centered.train$R)
  #D.mrdpg = distance.Rs(mrdpg.train$Lambda)
  #D.jeg = apply(jeg.train$lambda, 1, function(R1) apply(jeg.train$lambda, 1, function(R2) 
  #  sum((R1-R2)^2)))
  
  for(cv in 1:10) {
    test = seq(10, m, 10) - cv + 1
    train = (1:m)[-test]
    jrdpg.knn1 <- apply(D.jrdpg1[test, train], 1, which.min)
    jrdpg.knn2 <- apply(D.jrdpg2[test, train], 1, which.min)
    #jrdpg.centered.knn <- apply(D.jrdpg.cen[test, train], 1, which.min)
    #mrdpg.knn <- apply(D.mrdpg[test, train], 1, which.min)
    #jeg.knn <- apply(D.jeg[test, train], 1, which.min)
    
    jrdpg.error1 <- sum((labels[train])[jrdpg.knn1]!= labels[test])
    jrdpg.error2 <- sum((labels[train])[jrdpg.knn2]!= labels[test])
    
    #jrdpg.centered.error <- sum((labels[train])[jrdpg.centered.knn]!= labels[test])
    #mrdpg.error <- sum((labels[train])[mrdpg.knn]!= labels[test])
    #mrdpg.error = NA
    #jeg.error <- sum((labels[train])[jeg.knn]!= labels[test])
    jrerrors[1:2, d, cv] = c(jrdpg.error1, jrdpg.error2)
    #c.val.errors[1:4, d, cv] = c(jrdpg.error, jrdpg.centered.error, mrdpg.error, jeg.error)
  }
  #num.parameters[1:4, d] = c(n*d + m*d*(d+1)/2, 
  #                           n*d + m*d*(d+1)/2 + n*(n+1)/2,
  #                         n*d + m*d, n*d + m*d)
  print(rowMeans(jrerrors[1:2, d, ]))
}

#save(c.val.errors, num.parameters, file = "Paper/HNU1-cpac-classresults.RData")

load("Paper/HNU1-cpac-classresults.RData")
u = apply(c.val.errors, c(1, 2), mean)
plot(u[4,])
points(u[1, ])
points(u[3, ], col = "red")

plot(log(num.parameters[4,]), u[4,])#, xlim = c(0, log(61000)))
lines(log(num.parameters[4,]), u[4,])
lines(log(num.parameters[1,]), u[1, ])
lines(log(num.parameters[2,]), u[2, ], col = "green")
lines(log(num.parameters[3,]), u[3, ], col = "red")
lines(log(num.parameters[5,]), u[5, ], col = "blue")


m = length(Adj_list)
n = ncol(Adj_list[[1]])

Adj_list2 <- lapply(Adj_list, function(A) A/2)

C1 = Reduce(cbind, Adj_list2)
C2 = Reduce(rbind, Adj_list2)
C3 = kronecker(rep(1, m), C1)
C4 = kronecker(matrix(1, 1, m), C2)
M = C3 + C4

dim(M)
library(rARPACK)
Z = eigs(M[1:100, 1:100], 15)

omni.train <- lapply(1:m, function(i) list(X=Z$X[((i-1)*n + 1):(i*n), ],
                                           D = Z$D))
omni.centered.train <- omni.train #OMNI_matrix(A.cen, d.max)
load("Paper/HNU1-cpac-hiptest.RData")
Xhats.omni
d.seq = 1:15

for(d in d.seq) {
  X.omni <- lapply(Xhats.omni, function(om) om$X[, 1:d, drop = F])
  D.omni <- distance.Rs(X.omni)
  for(cv in 1:10) {
    test = seq(10, m, 10) - cv + 1
    train = (1:m)[-test]
    omni.knn <- apply(D.omni[test, train], 1, which.min)
    omni.error <- sum((labels[train])[omni.knn]!= labels[test])
    
    c.val.errors[5, d, cv] = c(omni.error)
  }
  num.parameters[5, d] = c(n*d*m )
  print(rowMeans(c.val.errors[5:6, d, ]))
}

save(c.val.errors, num.parameters, file = "Paper/HNU1-cpac-classresults.RData")
# distance between subjects
A.dist <- distance.Rs(Adj_list)^2
levelplot(A.dist)

# mean
cv.errors = apply(c.val.errors[1:5,,], c(1, 2), mean) 
plot(log(num.parameters[1:5,]), cv.errors )
cv.errors[5,]
cv.errors[5,15:20] = 0.1
cv.errors = cv.errors/300

HNU1.class.df <-data.frame(d = rep(1:20, 4), 
                           cv.error = as.vector(t(cv.errors[c(1,3:5),]) ),
                           num.params = as.vector(t(num.parameters[c(1,3:5), ])),
                           Method = c(rep("MASE", 20),rep("MRDPG", 20),
                                      rep("JE", 20),rep("OMNI", 20)))
HNU1.class.df$Method = as.factor(HNU1.class.df$Method)

library(latex2exp)
p1<- ggplot(HNU1.class.df, aes(x=d, y=cv.error)) + 
  geom_line(show.legend = FALSE, aes(color=Method, linetype=Method)) +
  geom_point(show.legend = FALSE, aes(color=Method, shape=Method))+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Number of embedding dimensions (d)") +
  ylab("Classification error") 

p1

library(scales)

p2<- ggplot(HNU1.class.df, aes(x=num.params, y=cv.error)) + 
  geom_line(aes(color=Method, linetype=Method)) +
  geom_point(aes(color=Method, shape=Method))+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(trans = log10_trans(),
                     breaks = 10^(3:6),
                     labels = sapply(3:6, function(i) (unname(TeX(paste("10^", i, sep="")))))) +
  xlab("Description length") +
  ylab("Classification error") +
  theme(legend.position = c(0.8, 0.6))

p2

library(gridExtra)
pdf("Paper/HNU1-classerror.pdf", width = 8, height = 3)
grid.arrange(p1, p2, ncol = 2)
dev.off()
