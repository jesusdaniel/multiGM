library(igraph)
library(iGraphMatch)
library(Matrix)
library(pracma)
library(tidyverse)
library(viridis)
library(hrbrthemes)

BG_Generator <- function(num=5, n=100, d=10){
  set.seed(2022)
  matlist = vector(mode = "list", length = num)
  adjlist = vector(mode = "list", length = num)
  for (j in 1:num) {
    matlist[[j]] = matrix(runif(n*n), nrow = n)
    graphj = sample_correlated_ieg_pair(n, matlist[[j]], matrix(1, nrow = n, ncol = n))
    adjlist[[j]] = as_adjacency_matrix(graphj[[1]], type = "both", names = F, sparse = F)
  }
  
  res = mase(adjlist, d=d)
  
  cosielist = vector(mode = "list", length = num)
  truematlist = vector(mode = "list", length = num)
  
  for (k in 1:num) {
    matk = res$V %*% res$R[[k]] %*% t(res$V)
    matk[matk<0]=0
    matk[matk>1]=1
    cosielist[[k]] = matk
    graphk = sample_correlated_ieg_pair(n, cosielist[[k]], matrix(1, nrow = n, ncol = n))
    truematlist[[k]] = as_adjacency_matrix(graphk[[1]], type = "both", names = F, sparse = F)
  }
  return(list(matlist, adjlist, cosielist, truematlist))
}

G_Generator <- function(G, p){
  # Generate G_m's
  # G is the adjacency matrix of the imput background graph
  # p is the probability of flipping an edge
  d <- nrow(G)
  Gm = matrix(nrow=d, ncol=d)
  for (i in 1:d^2){
    B = rbinom(1,1,p)
    Gm[i] = ifelse(B, abs(G[i]-1), G[i])
  }
  diag(Gm) <- 0
  return(Gm)
}


f1 <- function(n=100, m=10, matlist, pflip=0.1){
  set.seed(2022)
  n_bg = length(matlist)
  exist_graphs <- vector(mode = "list", length = n_bg + (m-1))
  new_graphs <- vector(mode = "list", length = n_bg)
  #avg_graph_group <- vector(mode = "list", length = n_bg)
  for (j in 1:n_bg) {
    exist_graphs[[j]] = G_Generator(matlist[[j]],pflip)
    new_graphs[[j]] = G_Generator(matlist[[j]],pflip)
  }
  for (j in (n_bg+1):length(exist_graphs)) {
    exist_graphs[[j]] = G_Generator(matlist[[1]],pflip)
  }
  
  A_Full = 0
  for (k in 1:length(exist_graphs)){
    A_Full = A_Full + exist_graphs[[k]]
  }
  A_Full = A_Full/length(exist_graphs)
  
  objvec = numeric(n_bg)
  errvec = numeric(n_bg)
  
  for (k in 1:n_bg) {
    gmres <- graph_match_FW(A_Full, new_graphs[[k]], start = "bari",
                            seeds = cbind(1:5,1:5))
    gm.P <- gmres$P
    objvec[k] = norm(A_Full - gm.P%*% new_graphs[[k]] %*%t(gm.P), type = "F")
    errvec[k] = 1-sum(diag(gm.P))/n
  }
  return(list(objvec, errvec))
}

simulate_list = BG_Generator(num=10)
result_list = f1(matlist = simulate_list[[4]])
# m>=4 okay!!

testfunc <- function(n=100, m1=10, m2=5, matlist, pflip=0.1){
  set.seed(2022)
  n_bg = length(matlist)
  exist_graphs <- vector(mode = "list", length = (m1+m2))
  new_graphs <- vector(mode = "list", length = 2)
  for (j in 1:(m1+m2)) {
    if (j<=m1){
      exist_graphs[[j]] = G_Generator(matlist[[1]],pflip)
    }else{
      exist_graphs[[j]] = G_Generator(matlist[[2]],pflip)
    }
  }
  new_graphs[[1]] = G_Generator(matlist[[1]],pflip)
  new_graphs[[2]] = G_Generator(matlist[[2]],pflip)
  
  A_Full = 0
  for (k in 1:length(exist_graphs)){
    A_Full = A_Full + exist_graphs[[k]]
  }
  A_Full = A_Full/length(exist_graphs)
  
  objvec = numeric(2)
  errvec = numeric(2)
  
  for (k in 1:2) {
    gmres <- graph_match_FW(A_Full, new_graphs[[k]], start = "bari",
                            seeds = cbind(1:5,1:5))
    gm.P <- gmres$P
    objvec[k] = norm(A_Full - gm.P%*% new_graphs[[k]] %*%t(gm.P), type = "F")
    errvec[k] = 1-sum(diag(gm.P))/n
  }
  return(list(objvec, errvec))
}


f2 <- function(n=100, m1=10, m2=5, matlist, pflip=0.1){
  set.seed(2022)
  n_bg = length(matlist)
  n_mat=n_bg-1
  exist_graphs <- vector(mode = "list", length = (m1+2*m2))
  for (k in 1:m1) {
    exist_graphs[[k]] = G_Generator(matlist[[1]],pflip)
  }
  new_graph <- G_Generator(matlist[[1]],pflip)
  outmat1 = matrix(nrow = n_mat, ncol = n_mat)
  outmat2 = matrix(nrow = n_mat, ncol = n_mat)
  for (ind1 in 1:n_mat) {
    for (ind2 in ind1:n_mat) {
      for (j in (m1+1):(m1+2*m2)) {
        if (j<=m1+m2){
          exist_graphs[[j]] = G_Generator(matlist[[ind1+1]],pflip)
        }else{
          exist_graphs[[j]] = G_Generator(matlist[[ind2+1]],pflip)
        }
      }
      A_Full = 0
      for (k in 1:length(exist_graphs)){
        A_Full = A_Full + exist_graphs[[k]]
      }
      A_Full = A_Full/length(exist_graphs)
      gmres <- graph_match_FW(A_Full, new_graph, start = "bari",
                                seeds = cbind(1:5,1:5))
      gm.P <- gmres$P
      outmat1[ind1,ind2] = norm(A_Full - gm.P%*% new_graph %*%t(gm.P), type = "F")
      outmat2[ind1,ind2] = 1-sum(diag(gm.P))/n
      }
  }
  return(list(outmat1, outmat2))
}

result_list_2 <- f2(matlist = simulate_list[[4]], m2=6)
resmat1 <- as.matrix(forceSymmetric(result_list_2[[1]]))
resmat2 <- as.matrix(forceSymmetric(result_list_2[[2]]))


x <- 1:9
y <- 1:9
data <- expand.grid(X=x, Y=y)
data$Obj <- as.vector(resmat1[9:1,])
data$Err <- as.vector(resmat2[9:1,])


ggplot(data, aes(X, Y, fill= Obj)) + 
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black")+
  theme_ipsum()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none", 
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggplot(data, aes(X, Y, fill= Err)) + 
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black")+
  theme_ipsum()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none", 
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
