## Multiple graph matching single background w/ ER graphs

## library installation and loading, data acquisition
library(igraph)
library(iGraphMatch) # Must be V1, not V2
library(Matrix)


## Graph generators
G_Background_Generator <- function(n, Bernp=0.5){
  # Generate a background ER model.
  # n is number of nodes, Bernp is edge probability
  G1 <- erdos.renyi.game(n, Bernp, type = c("gnp"))
  G1adj <- as_adjacency_matrix(G1, type = "both", names = F, sparse = F)
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


## G is background. 
## All G_1 to G_m are correlated with same correlation.

Doindex_R <- function(index, n, m, Bernp=0.5){
  # Apply clustered matching
  # Switch Bernp to change the propability parameter for ER graph G(n,p)
  G = G_Background_Generator(n, Bernp)
  Graphs <- list()
  for (j in 1:m){
    Graphs = append(Graphs, list(G_Generator(G,Pflip[index])))
  }
  A = 0
  for (k in 2:m){
    A = A+Graphs[[k]]
  }
  A = A/(m-1) # A is the average of in-sample networks
  
  res <- graph_match_FW(A, Graphs[[1]], start = "bari", seeds = cbind(1:5,1:5))
  P <- res$P
  M <- Graphs[[1]] - P%*%A%*%t(P)
  Error.rate = (n-sum(diag(P)))/n
  return(c(norm(M, type = "F")),Error.rate)
}

Doindex_50_10_R <- function(index){
  Doindex_R(index, 50, 10)
}

Doindex_50_100_R <- function(index){
  Doindex_R(index, 50, 100)
}

Doindex_50_1000_R <- function(index){
  Doindex_R(index, 50, 1000)
}

Doindex_100_10_R <- function(index){
  Doindex_R(index, 100, 10)
}

Doindex_100_100_R <- function(index){
  Doindex_R(index, 100, 100)
}

Doindex_100_1000_R <- function(index){
  Doindex_R(index, 100, 1000)
}

Pflip = seq(0,0.5,by=0.025)

set.seed(106)

length_num = length(Pflip)
Res50.10 <- sapply(1:length_num, Doindex_50_10_R)
Loss50.10 = Res50.10[1,]
Error50.10 = Res50.10[2,]

Res50.100 <- sapply(1:length_num, Doindex_50_100_R)
Loss50.100 = Res50.100[1,]
Error50.100 = Res50.100[2,]

Res50.1000 <- sapply(1:length_num, Doindex_50_1000_R)
Loss50.1000 = Res50.1000[1,]
Error50.1000 = Res50.1000[2,]

Res100.10 <- sapply(1:length_num, Doindex_100_10_R)
Loss100.10 = Res100.10[1,]
Error100.10 = Res100.10[2,]

Res100.100 <- sapply(1:length_num, Doindex_100_100_R)
Loss100.100 = Res100.100[1,]
Error100.100 = Res100.100[2,]

Res100.1000 <- sapply(1:length_num, Doindex_100_1000_R)
Loss100.1000 = Res100.1000[1,]
Error100.1000 = Res100.1000[2,]

plot(Pflip, Loss50.10, main = "50 nodes, 10 graphs, Bernp=0.5",
     xlab = "Prob. flip edge of background", ylab = "Loss function")
plot(Pflip, Loss50.100, main = "50 nodes, 100 graphs, Bernp=0.5",
     xlab = "Prob. flip edge of background", ylab = "Loss function")
plot(Pflip, Loss50.1000, main = "50 nodes, 1000 graphs, Bernp=0.5",
     xlab = "Prob. flip edge of background", ylab = "Loss function")
plot(Pflip, Loss100.10, main = "100 nodes, 10 graphs, Bernp=0.5",
     xlab = "Prob. flip edge of background", ylab = "Loss function")
plot(Pflip, Loss100.100, main = "100 nodes, 100 graphs, Bernp=0.5",
     xlab = "Prob. flip edge of background", ylab = "Loss function")
plot(Pflip, Loss100.1000, main = "100 nodes, 1000 graphs, Bernp=0.5",
     xlab = "Prob. flip edge of background", ylab = "Loss function")