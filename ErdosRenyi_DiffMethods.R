#### Multiple graph matching Different methods w/ ER graphs

## Library installation and loading, data acquisition
library(igraph)
library(iGraphMatch) # Must be v1 not v2
library(Matrix)
library(tidyverse)


## Usages
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

###############################################################
## G, B are two different background. 
## All G_1 to G_m are generated from G
## All B_1 to B_m are generated from B
## G_1 and B_1 are individually matched to M_bar: full avg and G_bar, B_bar

pflip = seq(0.1, 0.5, by=0.1)

LoopFunc <- function(i, n=80, m1=201, m2=2001, Bernp=c(0.2, 0.4)){
  G = G_Background_Generator(n, Bernp[1])
  B = G_Background_Generator(n, Bernp[2])
  Graphs <- vector(mode = "list", length = (m1+m2))
  for (j in 1:m1){
    Graphs[[j]] = G_Generator(G,pflip[i])
  }
  
  for (j in ((m1+1):(m1+m2))){
    Graphs[[j]] = G_Generator(B,pflip[i])
  }
  
  A_Full = 0
  for (k in 2:(m1+m2)){
    if (k != m1+1){A_Full = A_Full + Graphs[[k]]}
  }
  A_Full = A_Full/(m1+m2-2)

  A1 = 0
  for (k in 2:m1){
    A1 = A1 + Graphs[[k]]
  }
  A1 = A1/(m1-1)
  
  A2 = 0
  for (k in ((m1+2):(m1+m2))){
    A2 = A2 + Graphs[[k]]
  }
  A2 = A2/(m2-1)
  
  res_f_1 <- graph_match_FW(A_Full, Graphs[[1]], start = "bari", seeds = cbind(1:5,1:5))
  res_f_2 <- graph_match_FW(A_Full, Graphs[[m1+1]], start = "bari", seeds = cbind(1:5,1:5))
  res_1_1 <- graph_match_FW(A1, Graphs[[1]], start = "bari", seeds = cbind(1:5,1:5))
  res_1_2 <- graph_match_FW(A1, Graphs[[m1+1]], start = "bari", seeds = cbind(1:5,1:5))
  res_2_1 <- graph_match_FW(A2, Graphs[[1]], start = "bari", seeds = cbind(1:5,1:5))
  res_2_2 <- graph_match_FW(A2, Graphs[[m1+1]], start = "bari", seeds = cbind(1:5,1:5))
  
  P_f_1 <- sum(diag(res_f_1$P))/n
  P_f_2 <- sum(diag(res_f_2$P))/n
  P_1_1 <- sum(diag(res_1_1$P))/n
  P_2_2 <- sum(diag(res_2_2$P))/n
  P_1_2 <- sum(diag(res_1_2$P))/n
  P_2_1 <- sum(diag(res_2_1$P))/n
  
  P_f_G <- res_f_1$P
  P_f_B <- res_f_2$P
  P_G <- res_1_1$P
  P_B <- res_2_2$P
  P_B_G <- res_2_1$P
  P_G_B <- res_1_2$P
  
  obj_f_G <- norm(A_Full - P_f_G%*% Graphs[[1]] %*%t(P_f_G), type = "F")
  obj_f_B <- norm(A_Full - P_f_B%*% Graphs[[m1+1]] %*%t(P_f_B), type = "F")
  obj_G <- norm(A1 - P_G%*% Graphs[[1]] %*%t(P_G), type = "F")
  obj_B <- norm(A2 - P_B%*% Graphs[[m1+1]] %*%t(P_B), type = "F")
  obj_B_G <- norm(A2 - P_B_G%*% Graphs[[1]] %*%t(P_B_G), type = "F")
  obj_G_B <- norm(A1 - P_G_B%*% Graphs[[m1+1]] %*%t(P_G_B), type = "F")
  return(c(P_f_1, P_f_2, P_1_1, P_2_2, P_1_2, P_2_1,
           obj_f_G, obj_G, obj_B_G, obj_f_B, obj_B, obj_G_B))
}


result_list_MC50 = vector(mode = "list", length = 50)
set.seed(2022)
for (indexj in 1:50){
  result_list_MC50[[indexj]] <- lapply(1:5, LoopFunc)
}
#result_list <- lapply(1:5, LoopFunc)

result_list_MC3 = vector(mode = "list", length = 3)
set.seed(2022)
for (indexj in 1:3){
  result_list_MC50[[indexj]] <- lapply(1:5, LoopFunc)
}


r <- lapply(result_list, `[[`, 1)
result <- matrix(unlist(r), ncol=5)
colnames(result) <- c("q=0.1", "q=0.2", "q=0.3", "q=0.4", "q=0.5")
rownames(result) <- c("Full & G group", "Full & B group", "G group", 
                      "B group", "B with G group", "G with B group")
print(result)

pl <- unlist(lapply(result_list, `[[`, 2))
obj1 <- pl[seq(1, length(pl), 6)]
obj2 <- pl[seq(2, length(pl), 6)]
obj3 <- pl[seq(3, length(pl), 6)]
obj4 <- pl[seq(4, length(pl), 6)]
obj5 <- pl[seq(5, length(pl), 6)]
obj6 <- pl[seq(6, length(pl), 6)]
plot(pflip, obj1, type = "p",pch=15, col = "blue", main = "Objective function plot",
     xlab = "q(edge flip prob)", ylab = "Objective value", ylim = c(20,50))
points(pflip, obj2, col="deepskyblue4", pch=15)
points(pflip, obj3, col="blueviolet", pch=15)

points(pflip, obj4, col="darkorange", pch=17)
points(pflip, obj5, col="coral4", pch=17)
points(pflip, obj6, col="red", pch=17)

legend("bottomright", legend=c("Full & G group", "G to G avg", "G to B avg",
                               "Full & B group", "B to B avg", "B to G avg"),
       col = c("blue","deepskyblue4", "blueviolet", 
               "darkorange", "coral4", "red"), 
       pch=c(15,15,15,17,17,17), merge=FALSE)


results = data.frame(q=pflip, A1_to_Full_avg = obj1, 
                     A1_to_S_Avg = obj2, 
                     A1_to_T_Avg = obj3,
                     A2_to_Full_avg = obj4,
                     A2_to_T_avg = obj5, 
                     A2_to_S_avg = obj6)

results_New <- results %>% pivot_longer(cols = A1_to_Full_avg:A2_to_S_avg, 
                                    names_to = "Methods", 
                                    values_to = "Value")
ggplot(results_New,  aes(x = q, y = Value, group = Methods))+
  geom_line(aes(linetype = Methods), size=0.8)+
  geom_point(aes(pch = Methods), size=2.5)+
  ggtitle("Objective function plot for different averaging methods")+
  xlab("Edge Flip Probability")+ylab("Objective Function Value")+
  labs(fill = "Matching Performed")+
  theme_bw(base_size=15)
