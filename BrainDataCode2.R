library(Matrix)
library(igraph)
library(iGraphMatch)
library(viridis)
library(ggplot2)
library(hrbrthemes)

Make.mat <- function(filename){
  # From a ssv file to get a adjacency matrix
  folderpath = "~/BrainData/"
  ext = ".ssv"
  filepath = paste(folderpath, filename, ext, sep = "")
  obsdata <- read.table(filepath, quote="\"", comment.char="")
  Edgelist <- as.vector(t(cbind(obsdata[1],obsdata[2])))
  Mat = as_adjacency_matrix(make_graph(Edgelist, n = 70, directed = F))
  return(Mat)
}
patientID = paste("00", 25427:25456, sep = "")
Full.list <- function(){
  #Run this function to get a full adjacency matrix list
  filenamelist <- vector(mode = "list", length = 300)
  for (i in 1:300){
    PatID = patientID[[((i-1)%/%10)+1]]
    scanID = ifelse((i%%10)==0, 10, i%%10)
    name = paste("sub-", PatID, "_ses-", scanID, "_dwi_desikan", sep = "")
    filenamelist[[i]] = name
  }
  return(lapply(filenamelist, Make.mat))
}

Adj_Mat_list_full = Full.list()


# This section performs the following 2 tasks simultinuously:
# 1. Average all existing data, for a given new graph of the patient, 
#    match it to the averaged graph and give the error
# 2. Cluster all existing data, for a given new graph of the patient, 
#    match it to all clusters, find the cluster it belongs to, and, 
#    give weather the classification is correct, and the error within the 
#    class it is matched to.

DoSection <- function(n, seed=2021){
  # Perform the described task, n is the number of patient using. 
  # For the function to make sense, 1<n<31
  set.seed(seed)
  indexUse = sample(1:30, n)
  patientUse = sort(indexUse)
  UsedpatientID = patientID[patientUse]
  print("PatientID used:")
  print(UsedpatientID)
  print("-------------------------------------------------------------")
  matindex = c()
  for (j in 1:n) {
    matindex = c(matindex, (((patientUse[j]-1)*10)+1):(patientUse[j]*10))
  }
  Adj_Mat_list = Adj_Mat_list_full[matindex]
  holdindex = sample(1:10, n, replace = T)+10*(0:(n-1))
  hold_mat = Adj_Mat_list[holdindex]
  Ave_mat = Adj_Mat_list[-holdindex]
  Ave = 0
  for (k in 1:(9*n)){Ave = Ave + Ave_mat[[k]]}
  Ave = Ave/(9*n)
  
  
  Ave_mat_list = vector(mode = "list", length = n)
  for (l in 1:n){
    Ave_9 = Ave_mat[(1:9)+(l-1)*9]
    Ave = 0
    for (k in 1:9){Ave = Ave + Ave_9[[k]]}
    Ave = Ave/9
    Ave_mat_list[[l]] = Ave
  }
  
  # Truthvalue <- vector(mode="character", length=n)
  # Error2 <- numeric(n)
  # Objective_Group <- numeric(n)
  # Individual_Error <- numeric(n)
  # Individual_Error_2 <- numeric(n)
  # Individual_Error_3 <- numeric(n)
  # Individual_Error_4 <- numeric(n)
  HeatMapMat_Obj <- matrix(nrow = n, ncol = n)
  HeatMapMat_Err <- matrix(nrow = n, ncol = n)
  for (a in 1:n) {
    now = hold_mat[[a]]
    #rel_error = numeric(n)
    for (j in 1:n){
      res <- gm(Ave_mat_list[[j]], now, start = "bari", seeds = cbind(1:5,1:5))
      P <- res$soft
      M <- now- P%*%Ave_mat_list[[j]]%*%t(P)
      HeatMapMat_Obj[a,j] = (norm(M, type = "F"))
      HeatMapMat_Err[a,j] = sum(diag(P))/70
    }
    # matchclass <- which(rel_error==min(rel_error))
    # Truthvalue[a] = ifelse(matchclass==a, "Success", "Fail")
    # res.final <- gm(Ave_mat_list[[matchclass]], now, start = "bari", seeds = cbind(1:5,1:5))
    # P.final <- res.final$soft
    # M_now <- now - P.final%*%Ave_mat_list[[matchclass]]%*%t(P.final)
    # Objective_Group[a] <- (norm(M_now, type = "F"))
    # Error2[a] = 1-sum(diag(P.final))/70
  }
  return(list(HeatMapMat_Obj,HeatMapMat_Err))
}

result_list_heatmap <- DoSection(15)

mat_1 <- result_list_heatmap[[1]][15:1, ]
mat_2 <- result_list_heatmap[[2]][15:1, ]


x <- 1:15
y <- 1:15
data_new <- expand.grid(X=x, Y=y)
data_new$Obj <- as.vector(mat_1)
data_new$Err <- as.vector(mat_2)


ggplot(data_new, aes(X, Y, fill= Obj)) + 
  geom_tile() +
  scale_fill_viridis() +
  theme_ipsum()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none", 
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


ggplot(data_new, aes(X, Y, fill= Obj)) + 
  geom_tile() +
  scale_fill_gradient(low = "black", high = "white")+
  theme_ipsum()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none", 
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


ggplot(data_new, aes(X, Y, fill= Err)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  theme_ipsum()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none", 
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


ggplot(data_new, aes(X, Y, fill= Err)) + 
  geom_tile() +
  scale_fill_gradient(low = "black", high = "white")+
  theme_ipsum()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none", 
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
