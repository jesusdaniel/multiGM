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
  
  Error1 = numeric(n)
  obj1 = numeric(n)
  for (l in 1:n){
    res <- gm(Ave, hold_mat[[l]], start = "bari", seeds = cbind(1:5,1:5))
    P <- res$soft
    Error1[l] = 1-sum(diag(P))/70
    obj1[l] = norm(Ave - P%*%hold_mat[[l]]%*%t(P),type = "F")
  }
  print("Task 1 summary")
  res1 = t(rbind(UsedpatientID, obj1, Error1))
  colnames(res1) <- c("Patient ID", "Objective Value", "Error")
  print(res1)
  print("-------------------------------------------------------------")
  
  Ave_mat_list = vector(mode = "list", length = n)
  for (l in 1:n){
    Ave_9 = Ave_mat[(1:9)+(l-1)*9]
    Ave = 0
    for (k in 1:9){Ave = Ave + Ave_9[[k]]}
    Ave = Ave/9
    Ave_mat_list[[l]] = Ave
    }
  
  Truthvalue <- vector(mode="character", length=n)
  Error2 <- numeric(n)
  Objective_Group <- numeric(n)
  Individual_Error <- numeric(n)
  Individual_Error_2 <- numeric(n)
  Individual_Error_3 <- numeric(n)
  Individual_Error_4 <- numeric(n)
  for (a in 1:n) {
    now = hold_mat[[a]]
    rel_error = numeric(n)
    for (j in 1:n){
      res <- gm(Ave_mat_list[[j]], now, start = "bari", seeds = cbind(1:5,1:5))
      P <- res$soft
      M <- now- P%*%Ave_mat_list[[j]]%*%t(P)
      rel_error[j] = (norm(M, type = "F"))
    }
    matchclass <- which(rel_error==min(rel_error))
    Truthvalue[a] = ifelse(matchclass==a, "Success", "Fail")
    res.final <- gm(Ave_mat_list[[matchclass]], now, start = "bari", seeds = cbind(1:5,1:5))
    P.final <- res.final$soft
    M_now <- now - P.final%*%Ave_mat_list[[matchclass]]%*%t(P.final)
    #This.Group = Ave_mat[(1:9)+(matchclass-1)*9]
    #errorlist <- numeric(9)
    #errorlist_2 = numeric(9)
    #for (k in 1:9) {
    #  res.indi <- graph_match_FW(This.Group[[k]], now, start = "bari", seeds = cbind(1:5,1:5))
    #  P.indi <- res.indi$P
    #  M.indi <- now - P.indi%*%This.Group[[k]]%*%t(P.indi)
    #  errorlist[k] <- (norm(M.indi, type = "F"))
    #  errorlist_2[k] <- 1-sum(diag(P.indi))/70
    #}
    #Individual_Error[a] <- min(errorlist)
    #Individual_Error_2[a] <- errorlist_2[which(min(errorlist) == errorlist)]
    #Individual_Error_3[a] <- min(errorlist_2)
    #Individual_Error_4[a] <- errorlist[which(min(errorlist_2) == errorlist_2)]
    Objective_Group[a] <- (norm(M_now, type = "F"))
    Error2[a] = 1-sum(diag(P.final))/70
  }
  result = rbind(UsedpatientID, Truthvalue, Error2, Objective_Group)
  result = t(result)
  colnames(result) <- c("Patient ID", "Matching to true class", "Error of the match", 
                        "Objective function value")
  print("Task 2 summary")
  print(result)
  
  
  print("-------------------------------------------------------------")
  
  matrixinfo.obj = matrix(nrow = n, ncol = 9*n)
  matrixinfo.err = matrix(nrow = n, ncol = 9*n)
  
  for (ind1 in 1:n) {
    for (ind2 in 1:(9*n)) {
      #print(c(ind1,ind2))
      res.out <- gm(hold_mat[[ind1]], Ave_mat[[ind2]], 
                                start = "bari", seeds = cbind(1:5,1:5))
      P.this <- res.out$soft
      M <- hold_mat[[ind1]]- P.this%*%Ave_mat[[ind2]]%*%t(P.this)
      matrixinfo.obj[ind1, ind2] = (norm(M, type = "F"))
      matrixinfo.err[ind1, ind2] = 1-sum(diag(P.this))/70
    }
    
  }
  
  
  return(list(obj1, Error1, Error2, Objective_Group, matrixinfo.obj, matrixinfo.err))
  #print("-------------------------------------------")
  #result2 = rbind(UsedpatientID, Truthvalue, Error2, Individual_Error_3, 
  #               Objective_Group, Individual_Error_4)
  #result2 = t(result2)
  #colnames(result2) <- c("Patient ID", "Matching to true class", "Error of the match", "Error of indi match",
  #                      "Objective function from average", "Objective function from individual")
  #print("Task 2 summary, table 2: use error rate")
  #print(result2)
  #print("-------------------------------------------")
}

result_list <- DoSection(15)

mat_test_1 = cbind(result_list[[1]], result_list[[4]], result_list[[5]])
mat_test_2 = cbind(result_list[[2]], result_list[[3]], result_list[[6]])

mat_1 <- mat_test_1[15:1, ]
mat_2 <- mat_test_2[15:1, ]
vec1 <- list(mat_1[,1], mat_1[,2])
vec2 <- list(mat_2[,1], mat_2[,2])
repmat1 = matrix(data = vec1[[1]], nrow = 15, ncol = 10)
repmat2 = matrix(data = vec1[[2]], nrow = 15, ncol = 10)
repmat3 = matrix(data = vec2[[1]], nrow = 15, ncol = 10)
repmat4 = matrix(data = vec2[[2]], nrow = 15, ncol = 10)

mat1 <- cbind(repmat1, repmat2, mat_1[,3:137])
mat2 <- cbind(repmat3, repmat4, mat_2[,3:137])


x <- 1:15
y <- 1:155
data <- expand.grid(X=x, Y=y)
data$Obj <- as.vector(mat1)
data$Err <- as.vector(mat2)


ggplot(data, aes(X, Y, fill= Obj)) + 
  geom_tile() +
  scale_fill_viridis() +
  theme_ipsum()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none", 
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


ggplot(data, aes(X, Y)) + 
  geom_tile(aes(fill= Obj)) +
  scale_fill_gradient(low = "white", high = "black")+
  guides(fill = guide_colorbar(title = 'value'))+
  theme_ipsum()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


ggplot(data, aes(X, Y, fill= Err)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
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
  guides(fill = guide_colorbar(title = 'value'))+
  theme_ipsum()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
