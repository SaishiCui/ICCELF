rm(list=ls())
gc()

load(file="C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_coarse_1.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_coarse_label_1.RData")
load(file="C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_fine_1.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_fine_label_1.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Training_data/PBMC1coarse_trainingdata.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Training_data/PBMC1fine_trainingdata.RData")



PBMC_coarse_1_no_bal_no_comp<- t(PBMC_coarse_1)[,colSums(t(PBMC_coarse_1)) != 0]
PBMC_fine_1_no_bal_no_comp<- t(PBMC_fine_1)[,colSums(t(PBMC_fine_1)) != 0]




load(file = "C:/Users/css22/Desktop/Thesis2/Models/xgboost_PBMC1coarse_no_bal_no_comp.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Models/xgboost_PBMC1fine_no_bal_no_comp.RData")




### XGboost model

library(xgboost)


xgboost_train_compare <- function(data1, cell.info){
  

  label = as.integer(factor(cell.info))-1
    

  
  dtrain <- xgb.DMatrix(data = as.matrix(data1), label =  label)    
  params <- list(booster = "gbtree", objective = "multi:softprob",   eval_metric="mlogloss", num_class= length(unique(cell.info)),
                 eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample= 0.5, colsample_bytree= 0.5)
  
  xgbcv <- xgb.cv( params = params, data = dtrain, nrounds = 200, nfold = 10, showsd = T, 
                                              stratified = T, print_every_n = 1, early_stop_round = 50, maximize = F)
  
  nrounds = which(xgbcv$evaluation_log[,4] ==  min(xgbcv$evaluation_log[,4]) )
  

  model <- xgb.train (params = params, data = dtrain, nrounds = nrounds, print_every_n = 1, maximize = F)

  
  return(model)
}



xgboost_PBMC1coarse_no_bal_no_comp = xgboost_train_compare(data1=PBMC_coarse_1_no_bal_no_comp, cell.info = PBMC_coarse_label_1)
xgboost_PBMC1fine_no_bal_no_comp = xgboost_train_compare(data1=PBMC_fine_1_no_bal_no_comp, cell.info = PBMC_fine_label_1)






## PBMC model evaluation 

generating_PBMC_evaresult_xgb_comp <-function(type, model){
  overall_acc_list = c()
  avg_acc_list = c()
  avg_sensitivity_list= c()
  avg_specificity_list= c()
  avg_precision_list = c()
  avg_F1_list= c()
  avg_auc_list= c()
  avg_pr_list = c()
  
  
  for (i in 2:24) {
    
    if(type == "Coarse"){
      load(file = paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_coarse_",  i, ".RData"))
      load(file=paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_coarse_label_", i, ".RData"))
      data_matrix = get(paste0("PBMC_coarse_", i))
      model = model
      testing_label = get(paste0("PBMC_coarse_label_", i)) 
      
    }else{
      load(file = paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_fine_",  i, ".RData"))
      load(file=paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_fine_label_", i, ".RData"))
      data_matrix = get(paste0("PBMC_fine_", i))
      model =  model
      testing_label = get(paste0("PBMC_fine_label_", i)) 
    }
    
    
    

    if(type == "Coarse"){
        gene = colnames(PBMC_coarse_1_no_bal_no_comp)}else{
        gene = colnames(PBMC_fine_1_no_bal_no_comp)
        }
    testing_data = as.matrix(t(data_matrix[gene,]))

    pred = t(matrix( predict(model, testing_data, type = "class"), ncol =  dim(testing_data)[1]))

    
    if(type == "Coarse"){
      colnames(pred) = c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK", "other", "other T")
    }else{
      colnames(pred) = c("CD4 CTL", "CD4 Naive", "CD4 TCM",   
                                  "CD4 TEM",  "CD8 Naive", "CD8 TCM",
                                  "CD8 TEM",   "dnT", "gdT","MAIT","Treg") 
    }
    
    pred_label = as.factor(as.character(unlist(apply(pred,1, function(x) sample(names(which(x == max(x))),size=1, replace = F)))))
    cm=table(testing_label,  pred_label)
    
    if(dim(cm)[1] != dim(cm)[2]){
      cm_new = cbind(cm,rep(0, dim(cm)[1]))
      colnames(cm_new) = c(colnames(cm),setdiff(rownames(cm), colnames(cm)))
      cm = cm_new
    }
    
    ## Overall acc
    overall_acc = 0
    for (j in 1:dim(cm)[1]) {
      overall_acc = overall_acc +cm[j,j]/sum(cm)
    }
    
    overall_acc_list = append(overall_acc_list, overall_acc)
    
    ## Average acc
    
    avg_acc = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_acc = avg_acc + ((topleft + bottomright)/(topleft + bottomright+topright + bottomleft)) /dim(cm)[1]
    }
    
    avg_acc_list = append(avg_acc_list, avg_acc)
    
    
    ## Average sensitivity
    
    avg_sensitivity  = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_sensitivity = avg_sensitivity  + ((topleft)/(topleft+topright)) /dim(cm)[1]
    }
    
    avg_sensitivity_list = append(avg_sensitivity_list, avg_sensitivity)
    
    
    ## Average specificity
    
    avg_specificity  = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_specificity = avg_specificity  + ((bottomright)/(bottomleft+bottomright)) /dim(cm)[1]
    }
    
    avg_specificity_list = append(avg_specificity_list, avg_specificity)
    
    
    
    ## Average precision
    
    total_precision  = c() 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      if(topleft+bottomleft==0){
        next
      }
      total_precision = append(total_precision, ((topleft)/(topleft+bottomleft))) 
    }
    
    avg_precision_list = append(avg_precision_list, mean(total_precision))
    
    
    
    ## Average F1 score
    
    avg_F1  = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_F1 = avg_F1  + ((2*topleft)/(2*topleft+bottomleft+topright)) /dim(cm)[1]
    }
    
    avg_F1_list = append(avg_F1_list, avg_F1)
    
    
    
    ## Average AUC/PR
    
    library(pROC)
    library(PRROC)
    
    
    total_auc = c()
    total_pr = c()
    
    for (j in 1:dim(cm)[1]) {
      if(sum(pred[,j]) == 0){
        next
      }
      
      type_name <- names(table(testing_label))
      testing_label_01 = rep(0,length(testing_label))
      testing_label_01[testing_label == type_name[j]  ] = 1
      roc<- roc(testing_label_01, pred[,j])
      pr <- pr.curve(scores.class0 =pred[,j], weights.class0 = testing_label_01,curve = T)
      
      total_auc = append(total_auc, auc(roc))
      total_pr =  append(total_pr, pr$auc.integral) 
      
    }
    
    avg_auc_list = append(avg_auc_list, mean(total_auc))
    avg_pr_list = append(avg_pr_list, mean(total_pr))
    
    
    if(type == "Coarse"){
      rm(list =paste0("PBMC_coarse_", i) )
      rm(list =paste0("PBMC_coarse_label_", i) )
    }else{
      rm(list =paste0("PBMC_fine_", i) )
      rm(list =paste0("PBMC_fine_label_", i) )
    }
    
    gc()
    print(i)
  }
  
  
  
  
  return(list("Overall_Acc"=overall_acc_list,
              "Avg_acc" = avg_acc_list,
              "Avg_sensitivity" = avg_sensitivity_list,
              "Avg_specificity"= avg_specificity_list,
              "Avg_precision"= avg_precision_list,
              "Avg_F1"=avg_F1_list,
              "Avg_AUC"= avg_auc_list,
              "Avg_pr" = avg_pr_list))
  
  
}



result_PBMCcoarse_no_bal_no_comp = generating_PBMC_evaresult_xgb_comp(type = "Coarse", model = xgboost_PBMC1coarse_no_bal_no_comp)
result_PBMCfine_no_bal_no_comp = generating_PBMC_evaresult_xgb_comp(type = "fine", model = xgboost_PBMC1fine_no_bal_no_comp)





save(result_PBMCcoarse_no_bal_no_comp, file = "C:/Users/css22/Desktop/Thesis2/results3/result_PBMCcoarse_no_bal_no_comp.RData")
save(result_PBMCfine_no_bal_no_comp, file = "C:/Users/css22/Desktop/Thesis2/results3/result_PBMCfine_no_bal_no_comp.RData")





########## Only CLR transformation

clr_trans_function <-function(data_matrix = data_matrix_1, error = 0.0001 ){
  data_matrix = data_matrix + error
  center_matrix = matrix(rep(apply(data_matrix,2,function(x) prod(x^(1/length(x)))), each = dim(data_matrix)[1] ), ncol = dim(data_matrix)[2])
  logalt = log(data_matrix/center_matrix)
  return(logalt)
}


PBMC_coarse_1_no_bal_yes_comp<- t(clr_trans_function(data_matrix = PBMC_coarse_1[rowSums(PBMC_coarse_1) != 0,]))
PBMC_fine_1_no_bal_yes_comp<-  t(clr_trans_function(data_matrix = PBMC_fine_1[rowSums(PBMC_fine_1) != 0,]))



xgboost_PBMC1coarse_no_bal_yes_comp = xgboost_train_compare(data1=PBMC_coarse_1_no_bal_yes_comp, cell.info = PBMC_coarse_label_1)
xgboost_PBMC1fine_no_bal_yes_comp = xgboost_train_compare(data1=PBMC_fine_1_no_bal_yes_comp, cell.info = PBMC_fine_label_1)



## PBMC model evaluation 

generating_PBMC_evaresult_xgb_comp_only_clr <-function(type, model){
  overall_acc_list = c()
  avg_acc_list = c()
  avg_sensitivity_list= c()
  avg_specificity_list= c()
  avg_precision_list = c()
  avg_F1_list= c()
  avg_auc_list= c()
  avg_pr_list = c()
  
  
  for (i in 2:24) {
    
    if(type == "Coarse"){
      load(file = paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_coarse_",  i, ".RData"))
      load(file=paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_coarse_label_", i, ".RData"))
      data_matrix = get(paste0("PBMC_coarse_", i))
      model = model
      testing_label = get(paste0("PBMC_coarse_label_", i)) 
      
    }else{
      load(file = paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_fine_",  i, ".RData"))
      load(file=paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_fine_label_", i, ".RData"))
      data_matrix = get(paste0("PBMC_fine_", i))
      model =  model
      testing_label = get(paste0("PBMC_fine_label_", i)) 
    }
    
    
    
    
    if(type == "Coarse"){
      gene = colnames(PBMC_coarse_1_no_bal_no_comp)}else{
        gene = colnames(PBMC_fine_1_no_bal_no_comp)
      }
    data_matrix_clean = data_matrix[gene,]
    testing_data=t(clr_trans_function(data_matrix =  data_matrix_clean)) 
    pred = t(matrix( predict(model, testing_data, type = "class"), ncol =  dim(testing_data)[1]))
    
    
    if(type == "Coarse"){
      colnames(pred) = c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK", "other", "other T")
    }else{
      colnames(pred) = c("CD4 CTL", "CD4 Naive", "CD4 TCM",   
                         "CD4 TEM",  "CD8 Naive", "CD8 TCM",
                         "CD8 TEM",   "dnT", "gdT","MAIT","Treg") 
    }
    
    pred_label = as.factor(as.character(unlist(apply(pred,1, function(x) sample(names(which(x == max(x))),size=1, replace = F)))))
    cm=table(testing_label,  pred_label)
    
    if(dim(cm)[1] != dim(cm)[2]){
      cm_new = cbind(cm,rep(0, dim(cm)[1]))
      colnames(cm_new) = c(colnames(cm),setdiff(rownames(cm), colnames(cm)))
      cm = cm_new
    }
    
    ## Overall acc
    overall_acc = 0
    for (j in 1:dim(cm)[1]) {
      overall_acc = overall_acc +cm[j,j]/sum(cm)
    }
    
    overall_acc_list = append(overall_acc_list, overall_acc)
    
    ## Average acc
    
    avg_acc = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_acc = avg_acc + ((topleft + bottomright)/(topleft + bottomright+topright + bottomleft)) /dim(cm)[1]
    }
    
    avg_acc_list = append(avg_acc_list, avg_acc)
    
    
    ## Average sensitivity
    
    avg_sensitivity  = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_sensitivity = avg_sensitivity  + ((topleft)/(topleft+topright)) /dim(cm)[1]
    }
    
    avg_sensitivity_list = append(avg_sensitivity_list, avg_sensitivity)
    
    
    ## Average specificity
    
    avg_specificity  = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_specificity = avg_specificity  + ((bottomright)/(bottomleft+bottomright)) /dim(cm)[1]
    }
    
    avg_specificity_list = append(avg_specificity_list, avg_specificity)
    
    
    
    ## Average precision
    
    total_precision  = c() 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      if(topleft+bottomleft==0){
        next
      }
      total_precision = append(total_precision, ((topleft)/(topleft+bottomleft))) 
    }
    
    avg_precision_list = append(avg_precision_list, mean(total_precision))
    
    
    
    ## Average F1 score
    
    avg_F1  = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_F1 = avg_F1  + ((2*topleft)/(2*topleft+bottomleft+topright)) /dim(cm)[1]
    }
    
    avg_F1_list = append(avg_F1_list, avg_F1)
    
    
    
    ## Average AUC/PR
    
    library(pROC)
    library(PRROC)
    
    
    total_auc = c()
    total_pr = c()
    
    for (j in 1:dim(cm)[1]) {
      if(sum(pred[,j]) == 0){
        next
      }
      
      type_name <- names(table(testing_label))
      testing_label_01 = rep(0,length(testing_label))
      testing_label_01[testing_label == type_name[j]  ] = 1
      roc<- roc(testing_label_01, pred[,j])
      pr <- pr.curve(scores.class0 =pred[,j], weights.class0 = testing_label_01,curve = T)
      
      total_auc = append(total_auc, auc(roc))
      total_pr =  append(total_pr, pr$auc.integral) 
      
    }
    
    avg_auc_list = append(avg_auc_list, mean(total_auc))
    avg_pr_list = append(avg_pr_list, mean(total_pr))
    
    
    if(type == "Coarse"){
      rm(list =paste0("PBMC_coarse_", i) )
      rm(list =paste0("PBMC_coarse_label_", i) )
    }else{
      rm(list =paste0("PBMC_fine_", i) )
      rm(list =paste0("PBMC_fine_label_", i) )
    }
    
    gc()
    print(i)
  }
  
  
  
  
  return(list("Overall_Acc"=overall_acc_list,
              "Avg_acc" = avg_acc_list,
              "Avg_sensitivity" = avg_sensitivity_list,
              "Avg_specificity"= avg_specificity_list,
              "Avg_precision"= avg_precision_list,
              "Avg_F1"=avg_F1_list,
              "Avg_AUC"= avg_auc_list,
              "Avg_pr" = avg_pr_list))
  
  
}



result_PBMCcoarse_no_bal_yes_comp = generating_PBMC_evaresult_xgb_comp_only_clr(type = "Coarse", model = xgboost_PBMC1coarse_no_bal_yes_comp)
result_PBMCfine_no_bal_yes_comp = generating_PBMC_evaresult_xgb_comp_only_clr(type = "fine", model = xgboost_PBMC1fine_no_bal_yes_comp)





save(result_PBMCcoarse_no_bal_yes_comp, file = "C:/Users/css22/Desktop/Thesis2/results3/result_PBMCcoarse_no_bal_yes_comp.RData")
save(result_PBMCfine_no_bal_yes_comp, file = "C:/Users/css22/Desktop/Thesis2/results3/result_PBMCfine_no_bal_yes_comp.RData")




########## Only Balanced




PBMC_coarse_1_yes_bal_no_comp<- t(PBMC1coarse_trainingdata$data_list[[1]])[,colSums(t(PBMC1coarse_trainingdata$data_list[[1]])) != 0]
PBMC_fine_1_yes_bal_no_comp<- t(PBMC1fine_trainingdata$data_list[[1]])[,colSums(t(PBMC1fine_trainingdata$data_list[[1]])) != 0]




xgboost_PBMC1coarse_yes_bal_no_comp = xgboost_train_compare(data1=PBMC_coarse_1_yes_bal_no_comp, 
                                                            cell.info = PBMC1coarse_trainingdata$Cell_info)

xgboost_PBMC1fine_yes_bal_no_comp = xgboost_train_compare(data1=PBMC_fine_1_yes_bal_no_comp, 
                                                          cell.info = PBMC1fine_trainingdata$Cell_info)


save(xgboost_PBMC1coarse_yes_bal_no_comp, file = "C:/Users/css22/Desktop/Thesis2/Models/xgboost_PBMC1coarse_yes_bal_no_comp.RData")
save(xgboost_PBMC1fine_yes_bal_no_comp, file = "C:/Users/css22/Desktop/Thesis2/Models/xgboost_PBMC1fine_yes_bal_no_comp.RData")


generating_PBMC_evaresult_xgb_only_bal <-function(type, model){
  overall_acc_list = c()
  avg_acc_list = c()
  avg_sensitivity_list= c()
  avg_specificity_list= c()
  avg_precision_list = c()
  avg_F1_list= c()
  avg_auc_list= c()
  avg_pr_list = c()
  
  
  for (i in 2:24) {
    
    if(type == "Coarse"){
      load(file = paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_coarse_",  i, ".RData"))
      load(file=paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_coarse_label_", i, ".RData"))
      data_matrix = get(paste0("PBMC_coarse_", i))
      model = model
      testing_label = get(paste0("PBMC_coarse_label_", i)) 
      
    }else{
      load(file = paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_fine_",  i, ".RData"))
      load(file=paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_fine_label_", i, ".RData"))
      data_matrix = get(paste0("PBMC_fine_", i))
      model =  model
      testing_label = get(paste0("PBMC_fine_label_", i)) 
    }
    
    
    
    
    if(type == "Coarse"){
      gene = colnames(PBMC_coarse_1_yes_bal_no_comp)}else{
        gene = colnames(PBMC_fine_1_yes_bal_no_comp)
      }
    testing_data = as.matrix(t(data_matrix[gene,]))
    
    pred = t(matrix( predict(model, testing_data, type = "class"), ncol =  dim(testing_data)[1]))
    
    
    if(type == "Coarse"){
      colnames(pred) = c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK", "other", "other T")
    }else{
      colnames(pred) = c("CD4 CTL", "CD4 Naive", "CD4 TCM",   
                         "CD4 TEM",  "CD8 Naive", "CD8 TCM",
                         "CD8 TEM",   "dnT", "gdT","MAIT","Treg") 
    }
    
    pred_label = as.factor(as.character(unlist(apply(pred,1, function(x) sample(names(which(x == max(x))),size=1, replace = F)))))
    cm=table(testing_label,  pred_label)
    
    if(dim(cm)[1] != dim(cm)[2]){
      cm_new = cbind(cm,rep(0, dim(cm)[1]))
      colnames(cm_new) = c(colnames(cm),setdiff(rownames(cm), colnames(cm)))
      cm = cm_new
    }
    
    ## Overall acc
    overall_acc = 0
    for (j in 1:dim(cm)[1]) {
      overall_acc = overall_acc +cm[j,j]/sum(cm)
    }
    
    overall_acc_list = append(overall_acc_list, overall_acc)
    
    ## Average acc
    
    avg_acc = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_acc = avg_acc + ((topleft + bottomright)/(topleft + bottomright+topright + bottomleft)) /dim(cm)[1]
    }
    
    avg_acc_list = append(avg_acc_list, avg_acc)
    
    
    ## Average sensitivity
    
    avg_sensitivity  = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_sensitivity = avg_sensitivity  + ((topleft)/(topleft+topright)) /dim(cm)[1]
    }
    
    avg_sensitivity_list = append(avg_sensitivity_list, avg_sensitivity)
    
    
    ## Average specificity
    
    avg_specificity  = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_specificity = avg_specificity  + ((bottomright)/(bottomleft+bottomright)) /dim(cm)[1]
    }
    
    avg_specificity_list = append(avg_specificity_list, avg_specificity)
    
    
    
    ## Average precision
    
    total_precision  = c() 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      if(topleft+bottomleft==0){
        next
      }
      total_precision = append(total_precision, ((topleft)/(topleft+bottomleft))) 
    }
    
    avg_precision_list = append(avg_precision_list, mean(total_precision))
    
    
    
    ## Average F1 score
    
    avg_F1  = 0 
    for (j in 1:dim(cm)[1]) {
      topleft = cm[j,j]
      topright = sum(cm[j,-j])
      bottomleft = sum(cm[-j,j])
      bottomright = sum(cm) - topleft - topright - bottomleft
      avg_F1 = avg_F1  + ((2*topleft)/(2*topleft+bottomleft+topright)) /dim(cm)[1]
    }
    
    avg_F1_list = append(avg_F1_list, avg_F1)
    
    
    
    ## Average AUC/PR
    
    library(pROC)
    library(PRROC)
    
    
    total_auc = c()
    total_pr = c()
    
    for (j in 1:dim(cm)[1]) {
      if(sum(pred[,j]) == 0){
        next
      }
      
      type_name <- names(table(testing_label))
      testing_label_01 = rep(0,length(testing_label))
      testing_label_01[testing_label == type_name[j]  ] = 1
      roc<- roc(testing_label_01, pred[,j])
      pr <- pr.curve(scores.class0 =pred[,j], weights.class0 = testing_label_01,curve = T)
      
      total_auc = append(total_auc, auc(roc))
      total_pr =  append(total_pr, pr$auc.integral) 
      
    }
    
    avg_auc_list = append(avg_auc_list, mean(total_auc))
    avg_pr_list = append(avg_pr_list, mean(total_pr))
    
    
    if(type == "Coarse"){
      rm(list =paste0("PBMC_coarse_", i) )
      rm(list =paste0("PBMC_coarse_label_", i) )
    }else{
      rm(list =paste0("PBMC_fine_", i) )
      rm(list =paste0("PBMC_fine_label_", i) )
    }
    
    gc()
    print(i)
  }
  
  
  
  
  return(list("Overall_Acc"=overall_acc_list,
              "Avg_acc" = avg_acc_list,
              "Avg_sensitivity" = avg_sensitivity_list,
              "Avg_specificity"= avg_specificity_list,
              "Avg_precision"= avg_precision_list,
              "Avg_F1"=avg_F1_list,
              "Avg_AUC"= avg_auc_list,
              "Avg_pr" = avg_pr_list))
  
  
}


result_PBMCcoarse_yes_bal_no_comp = generating_PBMC_evaresult_xgb_only_bal(type = "Coarse", model = xgboost_PBMC1coarse_yes_bal_no_comp)
result_PBMCfine_yes_bal_no_comp = generating_PBMC_evaresult_xgb_only_bal(type = "fine", model = xgboost_PBMC1fine_yes_bal_no_comp)


round(unlist(lapply(result_PBMCfine_yes_bal_no_comp, function(x){mean(x, na.rm =T)} ) )*100,1)
round(unlist(lapply(result_PBMCfine_yes_bal_no_comp, function(x){sd(x, na.rm =T)} ) )*100,1)


