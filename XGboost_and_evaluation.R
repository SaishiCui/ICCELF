rm(list=ls())
gc()

load(file = "C:/Users/css22/Desktop/Thesis2/Training_data/Sim1coarse_trainingdata.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Training_data/Sim1coarse_clr_trainingdata.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Training_data/Sim1fine_trainingdata.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Training_data/Sim1fine_clr_trainingdata.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Training_data/PBMC1coarse_trainingdata.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Training_data/PBMC1coarse_clr_trainingdata.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Training_data/PBMC1fine_trainingdata.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Training_data/PBMC1fine_clr_trainingdata.RData")




load(file = "C:/Users/css22/Desktop/Thesis2/Models/xgboost_sim1coarse_finalstack.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Models/xgboost_sim1fine_finalstack.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Models/xgboost_PBMC1coarse_finalstack.RData")
load(file = "C:/Users/css22/Desktop/Thesis2/Models/xgboost_PBMC1fine_finalstack.RData")





### XGboost model (Multi-class classification)

library(xgboost)

N_layer = 5
xgboost_train <- function(data){
  
  if(data$Cell_info[1] == "Type 1"){
  label = as.integer(factor(data$Cell_info, levels= paste("Type", c(1:10)) ))-1
  }else{
    label = as.integer(factor(data$Cell_info))-1
  }
  
  dtrain <- lapply(data$data_list, function(x) xgb.DMatrix(data = as.matrix(t(x)), label =  label))    
  params <- list(booster = "gbtree", objective = "multi:softprob",   eval_metric="mlogloss", num_class= length(unique(data$Cell_info)),
                 eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample= 0.5, colsample_bytree= 0.5)
  
  xgbcv <- lapply(dtrain, function(x) xgb.cv( params = params, data = x, nrounds = 200, nfold = 10, showsd = T, 
                                              stratified = T, print_every_n = 100, early_stop_round = 50, maximize = F))
  
  nrounds = lapply(xgbcv, function(x) which(x$evaluation_log[,4] ==  min(x$evaluation_log[,4]) ))
  
  model = list()
  for(n in 1:N_layer){
    model[[n]] <- xgb.train (params = params, data = dtrain[[n]], nrounds = nrounds[[n]], print_every_n = 50, maximize = F)
  }
  
  return(model)
}



xgboost_PBMC1coarse_finalstack = xgboost_train(data = PBMC1coarse_trainingdata)
xgboost_PBMC1fine_finalstack = xgboost_train(data = PBMC1fine_trainingdata)

  ## Simulation (5 scenarios)
xgboost_Sim1_finalstack = xgboost_train(data = Sim_trainingdata1)
xgboost_Sim2_finalstack = xgboost_train(data = Sim_trainingdata2)
xgboost_Sim3_finalstack = xgboost_train(data = Sim_trainingdata3)
xgboost_Sim4_finalstack = xgboost_train(data = Sim_trainingdata4)
xgboost_Sim5_finalstack = xgboost_train(data = Sim_trainingdata5)


### XGboost model (Binary classification)

library(xgboost)

N_layer = 5
xgboost_train_binary <- function(data, cell_type){
  
  cell.info <- data$Cell_info
  cell.info[ cell.info != cell_type] = "Others"
  
  label = as.integer( as.factor(cell.info))-1
  
  
  dtrain <- lapply(data, function(x) xgb.DMatrix(data = as.matrix(t(x)), label =  label))    
  params <- list(booster = "gbtree", objective = "binary:logistic",   eval_metric="logloss",
                 eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample= 0.5, colsample_bytree= 0.5)
  
  xgbcv <- lapply(dtrain, function(x) xgb.cv( params = params, data = x, nrounds = 500, nfold = 10, showsd = T, 
                                              stratified = T, print_every_n = 100, early_stop_round = 20, maximize = F))
  
  nrounds = lapply(xgbcv, function(x) which(x$evaluation_log[,4] ==  min(x$evaluation_log[,4]) ))
  
  model = list()
  for(n in 1:N_layer){
    model[[n]] <- xgb.train (params = params, data = dtrain[[n]], nrounds = nrounds[[n]], print_every_n = 100, maximize = F)
  }
  
  return(model)
  
}

xgboost_PBMC1coarse_finalstack_binary = xgboost_train_binary(data = PBMC1coarse_trainingdata, cell_type = "B")
xgboost_PBMC1fine_finalstack_binary = xgboost_train_binary(data = PBMC1fine_trainingdata, cell_type = "gdT")




## PBMC model evaluation 

generating_PBMC_evaresult_xgb <-function(type, model){
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
    
    
    
    total_pred = matrix(0, nrow = dim(data_matrix)[2], ncol = length(unique(testing_label)))
    for(n in 1:5){
      if(type == "Coarse"){
        gene = rownames(PBMC1coarse_clr_trainingdata[[n]])}else{
          gene = rownames(PBMC1fine_clr_trainingdata[[n]])
        }
      data_matrix_clean = data_matrix[gene,]
      testing_data=t(clr_trans_function(data_matrix =  data_matrix_clean)) 
      pred = t(matrix( predict(model[[n]], testing_data, type = "class"), ncol =  dim(testing_data)[1]))
      total_pred = total_pred +pred
    }
    stacking_pred = total_pred / 5
    
    if(type == "Coarse"){
      colnames(stacking_pred) = c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK", "other", "other T")
    }else{
      colnames(stacking_pred) = c("CD4 CTL", "CD4 Naive", "CD4 TCM",   
                                  "CD4 TEM",  "CD8 Naive", "CD8 TCM",
                                  "CD8 TEM",   "dnT", "gdT","MAIT","Treg") 
    }
    
    pred_label = as.factor(as.character(unlist(apply(stacking_pred,1, function(x) sample(names(which(x == max(x))),size=1, replace = F)))))
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
      if(sum(stacking_pred[,j]) == 0){
        next
      }
      
      type_name <- names(table(testing_label))
      testing_label_01 = rep(0,length(testing_label))
      testing_label_01[testing_label == type_name[j]  ] = 1
      roc<- roc(testing_label_01, stacking_pred[,j])
      pr <- pr.curve(scores.class0 =stacking_pred[,j], weights.class0 = testing_label_01,curve = T)
      
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



result_PBMCcoarse_xgb = generating_PBMC_evaresult_xgb(type = "Coarse", model = xgboost_PBMC1coarse_finalstack)
result_PBMCfine_xgb = generating_PBMC_evaresult_xgb(type = "fine", model = xgboost_PBMC1fine_finalstack)


### For binary classification and simulation models evaluation, the procedure is similar as above. 


