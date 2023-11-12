




### RF models training





rf_train<-function(data, N_layer = 5, nn_pca_number = 100){
  
  
  library(randomForest)
  library(caret)
  library(irlba)
  
  
  pca_model =  lapply(data$data_list, function(x) prcomp_irlba(t(x), n = nn_pca_number, retx = TRUE, center = TRUE, scale. = TRUE)) 
  
  
  rf_model <- list()
  for(n in 1:N_layer){
    
    train_pca <- predict(pca_model[[n]], newdata =t(data1[[n]]))
    embedded_data <- as.data.frame(cbind(data$Cell_info,  as.data.frame(train_pca)))
    colnames(embedded_data)[1] = "CT"
    tuneGrid <- expand.grid(.mtry = seq(10,30,1) )
    rf_cv <- train(CT~.,
                   data = embedded_data,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   importance = F,
                   ntree = 500,
                   trControl =trainControl(method = "cv", number = 10, search ="grid"))
    rf_model[[n]] <- rf_cv$finalModel
  }
  
  out_list = list()
  out_list$pca_model <- pca_model
  out_list$rf_model <- rf_model
  
  return(out_list)
}
rf_train_binary<-function(data, N_layer = 5, nn_pca_number = 100, celltype){
  
  
  library(randomForest)
  library(caret)
  library(irlba)
  
  cell.info <- data$Cell_info
  cell.info[ cell.info != celltype] = "Others"

  pca_model =  lapply(data$data_list, function(x) prcomp_irlba(t(x), n = nn_pca_number, retx = TRUE, center = TRUE, scale. = TRUE)) 
  
  
  rf_model <- list()
  for(n in 1:N_layer){
    
    train_pca <- predict(pca_model[[n]], newdata =t(data1[[n]]))
    embedded_data <- as.data.frame(cbind(cell.info,  as.data.frame(train_pca)))
    colnames(embedded_data)[1] = "CT"
    tuneGrid <- expand.grid(.mtry = seq(10,30,1) )
    rf_cv <- train(CT~.,
                   data = embedded_data,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   importance = F,
                   ntree = 500,
                   trControl =trainControl(method = "cv", number = 10, search ="grid"))
    rf_model[[n]] <- rf_cv$finalModel
  }
  
  out_list = list()
  out_list$pca_model <- pca_model
  out_list$rf_model <- rf_model
  
  return(out_list)
}





### Shrinkage models training

library(glmnet)
glmnet_train <- function(data, alpha, celltype){
  centered_data = lapply(data$data_list, function(x) scale(t(x)))
  cell.info = data$Cell_info
  model = lapply(centered_data, function(x) cv.glmnet(x,cell.info, family = "multinomial", alpha =alpha,  type.measure = "class", nfolds = 10))
  return(model)
}
glmnet_train_binary <- function(data, alpha){
  centered_data = lapply(data$data_list, function(x) scale(t(x)))
  cell.info <- data$Cell_info
  cell.info[ cell.info != celltype] = "Others"
  model = lapply(centered_data, function(x) cv.glmnet(x,cell.info, family = "binomial", alpha =alpha,  type.measure = "auc", nfolds = 10))
  return(model)
}
  



## Kernel KNN for multi-class classification model training

kknn_get_best_paras <- function(data, N_layer = 5, nn_pca_number = 30){
  
  library(MASS)
  library(KernelKnn)
  library(irlba)
  
  
  acc = function (y_true, preds) {
    
    out = table(y_true, max.col(preds, ties.method = "random"))
    
    acc = sum(diag(out))/sum(out)
    
    acc
  }
  
  
  
  pca_model =  lapply(data$data_list, function(x) prcomp_irlba(t(x), n = nn_pca_number, retx = TRUE, center = TRUE, scale. = TRUE)) 
  
  
  

  label = as.integer(factor(data$Cell_info))

  
  best_para_list = list()
  for (n in 1:N_layer) {
    
    acc_max = 0
    k_list = c(3:20)
    method_list = c("manhattan", "euclidean", "minkowski", "chebyshev")
    kernel_list = c('uniform', 'triangular', 'epanechnikov', 'biweight', 'triweight', 
                    "tricube", 'gaussian', 'cosine', 'logistic', 'gaussianSimple', 
                    'silverman', 'inverse', 'exponential' ) 
    train_pca <- predict(pca_model[[n]], newdata =t(data1[[n]]))
    for(k in k_list){
      for( method in method_list){
        for(kernel in kernel_list){
          fit_cv = KernelKnnCV(train_pca, label, k = k , folds = 10, method = method, 
                               
                               weights_function = kernel , regression = F, Levels = unique(label), threads = 5)
          
          mean_acc = mean(unlist(lapply(1:length(fit_cv$preds), function(x) acc(label[fit_cv$folds[[x]]], fit_cv$preds[[x]]))))
          
          
          if(mean_acc>acc_max){
            acc_max = mean_acc
            best_para = c(k, method, kernel)
          }          
          
        }}}
    best_para_list[[n]] = best_para
    
  }
  
  out_list = list()
  out_list$pca_model <- pca_model
  out_list$best_para <- best_para_list
  
  return(out_list)
  
}
  
## Kernel SVM for binary classification model training



svm_train_binary <- function(data, N_layer = 5, nn_pca_number = 100, celltype){
  
  
  
  library(e1071)
  library(caret)
  
  pca_model =  lapply(data$data_list, function(x) prcomp_irlba(t(x), n = nn_pca_number, retx = TRUE, center = TRUE, scale. = TRUE)) 
  
  
  cell.info <- data$Cell_info
  cell.info[ cell.info != celltype] = "Others"
  cell.info = as.integer(as.factor(cell.info))
  svm_model = list()
  
  for(n in 1:N_layer){
    train_pca <- predict(pca_model[[n]], newdata =t(data$data_list[[n]]))
    train_data <- as.data.frame(cbind(cell.info,  as.data.frame(train_pca)))
    folds = createFolds(train_data$cell.info, k = 10)
    kernel = c("radial", "polynomial", "sigmoid")
    gamma  = c(0.001, 0.01, 0.1, 1, 10)
    C      = c(0.01,0.1,1,10,100)
    acc_max = 0
    for (k in 1:3) {
      for (g  in 1:5) {
        for(c in 1:5)  {
          
          cv = lapply(folds, function(x) { 
            training_fold = train_data[-x, ] 
            test_fold = train_data[x, ] 
            classifier = svm(formula = cell.info ~ .,
                             data = training_fold,
                             type = 'C-classification',
                             kernel = kernel[k],
                             gamma = gamma[g],
                             cost = C[c],
                             probability = TRUE)
            y_pred = predict(classifier, newdata = test_fold[,-1])
            cm = table(test_fold[, 1], y_pred)
            accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
            return(accuracy)
          })
          
          
          mean_acc = mean(unlist(cv))
          
          if(mean_acc>acc_max){
            acc_max = mean_acc
            best_para = c(kernel[k], gamma[g], C[c])
          }  
        }
      }
    }
    
    svm =   svm(formula = cell.info ~ .,
                data   = train_data,
                type   = 'C-classification',
                kernel = best_para[1],
                gamma  = best_para[2],
                cost   = best_para[3],
                probability = TRUE)
    
    svm_model[[n]] <- svm
  }
  
  
  
  out_list = list()
  out_list$pca_model <- pca_model
  out_list$svm_model <- svm_model
  
  return(out_list)
  
}




## Neural network (MLP) model training



MLP_train <- function( data, N_layer =5,  nn_pca_number = 100, hidden, rep = 3){
  
  
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }
  
  
  pca_model =  lapply(data$data_list, function(x) prcomp_irlba(t(x), n = nn_pca_number, retx = TRUE, center = TRUE, scale. = TRUE)) 
  
  
  library(neuralnet)
  
  nn_model <- list()
  for(n in 1:N_layer){
    
    train_pca <- predict(pca_model[[n]], newdata =t(data1[[n]]))
    train_data <- apply(train_pca,2,normalize)
    
    embedded_data <- as.data.frame(cbind(data$Cell_info,  as.data.frame(train_data)))
    colnames(embedded_data)[1] = "CT"
    nn_model[[n]] <- neuralnet(CT~., 
                               data =embedded_data, lifesign = "minimal",
                               stepmax = 3e+05, rep= rep, 
                               hidden= hidden, linear.output=F, threshold=0.01, 
                               act.fct = "logistic", err.fct = "ce")
  }
  
  
  
  out_list = list()
  out_list$pca_model = pca_model
  out_list$nn_model = nn_model
  
  return(out_list)
  
  
}
MLP_train_binary <- function(data, N_layer =5,  nn_pca_number = 100, hidden, rep = 3, celltype){
  
  cell.info <- data$Cell_info
  cell.info[cell.info != celltype] = "Others"
  
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }
  
  
  pca_model =  lapply(data$data_list, function(x) prcomp_irlba(t(x), n = nn_pca_number, retx = TRUE, center = TRUE, scale. = TRUE)) 
  
  
  library(neuralnet)
  
  nn_model <- list()
  for(n in 1:N_layer){
    
    train_pca <- predict(pca_model[[n]], newdata =t(data1[[n]]))
    train_data <- apply(train_pca,2,normalize)
    
    embedded_data <- as.data.frame(cbind(cell.info,  as.data.frame(train_data)))
    colnames(embedded_data)[1] = "CT"
    nn_model[[n]] <- neuralnet(CT~., 
                               data =embedded_data, lifesign = "minimal",
                               stepmax = 3e+05, rep= rep, 
                               hidden= hidden, linear.output=F, threshold=0.01, 
                               act.fct = "logistic", err.fct = "ce")
  }
  
  
  
  out_list = list()
  out_list$pca_model = pca_model
  out_list$nn_model = nn_model
  
  return(out_list)
  
  
}
