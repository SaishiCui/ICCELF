
library(compiler)


## data: scRNA-seq data matrix (or data frame), columns are cells, rows are gene, row names should be the gene name.
## simulated_N: the number for Monte-Carlo simulation, default is 5,000.
## balanced_N: the sample size of each balanced cell type.
## N_layer: the number of data layers to be generated.
## error: the error to be added in order to do centered log ratio.
## cell_label: cell labels corresponding to columns of the input dataset.


iccelf <- function(data, simulated_N = 5000, balanced_N, N_layer = 5, error = 0.0001, cell_label){

  
## Oversampling function
oversampling<-function(data, simulated_N, balanced_N){
  n = dim(data)[2]
  aug_data = matrix(NA, nrow = dim(data)[1], ncol = balanced_N - n  )
  for (gene in 1:dim(data)[1]) {
    
    if(gene %% 1000 == 0){
      print( paste0( round((gene/1000)/ (dim(data)[1]/1000) * 100,2), "%"))
    }
    
    
    lambda = mean(as.numeric(data[gene,]))
    if(lambda == 0){
      aug_data[gene,]= rep(0, balanced_N - n)
      next
    }
    
    sample_variance = sd(as.numeric(data[gene,]))^2
    
    sv_list = c()
    for(i in 1:balanced_N){
      simulated_sample = rpois(n,lambda)
      simulated_variance = sd(simulated_sample)^2
      sv_list = c(sv_list, simulated_variance )
    }
    p_value <- sum(sv_list > sample_variance)/balanced_N
    if(p_value >= 0.05){
      aug_data[gene,]= rpois(balanced_N - n , lambda )
    }else{
      
      p = lambda/sample_variance
      r= lambda^2/(sample_variance - lambda)
      aug_data[gene,]=rnbinom(balanced_N - n, r, p)
    }
    
  }
  output_data = cbind(data, aug_data)
  return(output_data)
}


## generate N_data layers' datasets 

generate_stacking_model_data<- function(input_data, simulated_N, balanced_N = balanced_N, 
                                        input_gene, cell.info, N_layer = N_layer){
  
  celltype = unique(cell.info)
  N_celltype = length(unique(cell.info))
  new_data_list = list()
  new_cell.info = c()
  for(j in 1:N_layer){
    new_data = matrix(NA, nrow = length(input_gene), ncol = balanced_N*N_celltype)
    for (i in 1:N_celltype) {
      if(sum(cell.info == celltype[i]) < balanced_N){
        data = input_data[input_gene, cell.info == celltype[i]]
        generate_data =  oversampling(data= data, simulated_N = simulated_N, balanced_N = balanced_N)
        new_data[, (1+balanced_N*(i-1)):(balanced_N*i)] = as.matrix(generate_data)
        if(j == 1){
          new_cell.info = append(new_cell.info, rep(celltype[i], balanced_N))
        }
      }else{
        data = input_data[input_gene, cell.info == celltype[i]]
        downsample_index <- sample(sum(cell.info == celltype[i]), balanced_N, replace = F)
        new_data[, (1+balanced_N*(i-1)):(balanced_N*i)] = as.matrix(data[,downsample_index])
        if(j == 1){
          new_cell.info = append(new_cell.info, rep(celltype[i], balanced_N))
        }
      }
      cat("Cell type generate:", i, "(total:", N_celltype, "types)", "\n")
    }
    
    rownames(new_data) = input_gene 
    colnames(new_data) = paste0("Cell ", 1:(balanced_N*N_celltype) )
    
    
    new_data_list[[j]] = new_data
    cat("Data layer generate:", j, "\n")
  }
  
  output_list = list()
  output_list$"Cell_info" = new_cell.info
  output_list$"data_list" = new_data_list
  return(output_list)
}


## Clr transformation

clr_trans_function <-function(data, error ){
  data = data + error
  center_matrix = matrix(rep(apply(data,2,function(x) prod(x^(1/length(x)))), each = dim(data)[1] ), ncol = dim(data)[2])
  logalt = log(data/center_matrix)
  return(logalt)
}




trainingdata = generate_stacking_model_data(input_data   =  data,
                                            simulated_N  =  simulated_N,
                                            balanced_N   =  balanced_N,                                                  
                                            input_gene   =  rownames(data),
                                            cell.info    =  cell_label,
                                            N_layer      =  N_layer)

## clr transformation for each part (true cells part and synthetic cells part and concatenate them)

clr_trainingdata = list()
for (N in 1:N_layer) {
  check_zero = apply(trainingdata$data_list[[N]],1,sum)
  zero_idx = as.numeric(which(check_zero == 0))
  clr_trainingdata[[N]] = trainingdata$data_list[[N]][-zero_idx,]
  
  true_cells_label = c()
  for (c in 1:length(unique(cell_label))) {
   cell_number <- as.numeric(table(cell_label)[c])
   if(cell_number < balanced_N){
     true_cell_position <- c((1+balanced_N*(c-1)):(balanced_N*(c-1) + cell_number))
     true_cells_label = append(true_cells_label, true_cell_position)
   }else{
     
     true_cell_position <- c((1+balanced_N*(c-1)):(balanced_N*c))
     true_cells_label = append(true_cells_label, true_cell_position)
   }
  }

  clr_trainingdata[[N]][,true_cells_label]   = clr_trans_function(data = clr_trainingdata[[N]][,true_cells_label], error = error)
  clr_trainingdata[[N]][,-true_cells_label]  = clr_trans_function(data = clr_trainingdata[[N]][,-true_cells_label], error = error)
}

out_list = list()
out_list$"Cell_info" = trainingdata$Cell_info
out_list$"data_list" = clr_trainingdata
return(out_list)

}


iccelf <- cmpfun(iccelf)