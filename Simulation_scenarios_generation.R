rm(list=ls())
gc()


library(devtools)
install_gitlab("sysbiobig/sparsim")


library(SPARSim)

data(Zheng_param_preset)

sparsim_gene_name<-NULL
for (i in 1:15000) {
  sparsim_gene_name<-append(sparsim_gene_name, paste0("Gene", i ))
}

sparsim_cell_name<-NULL
for (i in 1:6000) {
  sparsim_cell_name<-append(sparsim_cell_name, paste0("Cell", i ))
}



### new coarse
set.seed(1234)
downreg.min.set = rnorm(24, mean = 0.3, sd = 0.05)
downreg.max.set = rnorm(24, mean = 0.5, sd = 0.05)
upreg.min.set   = rnorm(24, mean = 1.8, sd = 0.05)
upreg.max.set   = rnorm(24, mean = 2.0, sd = 0.05)


### new fine
set.seed(1234)
downreg.min.set = rnorm(24, mean = 0.6, sd = 0.05)
downreg.max.set = rnorm(24, mean = 0.8, sd = 0.05)
upreg.min.set   = rnorm(24, mean = 1.4, sd = 0.05)
upreg.max.set   = rnorm(24, mean = 1.6, sd = 0.05)




## Generating results simulation datasets

 
 for (i in 1:20) {
   
 
  set.seed( i*10000 + 1 )
  DE_multiplier_1.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2)
  set.seed( i*10000 + 2 )
  DE_multiplier_1.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 3 )
  DE_multiplier_1.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 4 )
  DE_multiplier_1.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  fold_change_multiplier_1 <- c(DE_multiplier_1.1, rep(1, 180), DE_multiplier_1.2, rep(1,400),
                                DE_multiplier_1.3, rep(1,400), DE_multiplier_1.4, rep(1,17936) )
  
  
  set.seed( i*10000 + 5 )
  DE_multiplier_2.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
  set.seed( i*10000 + 6 )
  DE_multiplier_2.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 7 )
  DE_multiplier_2.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 8 )
  DE_multiplier_2.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 9 )
  DE_multiplier_2.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  fold_change_multiplier_2 <- c(rep(1,20), DE_multiplier_2.1, rep(1,360), DE_multiplier_2.2,
                                rep(1,200), DE_multiplier_2.3, rep(1,200), DE_multiplier_2.4, 
                                rep(1,200), DE_multiplier_2.5, rep(1,17736) )
  
  
  set.seed( i*10000 + 10 )
  DE_multiplier_3.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
  set.seed(  i*10000 + 11 )
  DE_multiplier_3.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 12 )
  DE_multiplier_3.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 13 )
  DE_multiplier_3.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 14 )
  DE_multiplier_3.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  fold_change_multiplier_3 <- c(rep(1,40), DE_multiplier_3.1, rep(1,140), DE_multiplier_3.2,
                                DE_multiplier_3.3, rep(1,800), DE_multiplier_3.4, 
                                rep(1,200), DE_multiplier_3.5, rep(1,17536) )
  
  
  
  set.seed( i*10000 + 15 )
  DE_multiplier_4.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
  set.seed( i*10000 + 16 )
  DE_multiplier_4.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 17 )
  DE_multiplier_4.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 18 )
  DE_multiplier_4.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 19 )
  DE_multiplier_4.5 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 20 )
  DE_multiplier_4.6 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  fold_change_multiplier_4 <- c(rep(1,60), DE_multiplier_4.1, rep(1,320), DE_multiplier_4.2, rep(1,400),
                                DE_multiplier_4.3, DE_multiplier_4.4, DE_multiplier_4.5, 
                                rep(1,200), DE_multiplier_4.6, rep(1,17536) )
  
  
  set.seed( i*10000 + 21 )
  DE_multiplier_5.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
  set.seed( i*10000 + 22 )
  DE_multiplier_5.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 23 )
  DE_multiplier_5.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 24 )
  DE_multiplier_5.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 25 )
  DE_multiplier_5.5 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  fold_change_multiplier_5 <- c(rep(1,80), DE_multiplier_5.1, rep(1,100), DE_multiplier_5.2,
                                rep(1,200), DE_multiplier_5.3, rep(1,800), DE_multiplier_5.4, 
                                DE_multiplier_5.5, rep(1,17536) )
  
  
  
  
  set.seed( i*10000 + 26 )
  DE_multiplier_6.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
  set.seed( i*10000 + 27 )
  DE_multiplier_6.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 28 )
  DE_multiplier_6.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 29 )
  DE_multiplier_6.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 30 )
  DE_multiplier_6.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  fold_change_multiplier_6 <- c(rep(1,100), DE_multiplier_6.1, rep(1,280), DE_multiplier_6.2,
                                rep(1,400), DE_multiplier_6.3, rep(1,200), DE_multiplier_5.4, 
                                rep(1,200), DE_multiplier_6.5, rep(1,17536) )
  
  
  
  set.seed( i*10000 + 31 )
  DE_multiplier_7.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
  set.seed( i*10000 + 32 )
  DE_multiplier_7.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 33 )
  DE_multiplier_7.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 34 )
  DE_multiplier_7.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  fold_change_multiplier_7 <- c(rep(1,120), DE_multiplier_7.1, rep(1,60), DE_multiplier_7.2,
                                rep(1,200), DE_multiplier_7.3, rep(1,400), DE_multiplier_7.4, rep(1,18136) )
  
  
  
  set.seed( i*10000 + 35 )
  DE_multiplier_8.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
  set.seed( i*10000 + 36 )
  DE_multiplier_8.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 37 )
  DE_multiplier_8.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 38 )
  DE_multiplier_8.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  fold_change_multiplier_8 <- c(rep(1,140), DE_multiplier_8.1, rep(1,440), DE_multiplier_8.2,
                                rep(1,200), DE_multiplier_8.3, rep(1,400), DE_multiplier_8.4, rep(1,17736) )
  
  
  
  
  
  set.seed( i*10000 + 39 )
  DE_multiplier_9.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
  set.seed( i*10000 + 40 )
  DE_multiplier_9.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 41 )
  DE_multiplier_9.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 42 )
  DE_multiplier_9.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 43 )
  DE_multiplier_9.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  fold_change_multiplier_9 <- c(rep(1,160), DE_multiplier_9.1, rep(1,20), DE_multiplier_9.2,
                                rep(1,800), DE_multiplier_9.3, rep(1,200), DE_multiplier_9.4, 
                                DE_multiplier_9.5, rep(1,17536) )
  
  
  
  
  set.seed( i*10000 + 44 )
  DE_multiplier_10.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
  set.seed( i*10000 + 45 )
  DE_multiplier_10.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 46 )
  DE_multiplier_10.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
  set.seed( i*10000 + 47 )
  DE_multiplier_10.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  set.seed( i*10000 + 48 )
  DE_multiplier_10.5 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
  fold_change_multiplier_10 <- c(rep(1,180), DE_multiplier_9.1, rep(1,600), DE_multiplier_10.2,
                                 DE_multiplier_10.3, rep(1,200), DE_multiplier_10.4, 
                                 DE_multiplier_10.5, rep(1,17736) )
  
  
  
  p1 = 60
  cell_1_list = c(p1, 2*p1, 4*p1, 6*p1, 10*p1)
  cell_2_list = c(3*p1, 4*p1, 6*p1, 7*p1, 10*p1)
  cell_3_list = c(4*p1, 6*p1, 7*p1, 8*p1, 10*p1)
  cell_4_list = c(6*p1, 8*p1, 8*p1, 9*p1, 10*p1)
  cell_5_list = c(7*p1, 8*p1, 8*p1, 10*p1, 10*p1)
  cell_6_list = c(8*p1, 10*p1, 10*p1, 10*p1, 10*p1)
  cell_7_list = c(10*p1, 10*p1, 10*p1, 11*p1, 10*p1)
  cell_8_list = c(16*p1, 14*p1, 12*p1, 12*p1, 10*p1)
  cell_9_list = c(20*p1, 16*p1, 15*p1, 13*p1, 10*p1)
  cell_10_list = c(25*p1, 22*p1, 20*p1, 14*p1, 10*p1)
  
  
  j = 5
  
  cell_type_1 <- SPARSim_create_DE_genes_parameter(
    sim_param =Zheng_param_preset$Zheng_C2, 
    fc_multiplier = fold_change_multiplier_1, 
    N_cells = cell_1_list[j],
    condition_name = "cell_1")
  
  
  cell_type_2 <- SPARSim_create_DE_genes_parameter(
    sim_param =Zheng_param_preset$Zheng_C2, 
    fc_multiplier = fold_change_multiplier_2, 
    N_cells = cell_2_list[j],
    condition_name = "cell_2")
  
  
  
  cell_type_3 <- SPARSim_create_DE_genes_parameter(
    sim_param =Zheng_param_preset$Zheng_C2 , 
    fc_multiplier = fold_change_multiplier_3, 
    N_cells = cell_3_list[j],
    condition_name = "cell_3")
  
  
  cell_type_4 <- SPARSim_create_DE_genes_parameter(
    sim_param =Zheng_param_preset$Zheng_C2 , 
    fc_multiplier = fold_change_multiplier_4, 
    N_cells = cell_4_list[j],
    condition_name = "cell_4")
  
  
  cell_type_5 <- SPARSim_create_DE_genes_parameter(
    sim_param =Zheng_param_preset$Zheng_C2 , 
    fc_multiplier = fold_change_multiplier_5, 
    N_cells = cell_5_list[j],
    condition_name = "cell_5")
  
  
  cell_type_6 <- SPARSim_create_DE_genes_parameter(
    sim_param =Zheng_param_preset$Zheng_C2 , 
    fc_multiplier = fold_change_multiplier_6, 
    N_cells = cell_6_list[j],
    condition_name = "cell_6")
  
  
  cell_type_7 <- SPARSim_create_DE_genes_parameter(
    sim_param =Zheng_param_preset$Zheng_C2 , 
    fc_multiplier = fold_change_multiplier_7, 
    N_cells = cell_7_list[j],
    condition_name = "cell_7")
  
  
  cell_type_8 <- SPARSim_create_DE_genes_parameter(
    sim_param =Zheng_param_preset$Zheng_C2 , 
    fc_multiplier = fold_change_multiplier_8, 
    N_cells = cell_8_list[j],
    condition_name = "cell_8")
  
  
  
  cell_type_9 <- SPARSim_create_DE_genes_parameter(
    sim_param =Zheng_param_preset$Zheng_C2 , 
    fc_multiplier = fold_change_multiplier_9, 
    N_cells = cell_9_list[j],
    condition_name = "cell_9")
  
  
  
  
  cell_type_10 <- SPARSim_create_DE_genes_parameter(
    sim_param =Zheng_param_preset$Zheng_C2 , 
    fc_multiplier = fold_change_multiplier_10, 
    N_cells = cell_10_list[j],
    condition_name = "cell_10")
  
  
  SPARSim_param_with_DE <- list(cell_type_1=cell_type_1,
                                cell_type_2=cell_type_2,
                                cell_type_3=cell_type_3,   
                                cell_type_4=cell_type_4,
                                cell_type_5=cell_type_5,
                                cell_type_6=cell_type_6,
                                cell_type_7=cell_type_7,   
                                cell_type_8=cell_type_8,
                                cell_type_9=cell_type_9,
                                cell_type_10=cell_type_10)
  
  ## Can not be seeded (SPARSIM)
  
  SPARSim_result <- SPARSim_simulation(SPARSim_param_with_DE )
  
  set.seed(1234*i)
  ind<-sample(c(2001:19536), size = 4536 ,replace = F)
  
  sparsim.data=SPARSim_result$count_matrix[-ind,]
  
  rownames(sparsim.data)=sparsim_gene_name
  colnames(sparsim.data)=sparsim_cell_name
  
  assign(paste0("Sim_test_", i), sparsim.data)
  
  save(list= paste0("Sim_test_",i), 
       file=paste0("C:/Users/css22/Desktop/Thesis2/SPARSim/Sim_test_", i,  ".RData" ))
  
  rm(list = paste0("Sim_test_", i))
  rm(sparsim.data)
  gc()
  
  print(i)
  
}
  

load(file = "C:/Users/css22/Desktop/Thesis2/SPARSim/Sim_train_1.RData")  
load(file = "C:/Users/css22/Desktop/Thesis2/SPARSim/Sim_train_2.RData")  
load(file = "C:/Users/css22/Desktop/Thesis2/SPARSim/Sim_train_3.RData")  
load(file = "C:/Users/css22/Desktop/Thesis2/SPARSim/Sim_train_4.RData")  
load(file = "C:/Users/css22/Desktop/Thesis2/SPARSim/Sim_train_5.RData")  



dim(Sim_train_1)
dim(Sim_train_2)
dim(Sim_train_3)
dim(Sim_train_4)
dim(Sim_train_5)

  
  













set.seed(1234)
seurat.obj<-CreateAssayObject(counts = SparSim_fine_aim2_1[, ] )
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)


m = as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)



umap<-uwot::umap(pc_df , n_components=2)
umap.df<-as.data.frame(umap)

cell.info_sim<-c(rep("Type 1",  60),
                 rep("Type 2",  180),
                 rep("Type 3",  240),
                 rep("Type 4",  360),
                 rep("Type 5",  420),
                 rep("Type 6",  480),
                 rep("Type 7",  600),
                 rep("Type 8",  960),
                 rep("Type 9",  1200),
                 rep("Type 10", 1500))

umap.df<-cbind(umap.df, cell.info_sim)
colnames(umap.df)[3]<-"group"


b=ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.9, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2", title = "Simulated Data (Coarse, Imbalanced)")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "right",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))






library(patchwork)
a+b








