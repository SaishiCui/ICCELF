rm(list=ls())
gc()

load(file = paste0("C:/Users/css22/Desktop/Thesis2/SPARSim/Sim_train_", 1,  ".RData" ))
load(file = paste0("C:/Users/css22/Desktop/Thesis2/SPARSim/Sim_train_", 2,  ".RData" ))
load(file = paste0("C:/Users/css22/Desktop/Thesis2/SPARSim/Sim_train_", 3,  ".RData" ))
load(file = paste0("C:/Users/css22/Desktop/Thesis2/SPARSim/Sim_train_", 4,  ".RData" ))
load(file = paste0("C:/Users/css22/Desktop/Thesis2/SPARSim/Sim_train_", 5,  ".RData" ))


load(file = paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_coarse_",  1, ".RData"))
load(file = paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_coarse_label_",  1, ".RData"))
load(file = paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_fine_",  1, ".RData"))
load(file = paste0("C:/Users/css22/Desktop/Thesis2/PBMC_data/PBMC_fine_label_",  1, ".RData"))




source("C:/Users/css22/Desktop/Thesis2/ICCELF.R")




### Creating training datasets for PBMC (Multi-class classification)



PBMC1coarse_trainingdata = iccelf(data = PBMC_coarse_1,
                                  simulated_N = 5000, 
                                  balanced_N = 400,
                                  N_layer = 5,
                                  error = 0.0001,
                                  cell_label = PBMC_coarse_label_1)




PBMC1fine_trainingdata = iccelf(data = PBMC_fine_1,
                                simulated_N = 5000, 
                                balanced_N = 400,
                                N_layer = 5,
                                error = 0.0001,
                                cell_label = PBMC_fine_label_1)








### Creating training datasets for simulation (5 scenarios)




Sim_trainingdata1 = iccelf(data = Sim_train_1,
                           simulated_N = 5000, 
                           balanced_N = 600,
                           N_layer = 5,
                           error = 0.0001,
                           cell_label = c(rep("Type 1",  60),
                                          rep("Type 2",  180),
                                          rep("Type 3",  240),
                                          rep("Type 4",  360),
                                          rep("Type 5",  420),
                                          rep("Type 6",  480),
                                          rep("Type 7",  600),
                                          rep("Type 8",  960),
                                          rep("Type 9",  1200),
                                          rep("Type 10", 1500)))







Sim_trainingdata2 = iccelf(data = Sim_train_2,
                           simulated_N = 5000, 
                           balanced_N = 600,
                           N_layer = 5,
                           error = 0.0001,
                           cell_label = c(rep("Type 1",  120),
                                          rep("Type 2",  240),
                                          rep("Type 3",  360),
                                          rep("Type 4",  480),
                                          rep("Type 5",  480),
                                          rep("Type 6",  600),
                                          rep("Type 7",  600),
                                          rep("Type 8",  840),
                                          rep("Type 9",  960),
                                          rep("Type 10", 1320)))








Sim_trainingdata3 = iccelf(data = Sim_train_3,
                           simulated_N = 5000, 
                           balanced_N = 600,
                           N_layer = 5,
                           error = 0.0001,
                           cell_label = c(rep("Type 1",  240),
                                          rep("Type 2",  360),
                                          rep("Type 3",  420),
                                          rep("Type 4",  480),
                                          rep("Type 5",  480),
                                          rep("Type 6",  600),
                                          rep("Type 7",  600),
                                          rep("Type 8",  720),
                                          rep("Type 9",  900),
                                          rep("Type 10", 1200)))







Sim_trainingdata4 = iccelf(data = Sim_train_4,
                           simulated_N = 5000, 
                           balanced_N = 600,
                           N_layer = 5,
                           error = 0.0001,
                           cell_label = c(rep("Type 1",  360),
                                          rep("Type 2",  420),
                                          rep("Type 3",  480),
                                          rep("Type 4",  540),
                                          rep("Type 5",  600),
                                          rep("Type 6",  600),
                                          rep("Type 7",  660),
                                          rep("Type 8",  720),
                                          rep("Type 9",  780),
                                          rep("Type 10", 840)))









Sim_trainingdata5 = iccelf(data = Sim_train_5,
                           simulated_N = 5000, 
                           balanced_N = 600,
                           N_layer = 5,
                           error = 0.0001,
                           cell_label = c(rep("Type 1",  600),
                                          rep("Type 2",  600),
                                          rep("Type 3",  600),
                                          rep("Type 4",  600),
                                          rep("Type 5",  600),
                                          rep("Type 6",  600),
                                          rep("Type 7",  600),
                                          rep("Type 8",  600),
                                          rep("Type 9",  600),
                                          rep("Type 10", 600)))


