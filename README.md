# ICCELF


1. The "ICCELF.R" file serves as the main R function file for our proposed framework ICCELF. It contains detailed information about the computational procedure
   and parameters used in the function. This "iccelf" function can generate the final data list for annotation analysis.
   
2. The "Simulation_scenarios_generation.R" is for creating 5 scenarios' simulation scRNA-seq data using the required package "SPARSim".

3. The "Generating_training_data.R" is to generate training data list using both real-world PBMC datasets and simulated datasets.

4. The "XGboost_and_evaluation.R" is how XGBoost is trained under ICCELF framework. It also includes the evaluation function for both PMBC and simulated datasets.

5. The "Other_Classifiers.R" is similar to "XGboost_and_evaluation.R", but it contains other machine learning and deep learning classifiers for training.

6. The "Validate_bal_and_comp_are_neccessary.R" is how we compare 4 scenarios to validate imbalanced and compositional correction are neccessary.


Each file serves a specific purpose within the framework and contributes to the overall analysis and evaluation of the ICCELF framework.
