# this script just add a new list with the true propensity score for
# both of the synthetic datasets from Imbens 2016 and Linden 2015

n_rep <- 100

# first dataset
dataname <- "Imbens_2016"

sim_results <- readRDS(paste0("sim_results/simulation_results_ovo_", dataname,"_", n_rep, ".rds"))


for (i in 1:n_rep){
  dataset = sim_results$outcome[[i]]
  
  # this code is from data_Imbens_2016.R
  X1 = dataset$X1
  X2 = dataset$X2
  X3 = dataset$X3
  X4 = dataset$X4
  X5 = dataset$X5
  X6 = dataset$X6
  
  xb2 <- 0.1*(X1^2+X2+X3+X4+X5+X6)
  xb3 <- 0.1*(X1+X2^2+X3^2+X4+X5+X6)
  exb2<-exp(xb2)
  exb3<-exp(xb3)
  pi1<-1/(1+exp(xb2)+exp(xb3))
  pi2<-exp(xb2)/(1+exp(xb2)+exp(xb3))
  pi3<-exp(xb3)/(1+exp(xb2)+exp(xb3))
  pi<-cbind(pi1,pi2,pi3)
  colnames(pi) <- c("1", "2", "3")
  pi <- as.data.frame(pi)
  sim_results[['true_e']][[i]] <- pi
}

saveRDS(sim_results, paste0("sim_results/simulation_results_final_", dataname,"_", n_rep, ".rds"))


# Second dataset
dataname <- "Linden_2015"

sim_results <- readRDS(paste0("sim_results/simulation_results_ovo_", dataname,"_", n_rep, ".rds"))


for (i in 1:n_rep){
  dataset = sim_results$outcome[[i]]
  
  # this code is from data_Linden_2015.R
  X1 = dataset$X1
  X2 = dataset$X2
  
  ex1 <- exp(1.5 * (-0.2 + X1 + X2))
  ex2 <- exp(1.2 * (-0.1 + X1 + X2))
  qi <- 1 + ex1 + ex2
  
  prob0 <- 1 / qi
  prob1 <- ex1 / qi
  prob2 <- ex2 / qi
  
  pi<-cbind(prob0, prob1,prob2)
  colnames(pi) <- c("1", "2", "3")
  pi <- as.data.frame(pi)
  sim_results[['true_e']][[i]] <- pi
}

saveRDS(sim_results, paste0("sim_results/simulation_results_final_", dataname,"_", n_rep, ".rds"))
