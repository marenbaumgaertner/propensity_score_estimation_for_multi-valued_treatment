library(readr)

# Import dataset
main_dataset <- read.csv2("data/Main_Dataset.csv")
main_dataset <- na.omit(main_dataset)

## Take Lat_Length = 2000ft is taken as treatment reference value 
raw_dataset <- main_dataset[,c("K_min", "K_max", "Por_min", "Por_max", "Pore_pressure", 
                               "Fracture_length_ft", "Fracture_height_ft", "Fracture_perm_md", "Fracture_width_in",
                               "Lateral_Length_ft" , "Spacing_ft",  "Spacing_efficency", "Heat_Perf_kJ")]

cov.names = c("K_min", "K_max", "Por_min", "Por_max", "Pore_pressure", 
              "Fracture_length_ft", "Fracture_height_ft", "Fracture_perm_md", "Fracture_width_in",
              "Spacing_ft", "Spacing_efficency")


## Identify covariates, treatment and outcome and preprocess data
X_raw = raw_dataset[, cov.names]
W_raw_norm = raw_dataset$Lateral_Length_ft
W_raw = (W_raw_norm - min(W_raw_norm))/1000

W.levels <- sort(unique(W_raw)); K = length(W.levels)

Y_raw_norm = raw_dataset$Heat_Perf_kJ
mean_raw = mean(log(Y_raw_norm)); sd_raw = sd(log(Y_raw_norm))
Y_raw = (log(Y_raw_norm) - mean_raw )/sd_raw

raw_df= data.frame(X_raw, W = W_raw, Y = Y_raw)
raw_df$Fracture_length_ft = (raw_df$Fracture_length_ft-min(raw_df$Fracture_length_ft))/(max(raw_df$Fracture_length_ft)-min(raw_df$Fracture_length_ft))
n_raw = nrow(X_raw)
d = ncol(X_raw)


table(raw_df$W)

colnames(raw_df)

random_subset_indices <- sample(nrow(raw_df), 5000, replace = FALSE)
X <- X_raw[random_subset_indices, ]
W <- W_raw[random_subset_indices]
raw_df <- raw_df[random_subset_indices, ]
#table(raw_df$W)
rownames(X) <- seq(1,5000)

X = X[, !(names(X) %in% c("K_min", "Por_min"))]
X = as.matrix(X)

dataset <- raw_df[, !(names(raw_df) %in% c("K_min", "Por_min"))] # remone due to Colinearity

rm(Y_raw, Y_raw_norm,main_dataset,raw_dataset,X_raw,W_raw,raw_df)




#raw_df_dt <- raw_df
