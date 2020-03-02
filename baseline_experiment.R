library(grf)

exp_df <- read.csv("experimental_data.csv")
obs_df <- read.csv("observational_data.csv")


# Train a causal forest on experimental data.
W <- rep(mean(exp_df$A), nrow(exp_df))
X <- as.matrix(exp_df[,-(1:2)])
Y <- as.matrix(exp_df[,1])
tau.forest <- causal_forest(X, Y, W)
