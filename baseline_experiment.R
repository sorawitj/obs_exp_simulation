library(grf)
library(ggplot2)

exp_df <- read.csv("data/experimental_data.csv")
obs_df <- read.csv("data/observational_data.csv")
test_df <- read.csv("test_data.csv")

sample_size = c(100, 400, 800, 1000)
rmse.obs <- c()
rmse.exp <- c()
rank_corr.obs <- c()
rank_corr.exp <- c()

get_eval_metric <- function(train_df, test_df) {
  # Train a causal forest on experimental data.
  X <- as.matrix(train_df[,-(1:2)])
  W <- as.matrix(train_df[,2])
  Y <- as.matrix(train_df[,1])
  tau.forest <- causal_forest(X, Y, W)
  
  # Estimate treatment effects for the test sample.
  tau.hat <- predict(tau.forest, as.matrix(test_df[,-(1:2)]))
  cate <- test_df[, 1]
  pred_cate <- tau.hat$predictions
  
  list(rmse=sqrt(mean((cate - pred_cate)**2)), rank_corr=cor(cate, pred_cate, method='spearman'))
}

for (n in sample_size){
  train_df.exp <- exp_df[1:n,]
  eval_metric <- get_eval_metric(train_df.exp, test_df)
  rmse.exp <- c(rmse.exp, eval_metric$rmse)
  rank_corr.exp <- c(rank_corr.exp, eval_metric$rank_corr)
  
  n.obs <- n*5
  train_df.obs <- obs_df[1:n.obs,]
  eval_metric <- get_eval_metric(train_df.obs, test_df)
  rmse.obs <- c(rmse.obs, eval_metric$rmse)
  rank_corr.obs <- c(rank_corr.obs, eval_metric$rank_corr)
}

df_res <- data.frame(rmse=c(rmse.exp, rmse.obs), rank_corr=c(rank_corr.exp, rank_corr.obs),
                     data_source=c(rep('exp', length(rmse.obs)), rep('obs', length(rmse.obs))),
                     sample_size=c(rep(sample_size, 2)))

p<-ggplot(df_res, aes(x=sample_size, y=rmse, group=data_source)) +
  geom_line(aes(color=data_source))+
  geom_point(aes(color=data_source))
aspect_ratio <- 1
p
ggsave('baseline_result.pdf',p, height = 4 , width = 4 * aspect_ratio)
