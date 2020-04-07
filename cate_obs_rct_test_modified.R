

## you will need the following packages:
## c("randomForest", "dplyr", "plyr", "glmnet")
## you can install them via:
## install.packages(c("randomForest", "dplyr", "plyr", "glmnet"))

library(randomForest)
library(dplyr)
library(grf)
library(data.table)

low_studentperf_cont_sim <- function(n_sample, is_obs=TRUE, pi.x=0.5){
    
    setwd("~/work/obs_exp_simulation")
    
    #source("calcMetrics.r")
    # read in dataset from https://archive.ics.uci.edu/ml/datasets/student+performance
    
    dat<-read.table("data/student-por.csv",sep=";",header=TRUE)
    
    dat <- dat[sample(nrow(dat), n_sample, replace = TRUE), ]
    
    #names(dat)
    for(i in 1:ncol(dat))
    {
        if(i!=15&i!=30)
        {
            dat[,i]<-as.numeric(dat[,i])-1 
        }
    }
    
    A<-dat[,4]
    Y<-apply(dat[,31:33],1,mean)
    dat<-cbind(Y,A,dat[,-c(4,31:33)])
    
    m.or <- lm(Y ~ ., data =dat)
    
    #fitting outcome model and pedicting counterfactual outcomes  
    
    Q <- cbind(QAW = predict(m.or),
               Q0W = predict(m.or, newdata = data.frame(A = 0, dat[,-2])),
               Q1W = predict(m.or, newdata = data.frame(A = 1, dat[,-2])))
    
    # Treatment effect in the original data
    #diff(colMeans(Q[,-1]))  # 0.1734114
    
    
    # Model the Pscore
    g <- glm(A ~ ., data = dat[,-1], family = "binomial")
    g1W <- predict(g, type = "response")
    #summary(g1W)  
    
    #hist(g1W)
    
    # Simulate the outcome and treatment assignment for each experimental dataset
    # Modification 1 : very close to the logistic regression models fitted to the 
    # actual data, however, higher prevalence of the exposure
    n <- nrow(dat)
    d.mod1 <- dat
    # Generate treatment so that it is a little less problematic
    beta <- coef(g)
    
    g1W.mod1 <- as.numeric(plogis(cbind(1, as.matrix(d.mod1[,-(1:2)])) %*% beta))*0.95
    #summary(g1W.mod1)
    #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    #0.05716 0.50660 0.73960 0.66160 0.83650 0.91780
    
    # new variables
    # risk.edu<-log((rowSums(dat[,5:22])+dat$higher*3+dat$health*5)+10) # risk based on education 
    risk.edu <- 3*dat$higher
    
    if (is_obs) {
    d.mod1$A <- rbinom(n, 1, g1W.mod1)
    }
    else{
        d.mod1$A <- rbinom(n, 1, pi.x)
    }
    # Use glm model as truth for this simulation
    beta.Q.mod1 <- coef(m.or)
    beta.Q.mod1[2] <- 0
    lp.drs.mod1 <-  cbind(1, as.matrix(dat[,-1])) %*% beta.Q.mod1
    psi0.mod1<- 2
    E_Y_given_X_trt <- function(lp.drs.mod1, risk.edu, trt)
    {
        main <- lp.drs.mod1
        interact <- 2*risk.edu
        main + trt*interact + trt*psi0.mod1
    }
    d.mod1$Y <- rnorm(n,E_Y_given_X_trt(lp.drs.mod1, risk.edu, d.mod1$A),sd(m.or$residuals))
    
    #rename remaining variables  (covariates / confounders)
    names(d.mod1)[which(names(d.mod1)!="A"&names(d.mod1)!="Y")]<-paste("V",1:(length(names(d.mod1))-2),sep="")
    
    rownames(d.mod1)<-1:nrow(d.mod1)
    
    #remove covariate to emulate hidden confounding scenario
    # hidden_cov <- c('Medu', 'famsize', 'guardian', 'failures', 'G3', 'absences', 'G2', 'G1', 'goout')
    # dat <- dat[ , !(names(dat) %in% hidden_cov)]
    
    list(x = d.mod1[, -(1:2)], 
         y = d.mod1[, 'Y'], 
         trt = d.mod1[, 'A'], 
         data = data.frame(y = d.mod1[, 'Y'], x = d.mod1[, -(1:2)], trt = d.mod1[, 'A']),
         delta = E_Y_given_X_trt(lp.drs.mod1, risk.edu, 1) - E_Y_given_X_trt(lp.drs.mod1, risk.edu, 0), 
         mu_1 = E_Y_given_X_trt(lp.drs.mod1, risk.edu, 1),
         mu_0 = E_Y_given_X_trt(lp.drs.mod1, risk.edu, 0))
}

## estimate value function
value_func <- function(mu_1, mu_0, recommended.trts)
{
    value_vec <- (recommended.trts == 1) * mu_1 + (recommended.trts != 1) * mu_0
    mean(value_vec)
}

## accuracy of treatment assignments
accuracy <- function(delta.true, recommended.trts)
{
    mean( (delta.true > 0) == recommended.trts)
}


## this is the basic version of the modified outcome method
## without the use of an offset from a preliminary estimate of the CATE.
## this function does have the option of using a modified residual (Y - \hat{E}[Y|X])*q,
## (where \hat{E}[Y|X] is specified through main.hat)
## instead of just the modified outcome Y*q.
#### NOTE: I use the term 'delta' here -- this is just another notation that's commonly
####       used for the CATE, ie Delta(X) = tau(X) = E[Y(1)|X] - E[Y(-1)|X]
modified_outcome_method <- function(x, y, trt,      ## RCT data (x is a matrix of covariates, y a vector of responses, trt vector for treatment statuses)
                                    x.test, y.test, trt.test, ## test data
                                    delta.hat,      ## estimate/prediction of the CATEs from an external data source for the RCT data
                                    main.hat,       ## estimate/prediction of the main effects  E(Y|X) from an external data source for the RCT data
                                    delta.hat.test, ## estimate/prediction of the CATEs from an external data source for the test data
                                    main.hat.test,  ## estimate/prediction of the main effects  E(Y|X) from an external data source for the test data
                                    outcome.type = c("y", "residual"),      ## should we use the response directly or a residual (Y - E(Y|X))?
                                    calibration.method = c("linear", "rf"), ## how to calibrate main.hat to the RCT data (linear works well)
                                    model.type = c("linear", "rf"))         ## linear model or RF for estimating the CATE?
{
    model.type         <- match.arg(model.type)
    calibration.method <- match.arg(calibration.method)
    outcome.type       <- match.arg(outcome.type)
    
    ## saving variable names 
    
    x.vnames <- colnames(x.vnames)
    if (is.null(x.vnames))
    {
        x.vnames <- paste0("x.", 1:ncol(x))
        colnames(x) <- x.vnames
    }
    
    trt.levels <- sort(unique(trt))
    
    if (length(trt.levels) == 2)
    {
        if (all(trt.levels == c(0, 1)))
        {
            trt[trt == 0] <- -1
        }
        trt.levels <- sort(unique(trt))
    }
    
    if (!(length(trt.levels) == 2) || !all(trt.levels == c(-1, 1)))
    {
        stop("'trt' must be a vector that takes two values, 1 for treated observations and -1 for untreated observations")
    }
    
    trt.levels <- sort(unique(trt.test))
    
    if (length(trt.levels) == 2)
    {
        if (all(trt.levels == c(0, 1)))
        {
            trt.test[trt.test == 0] <- -1
        }
        trt.levels <- sort(unique(trt.test))
    }
    
    if (!(length(trt.levels) == 2) || !all(trt.levels == c(-1, 1)))
    {
        stop("'trt.test' must be a vector that takes two values, 1 for treated observations and -1 for untreated observations")
    }
    
    stopifnot(ncol(x) == ncol(x.test))
    
    colnames(x.test) <- x.vnames
    x.plus <- paste(x.vnames, collapse = "+")
    
    ### make sure dimensions match
    stopifnot(nrow(x) == NROW(y))
    stopifnot(nrow(x) == NROW(trt))
    stopifnot(nrow(x) == NROW(main.hat))
    stopifnot(nrow(x) == NROW(delta.hat))
    
    
    stopifnot(nrow(x.test) == NROW(y.test))
    stopifnot(nrow(x.test) == NROW(trt.test))
    stopifnot(nrow(x.test) == NROW(main.hat.test))
    stopifnot(nrow(x.test) == NROW(delta.hat.test))
    
    data      <- data.frame(y = y, trt = trt, x)
    data.test <- data.frame(y = y.test, trt = trt.test, x.test)
    
    
    rct_dat_prop <- data
    rct_dat_prop$trt[rct_dat_prop$trt != 1] <- 0
    
    ##  -------------------
    ##  fit a propensity score model (even though we know the trt is randomized)
    ##  -------------------
    propens.mod <- glm(as.formula(paste("trt~", x.plus)), data = rct_dat_prop,
                       family = binomial())
    
    pi.x.hat <- drop(unname(fitted(propens.mod)))

    q_weight <- rct_dat_prop$trt / pi.x.hat - (1 - rct_dat_prop$trt) / (1 - pi.x.hat)
    
    q_weight_plus <- rct_dat_prop$trt / pi.x.hat + (1 - rct_dat_prop$trt) / (1 - pi.x.hat)
    
    q_weight <- q_weight / mean(q_weight_plus)
    
    if (calibration.method == "rf")
    {
        #data.rf <- data.frame(y = rct_dat_prop$y, main = main.hat, data[,x.vnames])
        data.rf <- data.frame(y = rct_dat_prop$y, main = main.hat, trt = data$trt)
        rf.main <- randomForest(y ~ ., data = data.rf, ntree = 2000, mtry = length(x.vnames) / 2)
        plot(rf.main)
        varImpPlot(rf.main)
        #main.hat.cal <- predict(rf.main, data.rf)
        
        data.rf1 <- data.rf0 <- data.rf
        data.rf1$trt <- 1
        data.rf0$trt <- -1
        main.hat.cal <- 0.5 * (unname(predict(rf.main, data.rf1)) + unname(predict(rf.main, data.rf0)))
    } else
    {
        data.rf <- data.frame(y = rct_dat_prop$y, main = main.hat, trt = data$trt)
        rf.main <- lm(y ~ .*trt, data = data.rf)
        # print(summary(rf.main))   
        
        data.rf1 <- data.rf0 <- data.rf
        data.rf1$trt <- 1
        data.rf0$trt <- -1
        main.hat.cal <- 0.5 * (unname(predict(rf.main, data.rf1)) + unname(predict(rf.main, data.rf0)))
    }
    
    rct_dat_prop$y.mod     <- rct_dat_prop$y * q_weight
    rct_dat_prop$resid.mod <- (rct_dat_prop$y - main.hat.cal) * q_weight
    rct_dat_prop$main.hat  <- main.hat
    
    if (outcome.type == "y")
    {
        formul_mod_resid_only_model <- as.formula(paste("y.mod ~", x.plus))
    } else
    {
        formul_mod_resid_only_model <- as.formula(paste("resid.mod ~", x.plus))
    }
    
    if (model.type == "linear")
    {
        mod_resid_only_model <- lm(formul_mod_resid_only_model, data = rct_dat_prop)
        
        delta_test <- unname(drop(predict(mod_resid_only_model, newdata = data.test)))
    } else
    {
        mod_resid_only_model <- randomForest(formul_mod_resid_only_model, data = rct_dat_prop,
                                             ntree = 1000)
        
        delta_test <- unname(drop(predict(mod_resid_only_model, newdata = data.test)))
    }
    
    list(model = mod_resid_only_model, delta_test = delta_test)
}


## this is the modified outcome method, where the estimate can either 
## be augmented with a preliminary estimate of the CATE via either:
##   i.  as an offset, in which case we need to add the preliminary estimate of the CATE back at the end
##   ii. as a predictor, where we would treate the prelim CATE est as a basis function of X
## this function has the option of using a modified residual (Y - \hat{E}[Y|X])*q,
## (where \hat{E}[Y|X] is specified through main.hat)
## instead of just the modified outcome Y*q
modified_outcome_aug_method <- function(x, y, trt,  ## RCT data (x is a matrix of covariates, y a vector of responses, trt vector for treatment statuses)
                                        x.test, y.test, trt.test, ## test data
                                        delta.hat, main.hat,
                                        delta.hat.test, main.hat.test,
                                        mu.1.hat, mu.1.hat.test,
                                        outcome.type = c("y", "residual"),
                                        aug.type = c("linear", "smooth", "rf", "offset.rf", "offset.linear"),
                                        calibration.method = c("linear", "rf"))
{
    outcome.type       <- match.arg(outcome.type)
    aug.type           <- match.arg(aug.type)
    calibration.method <- match.arg(calibration.method)
    
    ## saving variable names 
    
    x.vnames <- colnames(x.vnames)
    if (is.null(x.vnames))
    {
        x.vnames <- paste0("x.", 1:ncol(x))
        colnames(x) <- x.vnames
    }
    
    trt.levels <- sort(unique(trt))
    
    if (length(trt.levels) == 2)
    {
        if (all(trt.levels == c(0, 1)))
        {
            trt[trt == 0] <- -1
        }
        trt.levels <- sort(unique(trt))
    }
    
    if (!(length(trt.levels) == 2) || !all(trt.levels == c(-1, 1)))
    {
        stop("'trt' must be a vector that takes two values, 1 for treated observations and -1 for untreated observations")
    }
    
    trt.levels <- sort(unique(trt.test))
    
    if (length(trt.levels) == 2)
    {
        if (all(trt.levels == c(0, 1)))
        {
            trt.test[trt.test == 0] <- -1
        }
        trt.levels <- sort(unique(trt.test))
    }
    
    if (!(length(trt.levels) == 2) || !all(trt.levels == c(-1, 1)))
    {
        stop("'trt.test' must be a vector that takes two values, 1 for treated observations and -1 for untreated observations")
    }
    
    stopifnot(ncol(x) == ncol(x.test))
    
    colnames(x.test) <- x.vnames
    x.plus <- paste(x.vnames, collapse = "+")
    
    ### make sure dimensions match
    stopifnot(nrow(x) == NROW(y))
    stopifnot(nrow(x) == NROW(trt))
    stopifnot(nrow(x) == NROW(main.hat))
    stopifnot(nrow(x) == NROW(delta.hat))
    
    stopifnot(nrow(x.test) == NROW(y.test))
    stopifnot(nrow(x.test) == NROW(trt.test))
    stopifnot(nrow(x.test) == NROW(main.hat.test))
    stopifnot(nrow(x.test) == NROW(delta.hat.test))
    
    data      <- data.frame(y = y, trt = trt, x)
    data.test <- data.frame(y = y.test, trt = trt.test, x.test)
    
    
    
    ### fit propensity score model 
    rct_dat_prop <- data
    rct_dat_prop$trt[rct_dat_prop$trt != 1] <- 0
    propens.mod <- glm(as.formula(paste("trt~", x.plus)), data = rct_dat_prop,
                       family = binomial())
    
    pi.x.hat <- drop(unname(fitted(propens.mod)))
    
    ##  -------------------
    ##  calculate modified outcome
    ##  -------------------
    
    q_weight <- rct_dat_prop$trt / pi.x.hat - (1 - rct_dat_prop$trt) / (1 - pi.x.hat)
    
    ## calculate a normalizing factor (sort of like stabilized weights, ie a Hajek estimator vs a H-T estimator)
    q_weight_plus <- rct_dat_prop$trt / pi.x.hat + (1 - rct_dat_prop$trt) / (1 - pi.x.hat)
    
    q_weight <- q_weight / mean(q_weight_plus)
    
    
    
    ##  -------------------
    ##  calibrate the predicted main effects (E[Y|X])
    ##  to the RCT data
    ##  -------------------

    
    if (calibration.method == "rf")
    {
        #data.rf <- data.frame(y = rct_dat_prop$y, main = main.hat, data[,x.vnames])
        data.rf <- data.frame(y = rct_dat_prop$y, main = main.hat, trt = data$trt)
        rf.main <- randomForest(y ~ ., data = data.rf, ntree = 1000)
        plot(rf.main)
        varImpPlot(rf.main)
        
        data.rf1 <- data.rf0 <- data.rf
        data.rf1$trt <- 1
        data.rf0$trt <- -1
        main.hat.cal <- 0.5 * (unname(predict(rf.main, data.rf1)) + unname(predict(rf.main, data.rf0)))
    } else
    {
        data.rf <- data.frame(y = rct_dat_prop$y, main = main.hat, trt = data$trt)
        rf.main <- lm(y ~ .*trt, data = data.rf)
        # print(summary(rf.main))   
        
        data.rf1 <- data.rf0 <- data.rf
        data.rf1$trt <- 1
        data.rf0$trt <- -1
        main.hat.cal <- 0.5 * (unname(predict(rf.main, data.rf1)) + unname(predict(rf.main, data.rf0)))
    }
    
    
    
    rct_dat_prop$y.mod     <- rct_dat_prop$y * q_weight
    
    if (aug.type == "offset")
    {
        rct_dat_prop$resid.mod <- (rct_dat_prop$y - main.hat.cal) * q_weight - delta.hat
    } else
    {
        rct_dat_prop$resid.mod <- (rct_dat_prop$y - main.hat.cal) * q_weight
    }
    
    rct_dat_prop$main.hat  <- main.hat
    rct_dat_prop$delta.hat <- delta.hat
    
    data.test$delta.hat    <- delta.hat.test
    data.test$mu.1.hat     <- mu.1.hat.test
    
    ##  -------------------
    ##  estimate CATE depending on model
    ##  specifications
    ##  -------------------
    
    if (aug.type == "linear")
    {
        if (outcome.type == "y")
        {
            formul_mod_resid_model <- as.formula(paste("y.mod ~", x.plus, "+delta.hat+mu.1.hat"))
        } else
        {
            formul_mod_resid_model <- as.formula(paste("resid.mod ~", x.plus, "+delta.hat+mu.1.hat"))   
        }
        
        mod_resid_aug_model <- lm(formul_mod_resid_model, data = rct_dat_prop)
    } else if (aug.type == "smooth")
    {
        if (outcome.type == "y")
        {
            formul_mod_resid_model_gam <- as.formula(paste("y.mod ~", x.plus, "+s(delta.hat,mu.1.hat)"))
        } else
        {
            formul_mod_resid_model_gam <- as.formula(paste("resid.mod ~", x.plus, "+s(delta.hat,mu.1.hat)"))
        }
        
        
        mod_resid_aug_model <- gam(formul_mod_resid_model_gam, data = rct_dat_prop, family = gaussian())
    } else if (aug.type == "rf")
    {
        if (outcome.type == "y")
        {
            formul_mod_resid_model <- as.formula(paste("y.mod ~", x.plus, "+delta.hat+mu.1.hat"))
        } else
        {
            formul_mod_resid_model <- as.formula(paste("resid.mod ~", x.plus, "+delta.hat+mu.1.hat"))   
        }
        mod_resid_aug_model <- randomForest(formul_mod_resid_model, data = rct_dat_prop, ntree = 2000,
                                            mtry = length(x.vnames) / 2)
    } else if (aug.type == "offset.rf")
    {
        if (outcome.type == "y")
        {
            formul_mod_resid_model <- as.formula(paste("y.mod ~", x.plus))
        } else
        {
            formul_mod_resid_model <- as.formula(paste("resid.mod ~", x.plus))   
        }
        mod_resid_aug_model <- randomForest(formul_mod_resid_model, data = rct_dat_prop, ntree = 2000,
                                            mtry = length(x.vnames) / 2)
    } else if (aug.type == "offset.linear")
    {
        if (outcome.type == "y")
        {
            formul_mod_resid_model <- as.formula(paste("y.mod ~", x.plus))
        } else
        {
            formul_mod_resid_model <- as.formula(paste("resid.mod ~", x.plus))   
        }
        mod_resid_aug_model <- lm(formul_mod_resid_model, data = rct_dat_prop)
    }
    
    
    ## predicted CATE for test data
    delta_test <- unname(drop(predict(mod_resid_aug_model, newdata = data.test)))
    
    if (aug.type == "offset")
    {
        delta_test <- delta_test + delta.hat.test
    }
    
    list(model = mod_resid_aug_model, delta_test = delta_test)
}


## this function has the same arguments as `modified_outcome_aug_method()`
estimate_cate <- function(x, y, trt,
                          x.test, y.test, trt.test, ## test data
                          delta.hat, main.hat,
                          delta.hat.test, main.hat.test,
                          mu.1.hat, mu.1.hat.test,
                          outcome_type = c("y", "residual"),           ### outcome type. residual uses calibrated main effect estimate from Obs data
                          cate_functional_form = c("linear", "rf"),    ### functional form for the estimated CATE
                          prelim_cate_use = c("none", "predictor", "offset", "smooth.predictor")) ### how should the preliminary CATE be used?
{
    cate_functional_form <- match.arg(cate_functional_form)
    prelim_cate_use      <- match.arg(prelim_cate_use)
    outcome_type         <- match.arg(outcome_type)
    
    
    
    
    
    if (prelim_cate_use == "none")
    {
        result <- modified_outcome_method(x = x, y = y, trt = trt, 
                                          x.test = x.test, y.test = y.test, trt.test = trt.test, 
                                          delta.hat = delta.hat, main.hat = main.hat,
                                          delta.hat.test = delta.hat.test, main.hat.test = main.hat.test,
                                          model.type = cate_functional_form,
                                          outcome.type = outcome_type, calibration.method = "linear")
    } else
    {
        if (prelim_cate_use == "predictor")
        {
            aug.type <- cate_functional_form
        } else if (prelim_cate_use == "smooth.predictor")
        {
            aug.type <- "smooth"
        } else
        {
            aug.type <- paste(prelim_cate_use, cate_functional_form, sep=".")
        }
        
        result <- modified_outcome_aug_method(x = x, y = y, trt = trt, 
                                              x.test = x.test, y.test = y.test, trt.test = trt.test,
                                              delta.hat = delta.hat, main.hat = main.hat,
                                              delta.hat.test = delta.hat.test, main.hat.test = main.hat.test,
                                              mu.1.hat = mu.1.hat, mu.1.hat.test = mu.1.hat.test,
                                              outcome.type = outcome_type, 
                                              aug.type = aug.type, calibration.method = "linear")
    }
    result
}

## just a wrapper function to take arguments in a more convenient way for looping
estimate_cate_grid <- function(grid, 
                               x, y, trt,
                               x.test, y.test, trt.test, ## test data
                               delta.hat, main.hat,
                               delta.hat.test, main.hat.test,
                               mu.1.hat, mu.1.hat.test)
{
    print(grid)
    estimate_cate(x = x, y = y, trt = trt,
                  x.test = x.test, y.test = y.test, trt.test = trt.test, ## test data
                  delta.hat = delta.hat, main.hat = main.hat,
                  delta.hat.test = delta.hat.test, main.hat.test = main.hat.test,
                  mu.1.hat = mu.1.hat, mu.1.hat.test = mu.1.hat.test,
                  outcome_type = grid$outcome_type,
                  cate_functional_form = grid$cate_functional_form,
                  prelim_cate_use = grid$prelim_cate_use)
}

## this is a data frame with an enumeration of all of the 
## possible modified outcome regression methods for estimating the CATE.
## each row is a different modeling setup
model_setup_grid <- expand.grid(outcome_type = c("y", "residual"),  ### outcome type. residual uses calibrated main effect estimate from Obs data
                                cate_functional_form = c("linear", "rf"), ### functional form for the estimated CATE. linear or random forest
                                prelim_cate_use = c("none", "predictor", "offset"), ### how should the preliminary CATE estimate be used?
                                stringsAsFactors = FALSE)

## prelim_cate_use = "predictor" uses the preliminary CATE estimate as a new predictor (ie some basis function of X)
## prelim_cate_use = "offset" uses the preliminary CATE estimate as an offset in the model for the CATE (like the Kallus approach). 
## The final CATE estimate when using ``offset'' includes the offset as an additive term, like in the Kallus paper

model_setup_grid


## trial sample size
n.rct <- 250

## observational data sample size
n.obs <- 2500

## test set size for independent sample from trial population
n.rct.test <- 10000

## test set size for independent sample from observational study population
n.obs.test <- 10000

set.seed(42)

## generate training data from trial population
rct_dat <- low_studentperf_cont_sim(n.rct, FALSE)

## generate test data from trial population
rct_test <- low_studentperf_cont_sim(n.rct.test, FALSE)

## generate training data from observational study population
obs_dat <- low_studentperf_cont_sim(n.obs, TRUE)

## generate test data from observational study population
obs_test <- low_studentperf_cont_sim(n.obs.test, TRUE)

## proportion treated in the observational study:
mean(obs_test$trt)


##  -------------------
##  fit a random forest regression model to the observational
##  study population mean outcomes, ie a model for
##  E[Y|X,Trt,Study = Obs]
##  -------------------

x.vnames <- names(obs_dat$data)[grepl("x.", names(obs_dat$data))]
x.plus <- paste(x.vnames, collapse = "+")

outcome.rf.formula <- as.formula(paste("y ~", x.plus, "+trt"))

obs.outcome.rf.mod <- randomForest(outcome.rf.formula, data = obs_dat$data, ntree = 1000, do.trace = 100)
# obs.outcome.rf.mod <- glm(outcome.rf.formula, data = obs_dat$data, family = "gaussian")

plot(obs.outcome.rf.mod)

varImpPlot(obs.outcome.rf.mod)

cforest_obs <- causal_forest(obs_dat$x, obs_dat$y, obs_dat$trt)


##  -------------------
##  save predictions for RCT data
##  -------------------

rct_dat_1 <- rct_dat_0 <- rct_dat$data
rct_dat_1$trt <- 1
rct_dat_0$trt <- 0

preds_1 <- unname(predict(obs.outcome.rf.mod, newdata = rct_dat_1))
preds_0 <- unname(predict(obs.outcome.rf.mod, newdata = rct_dat_0))

## preliminary estimate of the CATE
delta.hat <- preds_1 - preds_0
# delta.hat <- predict(cforest_obs, rct_dat$x)[,1]

## preliminary estimate of the main effects 
main.hat <- 0.5 * (preds_1 + preds_0)


##  -------------------
##  save predictions for RCT test dataset
##  -------------------

rct_dat_test_1 <- rct_dat_test_0 <- rct_test$data
rct_dat_test_1$trt <- 1
rct_dat_test_0$trt <- 0

preds_test_rct_1 <- unname(predict(obs.outcome.rf.mod, newdata = rct_dat_test_1))
preds_test_rct_0 <- unname(predict(obs.outcome.rf.mod, newdata = rct_dat_test_0))

## preliminary estimate of the CATE
delta.hat.rct.test <- preds_test_rct_1 - preds_test_rct_0
# delta.hat.rct.test <- predict(cforest_obs, rct_test$x)[,1]

main.hat.rct.test <- 0.5 * (preds_test_rct_1 + preds_test_rct_0)

# ##  -------------------
# ##  save predictions for observational study test dataset
# ##  -------------------
# 
# obs_dat_test_1 <- obs_dat_test_0 <- obs_test$data
# obs_dat_test_1$trt <- 1
# obs_dat_test_0$trt <- -1
# 
# preds_test_obs_1 <- unname(predict(obs.outcome.rf.mod, newdata = obs_dat_test_1))
# preds_test_obs_0 <- unname(predict(obs.outcome.rf.mod, newdata = obs_dat_test_0))
# 
# ## preliminary estimate of the CATE
# delta.hat.obs.test <- preds_test_obs_1 - preds_test_obs_0
# # delta.hat.obs.test <- predict(cforest_obs, obs_test$x)[,1]
# 
# ## preliminary estimate of the main effects 
# main.hat.obs.test <- 0.5 * (preds_test_obs_1 + preds_test_obs_0)
# 
# ###


##  -------------------
##  estimate CATE with all methods
##  and save predictions for the independent RCT dataset
##  -------------------

result_list_rct <- plyr::alply(model_setup_grid, .margin = 1, .fun = estimate_cate_grid, 
                               x = rct_dat$x, y = rct_dat$y, trt = rct_dat$trt, 
                               x.test = rct_test$x, y.test = rct_test$y, trt.test = rct_test$trt, 
                               delta.hat = delta.hat, main.hat = main.hat,
                               delta.hat.test = delta.hat.rct.test, main.hat.test = main.hat.rct.test,
                               mu.1.hat = preds_1, mu.1.hat.test = preds_test_rct_1)

## evaluate performance on the test dataset.
## `RankCorrelation` is the rank correlation of the estimated
## CATE with the true CATE.
## 
## `Value` is the estimated value function under the treatment rule
## generated by the estimated CATE, ie I(CATE > 0) for treatment assignments
## 
## `AUC` is the area under the ROC curve to evaluate the estimated CATE
## in its ability to predict I(E[Y(1)|X = x] > E[Y(0)|X = x]) for new observations
results_rct <- plyr::ldply(result_list_rct, 
                           .fun = function(ml) c(RankCorrelation = cor(ml$delta_test, rct_test$delta, method = "spearman"),
                                                 Value = value_func(rct_test$mu_1, rct_test$mu_0, ml$delta_test > 0),
                                                 AUC = glmnet:::auc(rct_test$delta > 0, ml$delta_test),
                                                 rmse = sqrt(mean((ml$delta_test - rct_test$delta)^2))))

results_rct %>% arrange(desc(-rmse))


# ##  -------------------
# ##  estimate CATE with all methods
# ##  and save predictions for the independent observational dataset
# ##  -------------------
# 
# result_list_obs <- plyr::alply(model_setup_grid, .margin = 1, .fun = estimate_cate_grid, 
#                                x = rct_dat$x, y = rct_dat$y, trt = rct_dat$trt, 
#                                x.test = obs_test$x, y.test = obs_test$y, trt.test = obs_test$trt, 
#                                delta.hat = delta.hat, main.hat = main.hat,
#                                delta.hat.test = delta.hat.obs.test, main.hat.test = main.hat.obs.test,
#                                mu.1.hat = preds_1, mu.1.hat.test = preds_test_obs_1)
# 
# results_obs <- plyr::ldply(result_list_obs, 
#                            .fun = function(ml) c(RankCorrelation = cor(ml$delta_test, obs_test$delta, method = "spearman"),
#                                                  Value = value_func(obs_test$mu_1, obs_test$mu_0, ml$delta_test > 0),
#                                                  AUC = glmnet:::auc(obs_test$delta > 0, ml$delta_test),
#                                                  rmse = sqrt(mean((ml$delta_test - rct_test$delta)^2))))
# 
# results_obs %>% arrange(desc(-rmse))

cforest_rct <- causal_forest(rct_dat$x, rct_dat$y, rct_dat$trt)
cforest_obs <- causal_forest(obs_dat$x, obs_dat$y, obs_dat$trt)
delta.cforest.rct <- predict(cforest_rct, rct_test$x)
delta.cforest.obs <- predict(cforest_obs, rct_test$x)


rmse.rct.cf <- sqrt(colMeans((delta.cforest.rct - rct_test$delta)^2))[1]
rmse.obs.cf <- sqrt(colMeans((delta.cforest.obs - rct_test$delta)^2))[1]

cor(delta.cforest.rct[,1], as.vector(rct_test$delta), method='spearman')
cor(delta.cforest.obs[,1], as.vector(rct_test$delta), method='spearman')

print(rmse.rct.cf)
print(rmse.obs.cf)

rmse.obs <- sqrt(mean((delta.hat.rct.test - rct_test$delta)^2))
print(rmse.obs)
