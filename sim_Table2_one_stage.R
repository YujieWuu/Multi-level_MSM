#########################################################################
#########################################################################
#### 1. We only consider correlation within testing site (one layer),
#### 2. For the MSM, we use the stabilized weights, and we are using
####     the estimated weights 
#### 3. We use the geeglm to do the fitting


#### load library and useful functions
library(mvtnorm)
library(geepack)
library(lme4)
expit <- function(x){
  return( exp(x)/(1+exp(x)) )
}

###################################################
#### generate testing site for each individual ####
###################################################
gen.site <- function(M, size){
  #### param@M:    Number of testing sites.
  #### param@size: Total number of study participants.
  site <- rep(x= (1:M), each = (size/M))
  return(site)
}

###################################################
######### generate X0 for each individual #########
###################################################
gen.X0 <- function(mu0, sigma0,
                   size){
  #### param@mu0:    mean for X0
  #### param@sigma0: sd for X0
  #### param@size:   Total number of study participants
  X0 <- rnorm(size, mean = mu0, sd = sigma0)
  return(X0)
}

###################################################
######### generate A0 for each individual #########
###################################################
gen.A0 <- function(X0, 
                   beta0.0, beta0.1,
                   size){
  #### E(A0|X0)=expit( beta0.0 + beta0.1*X0 )
  prob <- expit( beta0.0 + beta0.1*X0 )
  A0   <- rbinom(size, 1, prob)
  return(list(A0 = A0, prob = prob))
}

###################################################
######### generate X1 for each individual #########
###################################################
gen.X1 <- function(X0, A0, site, 
                   gamma1.0, gamma1.1, gamma1.2, gamma1.site,
                   sigma.X1, size){
  #### X1 = gamma1.0 + gamma1.1*X0 + gamma1.2*A0 + gamma1.site*I(Site=s) 
  ####        + N(0, sigma.X1^2)
  mu.X1 <- gamma1.0 + gamma1.1*X0 + gamma1.2*A0 + gamma1.site[as.numeric(site)]
  X1    <- rnorm(size, mean = mu.X1, sd = sigma.X1)
  return(X1)
}

###################################################
######### generate Y1 for each individual #########
###################################################
#### Notice that Y_1=[Y_{1,l=0}, Y_{1,l=1}]^T
gen.Y1 <- function(X0, A0, site,
                   alpha0, alpha1.1, alpha2.1, alpha.site,
                   Sigma.Y1, size){
  #### E(Y1_|...) = alpha0 + alpha1.1*X0 + alpha2.1*A0 + alpha.site*I(Site=site)
  #### Sigma.Y1 is the variance-covariance matrix of Y1
  mu.Y1 <- matrix(alpha0 + alpha1.1*X0 + alpha2.1*A0 + alpha.site[site], ncol = 1)
  mu.Y1 <- cbind(mu.Y1, mu.Y1)
  Y1    <- matrix(NA, nrow = size, ncol = 2)
  for(i in 1:size){
    Y1[i,] <- rmvnorm(1, mean = mu.Y1[i,], sigma = Sigma.Y1)
  }
  return(Y1)
}

###################################################
######### generate A1 for each individual #########
###################################################
gen.A1 <- function(X0, A0, X1, Y1, site,
                   beta1.0, beta1.1, beta1.2, beta1.3, beta1.4, beta1.site,
                   size){
  #### E(A1|...) = expit( beta1.0 + beta1.1*X0 + beta1.2*X1 + beta1.3*mean(Y1) +
  ####               beta1.4*A0 + beta1.site*I(Site=site) )
  prob <- expit( beta1.0 + beta1.1*X0 + beta1.2*X1 + beta1.3*rowMeans(Y1) +
                   beta1.4*A0 + beta1.site[site] )
  A1   <- rbinom(size, 1, prob)
  return(list(A1 = A1, prob  = prob))
}

###################################################
######### generate X2 for each individual #########
###################################################
gen.X2 <- function(X0, A0, X1, A1, site,
                   gamma2.0, gamma2.1, gamma2.2, gamma2.3, gamma2.4, gamma2.site,
                   sigma.X2, size){
  #### X2 = gamma2.0 + gamma2.1*X0 + gamma2.2*A0 + gamma2.3*X1 + gamma2.4*A1 + 
  ####              + gamma2.site*I(Site=site) + N(0, sigma.X2^2)
  mu <- gamma2.0 + gamma2.1*X0 + gamma2.2*A0 + gamma2.3*X1 + gamma2.4*A1 + gamma2.site[site]
  X2 <- rnorm(size, mu, sigma.X2)
  return(X2)
}

###################################################
######### generate A2 for each individual #########
###################################################
gen.A2 <- function(X0, A0, X1, Y1, A1, X2, site,
                   beta2.0, beta2.1, beta2.2, beta2.3, 
                   beta2.4, beta2.5, beta2.6, beta2.site,
                   size){
  #### E(A2|...) = expit(beta2.0 + beta2.1*X0 + beta2.2*X1 + beta2.3*X2 +
  ####                beta2.4*A0 + beta2.5*A1 + beta2.6*rowMeans(Y1) + beta2.site*I(site=site))
  prob <- expit(beta2.0 + beta2.1*X0 + beta2.2*X1 + beta2.3*X2 +
                  beta2.4*A0 + beta2.5*A1 + beta2.6*rowMeans(Y1) + 
                  beta2.site[site] )
  A2   <- rbinom(size, 1, prob)
  return(list(A2 = A2, prob = prob))
}

###################################################
######### generate Y3 for each individual #########
###################################################
gen.Y3 <- function(X0, A0, X1, Y1, A1, X2, A2, site,
                   alpha0, alpha1.2, alpha1.3, alpha1.4, 
                   alpha2.2, alpha2.3, alpha2.4, 
                   alpha4, alpha.site, Sigma.Y3, size){
  #### E(Y3|...) = alpha0 + alpha1.2*X0 + alpha1.3*X1 + alpha1.4*X2 +
  ####             alpha2.2*A0 + alpha2.3*A1 + alpha2.4*A2 + 
  ####             alpha3.site*I(Site=site) + alpha4*mean(Y1)
  mu.Y3 <- matrix( alpha0 + alpha1.2*X0 + alpha1.3*X1 + alpha1.4*X2 +
                     alpha2.2*A0 + alpha2.3*A1 + alpha2.4*A2  + alpha.site[site] + alpha4*rowMeans(Y1), ncol = 1)
  mu.Y3 <- cbind(mu.Y3, mu.Y3)
  Y3    <- matrix(NA, nrow = size, ncol = 2)
  for(i in 1:size){
    Y3[i,] <- rmvnorm(1, mean = mu.Y3[i,], sigma = Sigma.Y3)
  }
  return(Y3)
}

###################################################
######### True causal parameters in MSM; ##########
######### G-formula for calculateing Y^{\bar{a}} ##
###################################################
MSM.param <- function(mu0,      sigma0,
                      beta0.0,  beta0.1,
                      gamma1.0, gamma1.1,  gamma1.2, gamma1.site,
                      alpha0,   alpha1.1,  alpha2.1, alpha2.2, alpha2.3, alpha2.4,
                      alpha.site,
                      beta1.0,  beta1.1,   beta1.2,  beta1.3,     beta1.4,  beta1.site,
                      gamma2.0, gamma2.1,  gamma2.2, gamma2.3,    gamma2.4, gamma2.site,
                      beta2.0,  beta2.1,   beta2.2,  beta2.3,     beta2.4,  beta2.5,     beta2.6, beta2.site,
                      alpha1.2,  alpha1.3, alpha1.4, alpha4){
  delta2 <- alpha1.3*gamma1.2 + alpha1.4*gamma2.2 + alpha1.4*gamma2.3*gamma1.2 + alpha2.2 + alpha4*alpha2.1
  delta3 <- alpha1.4*gamma2.4 + alpha2.3
  delta4 <- alpha2.4
  delta5 = delta6 = delta7 = delta8 <- 0
  delta9 <- alpha2.1
  delta0 <- -(-2*alpha0 + (alpha1.2 - 3*alpha1.1)*mu0 +alpha1.3*gamma1.0 + alpha1.3*gamma1.1*mu0 +
                alpha1.3*mean(gamma1.site) + alpha1.4*gamma2.0 + alpha1.4*gamma2.1*mu0 +
                alpha1.4*gamma2.3*( gamma1.0 + gamma1.1*mu0 + mean(gamma1.site) ) + 
                alpha1.4*mean(gamma2.site) - 2*mean(alpha.site) + 
                alpha4*( alpha0 + alpha1.1*mu0 + mean(alpha.site) ) )/2
  delta1 <- ( (alpha1.2 - alpha1.1)*mu0 +alpha1.3*gamma1.0 + alpha1.3*gamma1.1*mu0 +
                alpha1.3*mean(gamma1.site) + alpha1.4*gamma2.0 + alpha1.4*gamma2.1*mu0 +
                alpha1.4*gamma2.3*( gamma1.0 + gamma1.1*mu0 + mean(gamma1.site) ) + 
                alpha1.4*mean(gamma2.site) + 
                alpha4*( alpha0 + alpha1.1*mu0 + mean(alpha.site) ) )/2
  re      <- data.frame(delta0 = delta0, delta1 = delta1,
                        delta2 = delta2, delta3 = delta3,
                        delta4 = delta4, delta5 = delta5,
                        delta6 = delta6, delta7 = delta7,
                        delta8 = delta8, delta9 = delta9)
  return(re)
}
######################################################################
############### Use the true weights to fit the MSM ##################
##### Preliminary try, to see if the simulation setup is correct #####
######################################################################
sim.MSM.true.weight <- function(M=30, N_m=30,
                                mu0,      sigma0,
                                beta0.0,  beta0.1,
                                gamma1.0, gamma1.1,  gamma1.2, gamma1.site, sigma.X1,
                                alpha0,   alpha1.1,  alpha2.1, alpha2.2, alpha2.3, alpha2.4,
                                alpha.site,  Sigma.Y1,
                                beta1.0,  beta1.1,   beta1.2,  beta1.3,     beta1.4,  beta1.site,
                                gamma2.0, gamma2.1,  gamma2.2, gamma2.3,    gamma2.4, gamma2.site, sigma.X2,
                                beta2.0,  beta2.1,   beta2.2,  beta2.3,     beta2.4,  beta2.5,     beta2.6,  beta2.site,
                                alpha1.2, alpha1.3,  alpha1.4, alpha4,      Sigma.Y3){
  #### generate testing site 
  site    <- gen.site(M = M, size = M*N_m)
  #### generate X0
  X0      <- gen.X0(mu0 = mu0, sigma0 = sigma0, size = M*N_m)
  #### generate A0, and get the true weights
  temp    <- gen.A0(X0  = X0, beta0.0 = beta0.0, beta0.1 = beta0.1, size = M*N_m)
  A0      <- temp$A0
  A0.prob <- temp$prob
  #### generate X1
  X1 <- gen.X1(X0 = X0,      A0 = A0,      site = site, 
               gamma1.0 = gamma1.0,    gamma1.1 = gamma1.1, 
               gamma1.2 = gamma1.2, gamma1.site = gamma1.site,
               sigma.X1 = sigma.X1,        size = M*N_m)
  #### generate Y1 for both ears
  Y1 <- gen.Y1(X0  =  X0,    A0 = A0,   site = site,
               alpha0   = alpha0,   alpha1.1 = alpha1.1, 
               alpha2.1   = alpha2.1, alpha.site = alpha.site,
               Sigma.Y1 = Sigma.Y1,     size = M*N_m)
  #### generate A1 and true weights
  temp2 <- gen.A1(X0 = X0, A0 = A0, X1 = X1, Y1 = Y1, site= site,
                  beta1.0 = beta1.0, beta1.1 = beta1.1, 
                  beta1.2 = beta1.2, beta1.3 = beta1.3, 
                  beta1.4 = beta1.4, beta1.site = beta1.site,
                  size    = M*N_m)
  A1      <- temp2$A1
  A1.prob <- temp2$prob
  #### generate X2
  X2 <- gen.X2(      X0 = X0,   A0 = A0,    X1 = X1,  A1 = A1,  site = site,
                     gamma2.0 = gamma2.0,    gamma2.1 = gamma2.1, 
                     gamma2.2 = gamma2.2,    gamma2.3 = gamma2.3, 
                     gamma2.4 = gamma2.4, gamma2.site = gamma2.site,
                     sigma.X2 = sigma.X2,        size = M*N_m)
  #### generate A2
  temp3 <- gen.A2(X0 = X0, A0 = A0, X1 = X1, Y1 = Y1, A1 = A1, X2 = X2, site = site,
                  beta2.0 = beta2.0, beta2.1 = beta2.1, beta2.2 = beta2.2, 
                  beta2.3 = beta2.3, beta2.4 = beta2.4, beta2.5 = beta2.5, 
                  beta2.6 = beta2.6, beta2.site = beta2.site, size = M*N_m)
  A2      <- temp3$A2
  A2.prob <- temp3$prob
  #### generate Y3
  Y3 <- gen.Y3(      X0 = X0,   A0 = A0,     X1 = X1,   Y1 = Y1, 
                     A1 = A1,   X2 = X2,     A2 = A2, site = site,
                     alpha0   = alpha0,      alpha1.2 = alpha1.2, 
                     alpha1.3 = alpha1.3,    alpha1.4 = alpha1.4, 
                     alpha2.2   = alpha2.2, alpha2.3   = alpha2.3, alpha2.4   = alpha2.4,
                     alpha4   = alpha4,    alpha.site = alpha.site,   
                     Sigma.Y3 = Sigma.Y3,    size     = M*N_m)
  #### create data frame for weighted gee
  #### each individual will have 4 observations
  id.df      <- c()
  time.df    <- c()
  site.df    <- c()
  X0.df      <- c()
  A0.df      <- c()
  A0.prob.df <- c()
  X1.df      <- c()
  A1.df      <- c()
  A1.prob.df <- c()
  X2.df      <- c()
  A2.df      <- c()
  A2.prob.df <- c()
  Y.df       <- c()
  #weight.df  <- c()
  #weight.whole.df <- c()
  for(i in 1:(M*N_m)){
    id.df       <- append(id.df,   rep(i, 4))
    time.df     <- append(time.df, c(rep(1,2), rep(3,2)))
    site.df     <- append(site.df, rep(site[i], 4))
    X0.df       <- append(X0.df,   rep(X0[i],4))
    A0.df       <- append(A0.df,   rep(A0[i],4))
    A0.prob.df  <- append(A0.prob.df, rep(A0.prob[i], 4))
    X1.df       <- append(X1.df,   rep(X1[i],4))
    A1.df       <- append(A1.df,   rep(A1[i],4))
    A1.prob.df  <- append(A1.prob.df, rep(A1.prob[i], 4))
    X2.df       <- append(X2.df,   rep(X2[i],4))
    A2.df       <- append(A2.df,   rep(A2[i],4))
    A2.prob.df  <- append(A2.prob.df, rep(A2.prob[i], 4))
    Y.df        <- append(Y.df,    c(Y1[i,], Y3[i,]))
  }
  df.full  <- data.frame(id = id.df, time = time.df, site = site.df,
                         X0 = X0.df,   A0 = A0.df,     X1 = X1.df,    A1 = A1.df,
                         X2 = X2.df,   A2 = A2.df,      Y = Y.df, 
                         A0.prob = A0.prob.df,  A1.prob = A1.prob.df)
  #### create interaction with time
  df.full$A0.inter.t1 <- df.full$A0*(df.full$time==1)
  df.full$A0.inter.t3 <- df.full$A0*(df.full$time==3)
  df.full$A1.inter.t3 <- df.full$A1*(df.full$time==3)
  df.full$A2.inter.t3 <- df.full$A2*(df.full$time==3)
  
  #### Estimate the treatment history (Correctly specified model)
  #### Since at each time, the treatment is the same for both the left and right ear
  ####   we therefore only need half of the dataset, for the treatment, it is essentially
  ####   got duplicated.
  df.full$site.factor        <- factor(df.full$site)
  df.full$mean.Y1            <- rep(rowMeans(Y1), each = 4) 
  df.full                    <- df.full[order(df.full$site),]
  
  df.full.subset             <- df.full[seq(from = 1, to = nrow(df.full), by = 4),]
  #### A0
  logit.A0                   <- glm(A0 ~ X0 , family = "binomial", data = df.full.subset)
  df.full$est.A0.prob        <- predict.glm(logit.A0, df.full, type="response")
  
  logit.A0.num               <- glm(A0 ~ 1 , family = "binomial", data = df.full.subset)
  df.full$est.A0.prob.num    <- predict.glm(logit.A0.num, df.full, type="response")
  #### A1  
  logit.A1               <- glmer(A1 ~ X0 + X1 + mean.Y1 + A0 + (1 | site), 
                                family = "binomial", data = df.full.subset,
                                control = glmerControl(optimizer ="Nelder_Mead"))
  df.full$est.A1.prob    <- predict(logit.A1, df.full, type="response")
  
  logit.A1.num           <- glm(A1 ~  A0  , 
                                family = "binomial", data = df.full.subset)
  df.full$est.A1.prob.num    <- predict.glm(logit.A1.num, df.full, type="response")
  #### A2
  logit.A2               <- glmer(A2 ~ X0 + X1 + X2 + A0 + A1 + mean.Y1 + (1 | site),
                                family = "binomial", data = df.full.subset,
                                control = glmerControl(optimizer ="Nelder_Mead"))
  df.full$est.A2.prob    <- predict(logit.A2, df.full, type="response")
  
  logit.A2.num           <- glm(A2 ~ A0 + A1 ,
                                family = "binomial", data = df.full.subset)
  df.full$est.A2.prob.num    <- predict.glm(logit.A2.num, df.full, type="response")
  
  weight.time.spec       <- 1/( (df.full$time==1) * ((df.full$A0*df.full$est.A0.prob) + ( (1-df.full$A0)*(1-df.full$est.A0.prob) ) ) +
                                  (df.full$time==3) * ( ((df.full$A0*df.full$est.A0.prob) + ((1-df.full$A0)*(1-df.full$est.A0.prob))) *
                                                          ((df.full$A1*df.full$est.A1.prob) + ((1-df.full$A1)*(1-df.full$est.A1.prob))) *
                                                          ((df.full$A2*df.full$est.A2.prob) + ((1-df.full$A2)*(1-df.full$est.A2.prob))) ))
  est.num.time.spec      <-  (df.full$time==1) * ( (df.full$A0*df.full$est.A0.prob.num) + ((1-df.full$A0)*(1-df.full$est.A0.prob.num)) ) +
    (df.full$time==3) * ( ((df.full$A0*df.full$est.A0.prob.num) + ((1-df.full$A0)*(1-df.full$est.A0.prob.num))) *
                            ((df.full$A1*df.full$est.A1.prob.num) + ((1-df.full$A1)*(1-df.full$est.A1.prob.num))) *
                            ((df.full$A2*df.full$est.A2.prob.num) + ((1-df.full$A2)*(1-df.full$est.A2.prob.num))))
  df.full$weight.spec.stable       <- weight.time.spec*est.num.time.spec
  
  weight.time.fixed       <- 1/( ((df.full$A0*df.full$est.A0.prob) + ((1-df.full$A0)*(1-df.full$est.A0.prob))) *
                                   ((df.full$A1*df.full$est.A1.prob) + ((1-df.full$A1)*(1-df.full$est.A1.prob))) *
                                   ((df.full$A2*df.full$est.A2.prob) + ((1-df.full$A2)*(1-df.full$est.A2.prob))) )
  est.num.time.fixed      <- ((df.full$A0*df.full$est.A0.prob.num) + ((1-df.full$A0)*(1-df.full$est.A0.prob.num))) *
    ((df.full$A1*df.full$est.A1.prob.num) + ((1-df.full$A1)*(1-df.full$est.A1.prob.num))) *
    ((df.full$A2*df.full$est.A2.prob.num) + ((1-df.full$A2)*(1-df.full$est.A2.prob.num)))
  df.full$weight.fixed.stable       <- weight.time.fixed*est.num.time.fixed
  
  


  ###########################################################
  ############ Include the test site correlation #############
  ###########################################################
  #### working independence:
  #### time and individual-specific

  msm.site.fit.ind.spec       <- geeglm(Y ~ time  + A0.inter.t3 + A1.inter.t3 +
                                            A2.inter.t3 +  A0.inter.t1, id = id, 
                                            weights = weight.spec.stable, corstr = "independence",
                                            data = df.full)
  
  #### working exchangeable:
  msm.site.fit.exch      <- geeglm(Y ~ time  + A0.inter.t3 + A1.inter.t3 +
                                       A2.inter.t3 +  A0.inter.t1, id = id, 
                                       weights = weight.fixed.stable, corstr = "exchangeable",
                                       data = df.full)

  
  #### working AR1:
  msm.site.fit.ar1       <- geeglm(Y ~ time  + A0.inter.t3 + A1.inter.t3 +
                                       A2.inter.t3 +  A0.inter.t1, id = id, 
                                       weights = weight.fixed.stable, corstr = "ar1",
                                       data = df.full)
  
  df.full$wave <- rep(1:4, M*N_m)
  #### To model the correlation, we need a contrast matrix stacking
  #### all correlation coefficients from all clusters 
  ####   and a cofficient vector:[alpha1:2, alpha1:3,...,alpha(n_k-1):n_k]
  #### For this contrast matrix, the column represents the 
  #### correlation coefficients alpha1:2, alpha1:3, alpha1:4....
  #### Each row represents ij-th correlation coefficient of the k-th cluster
  #### for instance: row 1 represents alpha1:2 of the first cluster
  ####               row 2 represnets alpha1:3 of the first cluster
  #### This is following the GEE2.0 idea, where the correlation cofficients
  ####  is treated as the response variable.
  zcor <- genZcor(clusz = table(df.full$id), waves = df.full$wave, corstrv = 4)
  #### We want the correlation 1-3 1-4 2-3 2-4 to be the same
  #### Then zcor should have 3 columns
  z2 <- matrix(NA, nrow(zcor), 3)
  z2[,1] <- zcor[,1] 
  z2[,2] <- zcor[,2] + zcor[,3] + zcor[,4] + zcor[,5]
  z2[,3] <- zcor[,6]
  msm.individual.fit.defined   <- geeglm(Y ~ time  + A0.inter.t3 + A1.inter.t3 +
                                           A2.inter.t3 +  A0.inter.t1, id = id, 
                                         weights = weight.fixed.stable, corstr = "userdefined", zcor = z2,
                                         data = df.full)
  
 
  return(list( msm.site.ind.coef.spec     = msm.site.fit.ind.spec$coefficients,
               msm.site.ind.sd.spec       = summary(msm.site.fit.ind.spec)$coefficients[,2],
               msm.site.exch.coef    = msm.site.fit.exch$coefficients,
               msm.site.exch.sd      = summary(msm.site.fit.exch)$coefficients[,2],
               msm.site.ar1.coef     = msm.site.fit.ar1$coefficients,
               msm.site.ar1.sd       = summary(msm.site.fit.ar1)$coefficients[,2],
               msm.site.user.coef    = msm.individual.fit.defined$coefficients,
               msm.site.user.sd      = summary(msm.individual.fit.defined)$coefficients[,2]
  ))
}

### Set the true parameters
set.seed(20201120)
M <- 50
mu0 = 0; sigma0 = 1

beta0.0 = 0.1; beta0.1 = 0.1

gamma1.0 = 0.1; gamma1.1 = 0.1; gamma1.2 = 0.5; gamma1.site = rnorm(M, 0, 0.1)
sigma.X1 = 1

alpha0 = 0.5; alpha1.1 = 0.1
alpha2.1 = 1.8 ; alpha2.2 =  0.995 ; alpha2.3 = -0.55  ; alpha2.4 = 0.8
alpha.site = rnorm(M, 0, 1)
Sigma.Y1 = matrix(c(1,0.7,0.7,1), ncol = 2, byrow = TRUE)

beta1.0 = 0.1; beta1.1 = 0.1; beta1.2 = 0.1; beta1.3 =  0.5
beta1.4 = 0.5; beta1.site = rnorm(M, 0, sqrt(0.1)) 

gamma2.0 = 0.1; gamma2.1 = 0.1; gamma2.2 = 0.5; gamma2.3 = 0.1;
gamma2.4 = 0.5; gamma2.site = rnorm(M,0,1)
sigma.X2 = 1

beta2.0 = 0.1; beta2.1 = 0.1; beta2.2 = 0.1; beta2.3 = 0.1
beta2.4 = 0.2;  beta2.5 = 0.5;  beta2.6 = 0.5; 
beta2.site = rnorm(M, 0, sqrt(0.1)) 

alpha1.2 = 0.1; alpha1.3 = 0.1; alpha1.4 = 0.1
alpha4 = 0.5
Sigma.Y3 = matrix(c(1,0.7,0.7,1), ncol = 2, byrow = TRUE)
param <- MSM.param(mu0,      sigma0,
                   beta0.0,  beta0.1,
                   gamma1.0, gamma1.1,  gamma1.2, gamma1.site,
                   alpha0,   alpha1.1,  alpha2.1, alpha2.2, alpha2.3, alpha2.4,   
                   alpha.site,
                   beta1.0,  beta1.1,   beta1.2,  beta1.3,     beta1.4,  beta1.site,
                   gamma2.0, gamma2.1,  gamma2.2, gamma2.3,    gamma2.4, gamma2.site,
                   beta2.0,  beta2.1,   beta2.2,  beta2.3,     beta2.4,  beta2.5,     beta2.6, beta2.site,
                   alpha1.2, alpha1.3, alpha1.4, alpha4)
param <- param[-c(6:9)]
print(param)
#### Run simulations for 1000 times,
#### using geepack with the option weights
msm.site.ind.coef.spec     <- matrix(NA, nrow = 1000, ncol = 6)
msm.site.ind.sd.spec       <- matrix(NA, nrow = 1000, ncol = 6)

msm.site.exch.coef    <- matrix(NA, nrow = 1000, ncol = 6)
msm.site.exch.sd      <- matrix(NA, nrow = 1000, ncol = 6)

msm.site.ar1.coef     <- matrix(NA, nrow = 1000, ncol = 6)
msm.site.ar1.sd       <- matrix(NA, nrow = 1000, ncol = 6)

msm.site.user.coef     <- matrix(NA, nrow = 1000, ncol = 6)
msm.site.user.sd       <- matrix(NA, nrow = 1000, ncol = 6)
options(nwarnings = 10000)
options(warn=1)
warning.message        <- list()
for(sim in 1:1000){
  print(sim)
  run.msm <- function(M, N_m,
                      mu0,      sigma0,
                      beta0.0,  beta0.1,
                      gamma1.0, gamma1.1,  gamma1.2, gamma1.site, sigma.X1,
                      alpha0,   alpha1.1,  alpha2.1, alpha2.2, alpha2.3, alpha2.4,   
                      alpha.site,  Sigma.Y1,
                      beta1.0,  beta1.1,   beta1.2,  beta1.3,     beta1.4,  beta1.site,
                      gamma2.0, gamma2.1,  gamma2.2, gamma2.3,    gamma2.4, gamma2.site, sigma.X2,
                      beta2.0,  beta2.1,   beta2.2,  beta2.3,     beta2.4,  beta2.5,     beta2.6,  beta2.site,
                      alpha1.2, alpha1.3,  alpha1.4, alpha4,      Sigma.Y3){
    tryCatch(
      expr = {
        temp <- sim.MSM.true.weight(M = M, N_m = N_m,
                                    mu0,      sigma0,
                                    beta0.0,  beta0.1,
                                    gamma1.0, gamma1.1,  gamma1.2, gamma1.site, sigma.X1,
                                    alpha0,   alpha1.1,  alpha2.1, alpha2.2, alpha2.3, alpha2.4,   
                                    alpha.site,  Sigma.Y1,
                                    beta1.0,  beta1.1,   beta1.2,  beta1.3,     beta1.4,  beta1.site,
                                    gamma2.0, gamma2.1,  gamma2.2, gamma2.3,    gamma2.4, gamma2.site, sigma.X2,
                                    beta2.0,  beta2.1,   beta2.2,  beta2.3,     beta2.4,  beta2.5,     beta2.6,  beta2.site,
                                    alpha1.2, alpha1.3,  alpha1.4, alpha4,      Sigma.Y3)
        temp$warnings <- "No.warn"
        return(temp)
      },
      error = function(e){ 
        # (Optional)
        # Do this if an error is caught...
      },
      warning = function(w){
        temp <- sim.MSM.true.weight(M = M, N_m = N_m,
                                    mu0,      sigma0,
                                    beta0.0,  beta0.1,
                                    gamma1.0, gamma1.1,  gamma1.2, gamma1.site, sigma.X1,
                                    alpha0,   alpha1.1,  alpha2.1, alpha2.2, alpha2.3, alpha2.4,   
                                    alpha.site,  Sigma.Y1,
                                    beta1.0,  beta1.1,   beta1.2,  beta1.3,     beta1.4,  beta1.site,
                                    gamma2.0, gamma2.1,  gamma2.2, gamma2.3,    gamma2.4, gamma2.site, sigma.X2,
                                    beta2.0,  beta2.1,   beta2.2,  beta2.3,     beta2.4,  beta2.5,     beta2.6,  beta2.site,
                                    alpha1.2, alpha1.3,  alpha1.4, alpha4,      Sigma.Y3)
        temp$warnings <- conditionMessage(w)
        return(temp)
      },
      finally = {
        # (Optional)
        # Do this at the end before quitting the tryCatch structure...
      }
    )
    
  }
  temp <- run.msm(M = M, N_m = 30,
          mu0,      sigma0,
          beta0.0,  beta0.1,
          gamma1.0, gamma1.1,  gamma1.2, gamma1.site, sigma.X1,
          alpha0,   alpha1.1,  alpha2.1, alpha2.2, alpha2.3, alpha2.4,   
          alpha.site,  Sigma.Y1,
          beta1.0,  beta1.1,   beta1.2,  beta1.3,     beta1.4,  beta1.site,
          gamma2.0, gamma2.1,  gamma2.2, gamma2.3,    gamma2.4, gamma2.site, sigma.X2,
          beta2.0,  beta2.1,   beta2.2,  beta2.3,     beta2.4,  beta2.5,     beta2.6,  beta2.site,
          alpha1.2, alpha1.3,  alpha1.4, alpha4,      Sigma.Y3)
  msm.site.ind.coef.spec[sim,]     <- temp$msm.site.ind.coef.spec
  msm.site.ind.sd.spec[sim,]       <- temp$msm.site.ind.sd.spec
    
  msm.site.exch.coef[sim,]    <- temp$msm.site.exch.coef
  msm.site.exch.sd[sim,]      <- temp$msm.site.exch.sd
    
  msm.site.ar1.coef[sim,]     <- temp$msm.site.ar1.coef
  msm.site.ar1.sd[sim,]       <- temp$msm.site.ar1.sd
    
  msm.site.user.coef[sim,]     <- temp$msm.site.user.coef
  msm.site.user.sd[sim,]       <- temp$msm.site.user.sd
  warning.message[[sim]] <- temp$warnings
}
write.csv(msm.site.ind.coef.spec,"msm_site_ind_coef_spec.csv")
write.csv(msm.site.ind.sd.spec,"msm_site_ind_sd_spec.csv")
write.csv(msm.site.exch.coef,"msm_site_exch_coef.csv")
write.csv(msm.site.exch.sd,"msm_site_exch_sd.csv")
write.csv(msm.site.ar1.coef,"msm_site_ar1_coef.csv")
write.csv(msm.site.ar1.sd,"msm_site_ar1_sd.csv")
write.csv(msm.site.user.coef,"msm_site_user_coef.csv")
write.csv(msm.site.user.sd,"msm_site_user_sd.csv")
write.csv(warning.message, "warning_message.csv")
write.csv(param, "param.csv")
