library(readr)
library(lme4)
library(geepack)
#### The data we used is not publicly available. Therefore, we simulated the data set
oneres <- data.frame(read_csv("Simulated_data.csv"))[,-1]

#ibu: Ibuprofen
#ra: Rheumatoid arthritis 
#asth: Asthma
#hbp: High blood pressure
#db: 	Diabetes mellitus
oneres.sub <- oneres[,c(
  "left4000", "right4000", "left4000_3", "right4000_3",
  "wt09", "wt11", "wt13", "wt15",
  "smk09", "smk11", "smk13", "smk15",
  "bab09", "bab11", "bab13", "bab15",
  "asp11", "asp13", "asp15",
  "bmi11", "bmi13", "bmi15",
  "age11", "age13", "age15",
  "ibu09","ibu11", "ibu13", "ibu15",
  "race", "site", "site_3", "id")]


#### USE 4000
index.outcome.miss <- which(rowSums(is.na(oneres.sub[,1:4]))!=0)
if(length(index.outcome.miss)>0){
  oneres.sub <- oneres.sub[-index.outcome.miss,]
  
  for(j in 1:ncol(oneres.sub)){
    print(paste0(colnames(oneres.sub)[j],":",  sum(is.na(oneres.sub[,j])) ))
  }
  #### impute weight using weight information from previous questionnaire
  index.weight.miss <- which(rowSums(is.na(oneres.sub[,c("wt11", "wt13", "wt15")]))!=0)
  col.index <- which(colnames(oneres.sub) %in% c("wt11", "wt13", "wt15"))
  for(i in index.weight.miss){
    temp       <- oneres.sub[i,]
    which.miss <- col.index[is.na(temp[,col.index])]
    
    for(j in which.miss){
      oneres.sub[i, j] <- oneres.sub[i, (j-1)]
    }
    
  }
  
}

#### impute smoke using smk information from previous questionnaire
smk11 <- oneres.sub$smk11
smk11 <- ifelse(smk11=="current", "1", smk11)
smk11 <- ifelse(smk11=="never", "0", smk11)
smk11 <- ifelse(smk11=="past", "0", smk11)
smk11 <- ifelse(smk11=="missing", NA, smk11)
smk11 <- as.numeric(smk11)
oneres.sub$smk11 <- smk11

oneres.sub$smk09 <- ifelse(oneres.sub$smk09=="No",0, oneres.sub$smk09)
oneres.sub$smk09 <- ifelse(oneres.sub$smk09=="Yes",1, oneres.sub$smk09)
oneres.sub$smk09 <- as.numeric(oneres.sub$smk09)

oneres.sub$smk13 <- ifelse(oneres.sub$smk13=="No",0, oneres.sub$smk13)
oneres.sub$smk13 <- ifelse(oneres.sub$smk13=="Yes",1, oneres.sub$smk13)
oneres.sub$smk13 <- as.numeric(oneres.sub$smk13)

oneres.sub$smk15 <- ifelse(oneres.sub$smk15=="No",0, oneres.sub$smk15)
oneres.sub$smk15 <- ifelse(oneres.sub$smk15=="Yes",1, oneres.sub$smk15)
oneres.sub$smk15 <- as.numeric(oneres.sub$smk15)

index.smk.miss <- which(rowSums(is.na(oneres.sub[,c("smk11", "smk13", "smk15")]))!=0)
if(length(index.smk.miss)>0){
  col.index <- which(colnames(oneres.sub) %in% c("smk11", "smk13", "smk15"))
  for(i in index.smk.miss){
    temp       <- oneres.sub[i,]
    which.miss <- col.index[is.na(temp[,col.index])]
    for(j in which.miss){
      oneres.sub[i, j] <- oneres.sub[i, (j-1)]
    }
  }
}


#### impute asp

oneres.sub$asp11 <- ifelse(oneres.sub$asp11=="No", 0, oneres.sub$asp11)
oneres.sub$asp11 <- ifelse(oneres.sub$asp11=="Yes", 1, oneres.sub$asp11)
oneres.sub$asp11 <- as.numeric(oneres.sub$asp11)
oneres.sub$asp13 <- ifelse(oneres.sub$asp13=="No", 0, oneres.sub$asp13)
oneres.sub$asp13 <- ifelse(oneres.sub$asp13=="Yes", 1, oneres.sub$asp13)
oneres.sub$asp13 <- as.numeric(oneres.sub$asp13)
oneres.sub$asp15 <- ifelse(oneres.sub$asp15=="No", 0, oneres.sub$asp15)
oneres.sub$asp15 <- ifelse(oneres.sub$asp15=="Yes", 1, oneres.sub$asp15)
oneres.sub$asp15 <- as.numeric(oneres.sub$asp15)

oneres.sub$bab09 <- ifelse(oneres.sub$bab09=="No", 0, oneres.sub$bab09)
oneres.sub$bab09 <- ifelse(oneres.sub$bab09=="Yes", 1, oneres.sub$bab09)
oneres.sub$bab09 <- as.numeric(oneres.sub$bab09)
oneres.sub$bab11 <- ifelse(oneres.sub$bab11=="No", 0, oneres.sub$bab11)
oneres.sub$bab11 <- ifelse(oneres.sub$bab11=="Yes", 1, oneres.sub$bab11)
oneres.sub$bab11 <- as.numeric(oneres.sub$bab11)
oneres.sub$bab13 <- ifelse(oneres.sub$bab13=="No", 0, oneres.sub$bab13)
oneres.sub$bab13 <- ifelse(oneres.sub$bab13=="Yes", 1, oneres.sub$bab13)
oneres.sub$bab13 <- as.numeric(oneres.sub$bab13)
oneres.sub$bab15 <- ifelse(oneres.sub$bab15=="No", 0, oneres.sub$bab15)
oneres.sub$bab15 <- ifelse(oneres.sub$bab15=="Yes", 1, oneres.sub$bab15)
oneres.sub$bab15 <- as.numeric(oneres.sub$bab15)

index.asp.miss <- which(rowSums(is.na(oneres.sub[,c("asp11", "asp13", "asp15")]))!=0)
col.index <- which(colnames(oneres.sub) %in% c("asp11", "asp13", "asp15"))
if(length(index.asp.miss)>0){
  for(i in index.asp.miss){
    temp       <- oneres.sub[i,]
    which.miss <- col.index[is.na(temp[,col.index])]
    for(j in which.miss){
      oneres.sub[i, j] <- oneres.sub[i, (j-1)]
    }
  }
  
}

index.bab.miss <- which(rowSums(is.na(oneres.sub[,c("bab11", "bab13", "bab15")]))!=0)
col.index <- which(colnames(oneres.sub) %in% c("bab11", "bab13", "bab15"))

if(length(index.bab.miss)>0){
  for(i in index.bab.miss){
    temp       <- oneres.sub[i,]
    which.miss <- col.index[is.na(temp[,col.index])]
    for(j in which.miss){
      oneres.sub[i, j] <- oneres.sub[i, (j-1)]
    }
  }
  
}


#### impute bmi
index.bmi.miss <- which(rowSums(is.na(oneres.sub[,c("bmi11", "bmi13", "bmi15")]))!=0)
col.index <- which(colnames(oneres.sub) %in% c("bmi11", "bmi13", "bmi15"))
if(length(index.bmi.miss)>0){
  for(i in index.bmi.miss){
    temp       <- oneres.sub[i,]
    which.miss <- col.index[is.na(temp[,col.index])]
    if(min(which.miss) > col.index[1]){
      for(j in which.miss){
        oneres.sub[i, j] <- oneres.sub[i, (j-1)]
      }
    }
  }
}


#### impute buluofen
oneres.sub$ibu09 <- ifelse(oneres.sub$ibu09=="No", 0, oneres.sub$ibu09)
oneres.sub$ibu09 <- ifelse(oneres.sub$ibu09=="Yes", 1, oneres.sub$ibu09)
oneres.sub$ibu09 <- as.numeric(oneres.sub$ibu09)
oneres.sub$ibu11 <- ifelse(oneres.sub$ibu11=="No", 0, oneres.sub$ibu11)
oneres.sub$ibu11 <- ifelse(oneres.sub$ibu11=="Yes", 1, oneres.sub$ibu11)
oneres.sub$ibu11 <- as.numeric(oneres.sub$ibu11)
oneres.sub$ibu13 <- ifelse(oneres.sub$ibu13=="No", 0, oneres.sub$ibu13)
oneres.sub$ibu13 <- ifelse(oneres.sub$ibu13=="Yes", 1, oneres.sub$ibu13)
oneres.sub$ibu13 <- as.numeric(oneres.sub$ibu13)
oneres.sub$ibu15 <- ifelse(oneres.sub$ibu15=="No", 0, oneres.sub$ibu15)
oneres.sub$ibu15 <- ifelse(oneres.sub$ibu15=="Yes", 1, oneres.sub$ibu15)
oneres.sub$ibu15 <- as.numeric(oneres.sub$ibu15)

index.ibu.miss <- which(rowSums(is.na(oneres.sub[,c("ibu11", "ibu13", "ibu15")]))!=0)
col.index <- which(colnames(oneres.sub) %in% c("ibu11", "ibu13", "ibu15"))
if(length(index.ibu.miss)>0){
  for(i in index.ibu.miss){
    temp       <- oneres.sub[i,]
    which.miss <- col.index[is.na(temp[,col.index])]
    for(j in which.miss){
      oneres.sub[i, j] <- oneres.sub[i, (j-1)]
    }
  }
}

index.still.miss <- which(rowSums(is.na(oneres.sub))!=0)
#### further remove 6
if(length(index.still.miss)>0){
  oneres.sub <- oneres.sub[-index.still.miss,]
}

#### cross classify sites
site <- paste0(oneres.sub$site, oneres.sub$site_3)
unique.sites <- unique(site)
oneres.sub$site <- site
site.num <- c()
for(i in unique.sites){
  site.num <- append(site.num, sum(oneres.sub$site==i))
}
#### sites with smaller than 5
small.sites <- unique.sites[which(site.num<5)]
if(length(small.sites) >0){
  oneres.sub2  <- oneres.sub[!(oneres.sub$site%in%small.sites),]
}else{
  oneres.sub2 <- oneres.sub
}
nrow(oneres.sub2)
nrow(oneres)

oneres.sub$asp11 <- oneres.sub$asp11 + oneres.sub$bab11
oneres.sub$asp11 <- ifelse(oneres.sub$asp11 > 0, 1, 0)
oneres.sub$asp13 <- oneres.sub$asp13 + oneres.sub$bab13
oneres.sub$asp13 <- ifelse(oneres.sub$asp13 > 0, 1, 0)
oneres.sub$asp15 <- oneres.sub$asp15 + oneres.sub$bab15
oneres.sub$asp15 <- ifelse(oneres.sub$asp15 > 0, 1, 0)


### Estimate A0
####
logit.A0.denom             <- glm(asp11 ~ wt11 + smk11 + race + age11 + bmi11 + ibu11 , 
                                  family = "binomial", data = oneres.sub2)
oneres.sub2$prob.A0.denom  <- predict(logit.A0.denom, oneres.sub2, type = "response")

logit.A0.numer             <- glm(asp11 ~ 1, 
                                  family = "binomial", data = oneres.sub2)
oneres.sub2$prob.A0.numer  <- predict(logit.A0.numer, oneres.sub2, type = "response")

### Estimate A1
####
oneres.sub2$mean.y <- (oneres.sub2$left4000 + oneres.sub2$right4000)/2

logit.A1.denom             <- glmer(asp13 ~ mean.y + wt13 + smk13 + race + asp11 + age13 + bmi13 + ibu13 + (1|site) ,
                                  family = "binomial", data = oneres.sub2)
oneres.sub2$prob.A1.denom  <- predict(logit.A1.denom, oneres.sub2, type = "response")

logit.A1.numer            <- glm(asp13 ~ asp11 , 
                                 family = "binomial", data = oneres.sub2)
oneres.sub2$prob.A1.numer  <- predict(logit.A1.numer, oneres.sub2, type = "response")

### Estimate A2
####
logit.A2.denom             <- glmer(asp15 ~ mean.y + wt15 + smk15 + race + site + asp11 + asp13 + age15 + bmi15 + ibu15 + (1|site), 
                                  family = "binomial", data = oneres.sub2)
oneres.sub2$prob.A2.denom  <- predict(logit.A2.denom, oneres.sub2, type = "response")

logit.A2.numer            <- glm(asp15 ~ asp13 + asp11 , 
                                 family = "binomial", data = oneres.sub2)
oneres.sub2$prob.A2.numer  <- predict(logit.A2.numer, oneres.sub2, type = "response")

#### Time and individual specific Weights
weights.base  <- (oneres.sub2$asp11)*(oneres.sub2$prob.A0.numer/oneres.sub2$prob.A0.denom) + 
  (1-oneres.sub2$asp11)*((1-oneres.sub2$prob.A0.numer)/(1-oneres.sub2$prob.A0.denom))
weights.foll  <- ((oneres.sub2$asp11)*(oneres.sub2$prob.A0.numer/oneres.sub2$prob.A0.denom) + 
                    (1-oneres.sub2$asp11)*((1-oneres.sub2$prob.A0.numer)/(1-oneres.sub2$prob.A0.denom))) * 
  ((oneres.sub2$asp13)*(oneres.sub2$prob.A1.numer/oneres.sub2$prob.A1.denom) + 
     (1-oneres.sub2$asp13)*((1-oneres.sub2$prob.A1.numer)/(1-oneres.sub2$prob.A1.denom))) * 
  ((oneres.sub2$asp15)*(oneres.sub2$prob.A2.numer/oneres.sub2$prob.A2.denom) + 
     (1-oneres.sub2$asp15)*((1-oneres.sub2$prob.A2.numer)/(1-oneres.sub2$prob.A2.denom)))
oneres.sub2$weights.base <- weights.base
oneres.sub2$weights.foll <- weights.foll

#### Into long format
ids <- oneres.sub2$id
df.full <- data.frame()
for(id in ids){
  temp <- oneres.sub2[oneres.sub2$id==id,]
  df.temp <- data.frame(
    Y4000 = unlist(temp[c("left4000", "right4000", "left4000_3", "right4000_3")]),
    A0 = unlist(rep(temp["asp11"], 4)),
    A1 = unlist(rep(temp["asp13"], 4)),
    A2 = unlist(rep(temp["asp15"], 4)),
    weights.spec = c(rep(unlist(temp["weights.base"]), 2),
                     rep(unlist(temp["weights.foll"]), 2)),
    weights.fixed = rep(unlist(temp["weights.foll"]), 4),
    time = rep(c(1,3), each = 2),
    id   = rep(unlist(temp["id"]), 4),
    site = rep(unlist(temp["site"]), 4))
  df.full <- rbind(df.full, df.temp)
}
rownames(df.full) <- c()

df.full$A0.inter.t1 <- df.full$A0*(df.full$time==1)
df.full$A0.inter.t3 <- df.full$A0*(df.full$time==3)
df.full$A1.inter.t3 <- df.full$A1*(df.full$time==3)
df.full$A2.inter.t3 <- df.full$A2*(df.full$time==3)

df.full <- df.full[order(df.full$site, df.full$id, df.full$time), ]

weight.fixed.site   <- c()
for(i in unique(df.full$site)){
  df.temp               <- df.full[df.full$site==i,]
  if(nrow(df.temp) > 0){
    df.temp2              <- df.temp[seq(from = 1, to = nrow(df.temp), by = 4),]
    weight.temp.stable    <- rep(prod(df.temp2$weights.fixed), nrow(df.temp))
    weight.fixed.site     <- append(weight.fixed.site, weight.temp.stable)
  }
}
df.full$weight.fixed.site <- weight.fixed.site


#########################
#### Frequency @4000 ####
#########################
msm.ind.cov.4000       <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                   A2.inter.t3 +  A0.inter.t1, id = id, 
                                 weights = weights.spec, corstr = "independence",
                                 data = df.full)
msm.ind.cov.4000re <- round(summary(msm.ind.cov.4000)$coef,3)[,1:2]
msm.ind.cov.4000re <- paste0(msm.ind.cov.4000re[,1]," (", msm.ind.cov.4000re[,2], ")")

msm.exch.cov.4000      <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                   A2.inter.t3 +  A0.inter.t1, id = id, 
                                 weights = weights.fixed, corstr = "exchangeable",
                                 data = df.full)
msm.exch.cov.4000re <- round(summary(msm.exch.cov.4000)$coef,3)[,1:2]
msm.exch.cov.4000re <- paste0(msm.exch.cov.4000re[,1]," (", msm.exch.cov.4000re[,2], ")")

point <- sum(summary(msm.exch.cov.4000)$coef[3:5,1])
var.cum <- matrix(1,ncol=3,nrow=1)%*%vcov(msm.exch.cov.4000)[3:5,3:5]%*%matrix(1,ncol=1,nrow=3)
point-1.96*sqrt(var.cum)
point+1.96*sqrt(var.cum)

msm.ar1.cov.4000      <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                  A2.inter.t3 +  A0.inter.t1, id = id, 
                                weights = weights.fixed, corstr = "ar1",
                                data = df.full)
msm.ar1.cov.4000re <- round(summary(msm.ar1.cov.4000)$coef,3)[,1:2]
msm.ar1.cov.4000re <- paste0(msm.ar1.cov.4000re[,1]," (", msm.ar1.cov.4000re[,2], ")")

##### Sigma_1
df.full$wave <- rep(1:4, length(ids))
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
msm.user.cov.4000   <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                A2.inter.t3 +  A0.inter.t1, id = id, 
                              weights = weights.fixed, corstr = "userdefined", zcor = z2,
                              data = df.full)
msm.user.cov.4000re <- round(summary(msm.user.cov.4000)$coef,3)[,1:2]
msm.user.cov.4000re <- paste0(msm.user.cov.4000re[,1]," (", msm.user.cov.4000re[,2], ")")

z2.1 <- z2

##### Sigma_2
df.full$wave <- rep(1:4, length(ids))
zcor <- genZcor(clusz = table(df.full$id), waves = df.full$wave, corstrv = 4)
#### We want the correlation 1-3 1-4 2-3 2-4 to be the same
#### Then zcor should have 3 columns
z2 <- matrix(NA, nrow(zcor), 2)
z2[,1] <- zcor[,1] + zcor[,6]
z2[,2] <- zcor[,2] + zcor[,3] + zcor[,4] + zcor[,5]
msm.user.cov.4000.sigma2   <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                       A2.inter.t3 +  A0.inter.t1, id = id, 
                                     weights = weights.fixed, corstr = "userdefined", zcor = z2,
                                     data = df.full)
msm.user2.cov.4000re <- round(summary(msm.user.cov.4000.sigma2)$coef,3)[,1:2]
msm.user2.cov.4000re <- paste0(msm.user2.cov.4000re[,1]," (", msm.user2.cov.4000re[,2], ")")
z2.2 <- z2

##### Sigma_3
df.full$wave <- rep(1:4, length(ids))
zcor <- genZcor(clusz = table(df.full$id), waves = df.full$wave, corstrv = 4)
#### We want the correlation 1-3 1-4 2-3 2-4 to be the same
#### Then zcor should have 4 columns
z2 <- matrix(NA, nrow(zcor), 4)
z2[,1] <- zcor[,1] 
z2[,2] <- zcor[,2] + zcor[,5]
z2[,3] <- zcor[,3] + zcor[,4]
z2[,4] <- zcor[,6] 

msm.user.cov.4000.sigma3   <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                       A2.inter.t3 +  A0.inter.t1, id = id, 
                                     weights = weights.fixed, corstr = "userdefined", zcor = z2,
                                     data = df.full)
msm.user3.cov.4000re <- round(summary(msm.user.cov.4000.sigma3)$coef,3)[,1:2]
msm.user3.cov.4000re <- paste0(msm.user3.cov.4000re[,1]," (", msm.user3.cov.4000re[,2], ")")
z2.3 <- z2

##### Sigma_4
df.full$wave <- rep(1:4, length(ids))
zcor <- genZcor(clusz = table(df.full$id), waves = df.full$wave, corstrv = 4)
#### We want the correlation 1-3 1-4 2-3 2-4 to be the same
#### Then zcor should have 3 columns
z2 <- matrix(NA, nrow(zcor), 3)
z2[,1] <- zcor[,1] + zcor[,6] 
z2[,2] <- zcor[,2] + zcor[,5]
z2[,3] <- zcor[,3] + zcor[,4]

msm.user.cov.4000.sigma4   <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                       A2.inter.t3 +  A0.inter.t1, id = id, 
                                     weights = weights.fixed, corstr = "userdefined", zcor = z2,
                                     data = df.full)
msm.user4.cov.4000re <- round(summary(msm.user.cov.4000.sigma4)$coef,3)[,1:2]
msm.user4.cov.4000re <- paste0(msm.user4.cov.4000re[,1]," (", msm.user4.cov.4000re[,2], ")")
z2.4 <- z2

#### Unstructured
msm.user.cov.4000.un   <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                   A2.inter.t3 +  A0.inter.t1, id = id, 
                                 weights = weights.fixed, corstr = "unstructured",
                                 data = df.full)
msm.un.cov.4000re <- round(summary(msm.user.cov.4000.un)$coef,3)[,1:2]
msm.un.cov.4000re <- paste0(msm.un.cov.4000re[,1]," (", msm.un.cov.4000re[,2], ")")

###### meta analysis

M <- length(unique(df.full$site))

param.temp.ind   <- matrix(NA, nrow = M, ncol = 6)
var.temp.ind     <- matrix(NA, nrow = M, ncol = 6)
param.temp.exch  <- matrix(NA, nrow = M, ncol = 6)
var.temp.exch    <- matrix(NA, nrow = M, ncol = 6)
param.temp.ar1   <- matrix(NA, nrow = M, ncol = 6)
var.temp.ar1     <- matrix(NA, nrow = M, ncol = 6)
param.temp.user.1  <- matrix(NA, nrow = M, ncol = 6)
var.temp.user.1    <- matrix(NA, nrow = M, ncol = 6)
param.temp.user.2  <- matrix(NA, nrow = M, ncol = 6)
var.temp.user.2    <- matrix(NA, nrow = M, ncol = 6)
param.temp.user.3  <- matrix(NA, nrow = M, ncol = 6)
var.temp.user.3    <- matrix(NA, nrow = M, ncol = 6)
param.temp.user.4  <- matrix(NA, nrow = M, ncol = 6)
var.temp.user.4    <- matrix(NA, nrow = M, ncol = 6)
param.temp.unstruc  <- matrix(NA, nrow = M, ncol = 6)
var.temp.unstruc    <- matrix(NA, nrow = M, ncol = 6)
which.rank.suff  <- c()


unique.site <- unique(df.full$site)
for(i in 1:M){
  df.temp                <- df.full[df.full$site==unique.site[i],]
  ids                    <- unique(df.temp$id)
  tryCatch(
    expr = {
      #### working independence:
      msm.site.fit.ind.spec  <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                         A2.inter.t3 +  A0.inter.t1, id = id, 
                                       weights = weights.spec, corstr = "independence",
                                       data = df.temp)
      #### working exchangeable:
      msm.site.fit.exch      <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                         A2.inter.t3 +  A0.inter.t1, id = id, 
                                       weights = weights.fixed, corstr = "exchangeable",
                                       data = df.temp)
      #### working AR1:
      msm.site.fit.ar1       <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                         A2.inter.t3 +  A0.inter.t1, id = id, 
                                       weights = weights.fixed, corstr = "ar1",
                                       data = df.temp)
      
      ##### Sigma_1
      df.temp$wave <- rep(1:4, length(ids))
      zcor <- genZcor(clusz = table(df.temp$id), waves = df.temp$wave, corstrv = 4)
      z2 <- matrix(NA, nrow(zcor), 3)
      z2[,1] <- zcor[,1] 
      z2[,2] <- zcor[,2] + zcor[,3] + zcor[,4] + zcor[,5]
      z2[,3] <- zcor[,6]
      
      msm.site.fit.user.1      <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                           A2.inter.t3 +  A0.inter.t1, id = id, 
                                         weights = weights.fixed, corstr = "userdefined", zcor = z2,
                                         data = df.temp)
      
      ##### Sigma_2
      df.temp$wave <- rep(1:4, length(ids))
      zcor <- genZcor(clusz = table(df.temp$id), waves = df.temp$wave, corstrv = 4)
      z2 <- matrix(NA, nrow(zcor), 2)
      z2[,1] <- zcor[,1] + zcor[,6]
      z2[,2] <- zcor[,2] + zcor[,3] + zcor[,4] + zcor[,5]
      
      msm.site.fit.user.2      <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                           A2.inter.t3 +  A0.inter.t1, id = id, 
                                         weights = weights.fixed, corstr = "userdefined", zcor = z2,
                                         data = df.temp)
      
      ##### Sigma_3
      df.temp$wave <- rep(1:4, length(ids))
      zcor <- genZcor(clusz = table(df.temp$id), waves = df.temp$wave, corstrv = 4)
      z2 <- matrix(NA, nrow(zcor), 4)
      z2[,1] <- zcor[,1] 
      z2[,2] <- zcor[,2] + zcor[,5]
      z2[,3] <- zcor[,3] + zcor[,4]
      z2[,4] <- zcor[,6] 
      msm.site.fit.user.3      <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                           A2.inter.t3 +  A0.inter.t1, id = id, 
                                         weights = weights.fixed, corstr = "userdefined", zcor = z2,
                                         data = df.temp)
      
      ##### Sigma_4
      df.temp$wave <- rep(1:4, length(ids))
      zcor <- genZcor(clusz = table(df.temp$id), waves = df.temp$wave, corstrv = 4)
      z2 <- matrix(NA, nrow(zcor), 3)
      z2[,1] <- zcor[,1] + zcor[,6] 
      z2[,2] <- zcor[,2] + zcor[,5]
      z2[,3] <- zcor[,3] + zcor[,4]
      msm.site.fit.user.4      <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                           A2.inter.t3 +  A0.inter.t1, id = id, 
                                         weights = weights.fixed, corstr = "userdefined", zcor = z2,
                                         data = df.temp)
      
      msm.unstruc              <- geeglm(Y4000 ~ time  + A0.inter.t3 + A1.inter.t3 +
                                           A2.inter.t3 +  A0.inter.t1, id = id, 
                                         weights = weights.fixed, corstr = "unstructured",
                                         data = df.temp)
      
      #### collect the estimates and variances
      param.temp.ind[i,] <- msm.site.fit.ind.spec$coefficients
      var.temp.ind[i,]   <- (summary(msm.site.fit.ind.spec)$coef[,2])
      
      param.temp.exch[i,] <- msm.site.fit.exch$coefficients
      var.temp.exch[i,]   <- (summary(msm.site.fit.exch)$coef[,2])
      
      param.temp.ar1[i,] <- msm.site.fit.ar1$coefficients
      var.temp.ar1[i,]   <- (summary(msm.site.fit.ar1)$coef[,2])
      
      param.temp.user.1[i,] <- msm.site.fit.user.1$coefficients
      var.temp.user.1[i,]   <- (summary(msm.site.fit.user.1)$coef[,2])
      
      param.temp.user.2[i,] <- msm.site.fit.user.2$coefficients
      var.temp.user.2[i,]   <- (summary(msm.site.fit.user.2)$coef[,2])
      
      param.temp.user.3[i,] <- msm.site.fit.user.3$coefficients
      var.temp.user.3[i,]   <- (summary(msm.site.fit.user.3)$coef[,2])
      
      param.temp.user.4[i,] <- msm.site.fit.user.4$coefficients
      var.temp.user.4[i,]   <- (summary(msm.site.fit.user.4)$coef[,2])
      
      param.temp.unstruc[i,] <- msm.unstruc$coefficients
      var.temp.unstruc[i,]   <- (summary(msm.unstruc)$coef[,2])
      which.rank.suff      <- append(which.rank.suff, i)
    },
    error = function(e){ 
      param.temp.ind[i,] <- NA
      var.temp.ind[i,]   <- NA
      
      param.temp.exch[i,] <- NA
      var.temp.exch[i,]   <- NA
      
      param.temp.ar1[i,] <- NA
      var.temp.ar1[i,]   <- NA
      
      param.temp.user.1[i,] <- NA
      var.temp.user.1[i,]   <- NA
      
      param.temp.user.2[i,] <- NA
      var.temp.user.2[i,]   <- NA
      
      param.temp.user.3[i,] <- NA
      var.temp.user.3[i,]   <- NA
      
      param.temp.user.4[i,] <- NA
      var.temp.user.4[i,]   <- NA
      
      param.temp.unstruc[i,] <- NA
      var.temp.unstruc[i,]   <- NA
    },
    warning = function(w){
      param.temp.ind[i,] <- NA
      var.temp.ind[i,]   <- NA
      
      param.temp.exch[i,] <- NA
      var.temp.exch[i,]   <- NA
      
      param.temp.ar1[i,] <- NA
      var.temp.ar1[i,]   <- NA
      
      param.temp.user.1[i,] <- NA
      var.temp.user.1[i,]   <- NA
      
      param.temp.user.2[i,] <- NA
      var.temp.user.2[i,]   <- NA
      
      param.temp.user.3[i,] <- NA
      var.temp.user.3[i,]   <- NA
      
      param.temp.user.4[i,] <- NA
      var.temp.user.4[i,]   <- NA
      
      param.temp.unstruc[i,] <- NA
      var.temp.unstruc[i,]   <- NA
    },
    finally = {
      # (Optional)
      # Do this at the end before quitting the tryCatch structure...
    }
  )
}
library(meta)
#### Meta analysis on working independence matrix
colnames(param.temp.ind) <- c("delta1","delta2","delta3","delta4","delta5","delta9")
colnames(var.temp.ind)   <- c("delta1.sd","delta2.sd","delta3.sd","delta4.sd","delta5.sd","delta9.sd")
df.ind <- data.frame(cbind(param.temp.ind, var.temp.ind))[which.rank.suff,]
meta.ind.delta1 <- metagen(TE = delta1, seTE = delta1.sd, data = df.ind, method.tau = "REML" )
meta.ind.delta1$TE.random
meta.ind.delta1$seTE.random
meta.ind.delta1$pval.random

meta.ind.delta2 <- metagen(TE = delta2, seTE = delta2.sd, data = df.ind, method.tau = "REML" )
meta.ind.delta2$TE.random
meta.ind.delta2$seTE.random
meta.ind.delta2$pval.random


meta.ind.delta3 <- metagen(TE = delta3, seTE = delta3.sd, data = df.ind, method.tau = "REML" )
meta.ind.delta3$TE.random
meta.ind.delta3$seTE.random
meta.ind.delta3$pval.random


meta.ind.delta4 <- metagen(TE = delta4, seTE = delta4.sd, data = df.ind, method.tau = "REML" )
meta.ind.delta4$TE.random
meta.ind.delta4$seTE.random
meta.ind.delta4$pval.random

meta.ind.delta5 <- metagen(TE = delta5, seTE = delta5.sd, data = df.ind, method.tau = "REML" )
meta.ind.delta5$TE.random
meta.ind.delta5$seTE.random
meta.ind.delta5$pval.random

meta.ind.delta9 <- metagen(TE = delta9, seTE = delta9.sd, data = df.ind, method.tau = "REML" )
meta.ind.delta9$TE.random
meta.ind.delta9$seTE.random
meta.ind.delta9$pval.random

meta.ind.est <- c(round(meta.ind.delta1$TE.random, digits = 3),
                  round(meta.ind.delta2$TE.random, digits = 3),
                  round(meta.ind.delta3$TE.random, digits = 3),
                  round(meta.ind.delta4$TE.random, digits = 3),
                  round(meta.ind.delta5$TE.random, digits = 3),
                  round(meta.ind.delta9$TE.random, digits = 3))

meta.ind.sd  <- c(round(meta.ind.delta1$seTE.random, digits = 3),
                  round(meta.ind.delta2$seTE.random, digits = 3),
                  round(meta.ind.delta3$seTE.random, digits = 3),
                  round(meta.ind.delta4$seTE.random, digits = 3),
                  round(meta.ind.delta5$seTE.random, digits = 3),
                  round(meta.ind.delta9$seTE.random, digits = 3))
mata.ind     <- paste0(meta.ind.est," (", meta.ind.sd, ")")



#### Meta analysis on working exchangeavle matrix
colnames(param.temp.exch) <- c("delta1","delta2","delta3","delta4","delta5","delta9")
colnames(var.temp.exch)   <- c("delta1.sd","delta2.sd","delta3.sd","delta4.sd","delta5.sd","delta9.sd")
df.exch <- data.frame(cbind(param.temp.exch, var.temp.exch))[which.rank.suff,]
meta.exch.delta1 <- metagen(TE = delta1, seTE = delta1.sd, data = df.exch, method.tau = "REML" )
meta.exch.delta1$TE.random
meta.exch.delta1$seTE.random
meta.exch.delta1$pval.random

meta.exch.delta2 <- metagen(TE = delta2, seTE = delta2.sd, data = df.exch, method.tau = "REML" )
meta.exch.delta2$TE.random
meta.exch.delta2$seTE.random
meta.exch.delta2$pval.random

meta.exch.delta3 <- metagen(TE = delta3, seTE = delta3.sd, data = df.exch, method.tau = "REML" )
meta.exch.delta3$TE.random
meta.exch.delta3$seTE.random
meta.exch.delta3$pval.random

meta.exch.delta4 <- metagen(TE = delta4, seTE = delta4.sd, data = df.exch, method.tau = "REML" )
meta.exch.delta4$TE.random
meta.exch.delta4$seTE.random
meta.exch.delta4$pval.random

meta.exch.delta5 <- metagen(TE = delta5, seTE = delta5.sd, data = df.exch, method.tau = "REML" )
meta.exch.delta5$TE.random
meta.exch.delta5$seTE.random
meta.exch.delta5$pval.random

meta.exch.delta9 <- metagen(TE = delta9, seTE = delta9.sd, data = df.exch, method.tau = "REML" )
meta.exch.delta9$TE.random
meta.exch.delta9$seTE.random
meta.exch.delta9$pval.random

meta.exch.est <- c(round(meta.exch.delta1$TE.random, digits = 3),
                   round(meta.exch.delta2$TE.random, digits = 3),
                   round(meta.exch.delta3$TE.random, digits = 3),
                   round(meta.exch.delta4$TE.random, digits = 3),
                   round(meta.exch.delta5$TE.random, digits = 3),
                   round(meta.exch.delta9$TE.random, digits = 3))

meta.exch.sd  <- c(round(meta.exch.delta1$seTE.random, digits = 3),
                   round(meta.exch.delta2$seTE.random, digits = 3),
                   round(meta.exch.delta3$seTE.random, digits = 3),
                   round(meta.exch.delta4$seTE.random, digits = 3),
                   round(meta.exch.delta5$seTE.random, digits = 3),
                   round(meta.exch.delta9$seTE.random, digits = 3))
mata.exch     <- paste0(meta.exch.est," (", meta.exch.sd, ")")


#### Meta analysis on working AR1 matrix
colnames(param.temp.ar1) <- c("delta1","delta2","delta3","delta4","delta5","delta9")
colnames(var.temp.ar1)   <- c("delta1.sd","delta2.sd","delta3.sd","delta4.sd","delta5.sd","delta9.sd")
df.ar1 <- data.frame(cbind(param.temp.ar1, var.temp.ar1))[which.rank.suff,]
meta.ar1.delta1 <- metagen(TE = delta1, seTE = delta1.sd, data = df.ar1, method.tau = "REML" )
meta.ar1.delta1$TE.random
meta.ar1.delta1$seTE.random
meta.ar1.delta1$pval.random

meta.ar1.delta2 <- metagen(TE = delta2, seTE = delta2.sd, data = df.ar1, method.tau = "REML" )
meta.ar1.delta2$TE.random
meta.ar1.delta2$seTE.random
meta.ar1.delta2$pval.random


meta.ar1.delta3 <- metagen(TE = delta3, seTE = delta3.sd, data = df.ar1, method.tau = "REML" )
meta.ar1.delta3$TE.random
meta.ar1.delta3$seTE.random
meta.ar1.delta3$pval.random

meta.ar1.delta4 <- metagen(TE = delta4, seTE = delta4.sd, data = df.ar1, method.tau = "REML" )
meta.ar1.delta4$TE.random
meta.ar1.delta4$seTE.random
meta.ar1.delta4$pval.random

meta.ar1.delta5 <- metagen(TE = delta5, seTE = delta5.sd, data = df.ar1, method.tau = "REML" )
meta.ar1.delta5$TE.random
meta.ar1.delta5$seTE.random
meta.ar1.delta5$pval.random

meta.ar1.delta9 <- metagen(TE = delta9, seTE = delta9.sd, data = df.ar1, method.tau = "REML" )
meta.ar1.delta9$TE.random
meta.ar1.delta9$seTE.random
meta.ar1.delta9$pval.random

meta.ar1.est <- c(round(meta.ar1.delta1$TE.random, digits = 3),
                  round(meta.ar1.delta2$TE.random, digits = 3),
                  round(meta.ar1.delta3$TE.random, digits = 3),
                  round(meta.ar1.delta4$TE.random, digits = 3),
                  round(meta.ar1.delta5$TE.random, digits = 3),
                  round(meta.ar1.delta9$TE.random, digits = 3))

meta.ar1.sd  <- c(round(meta.ar1.delta1$seTE.random, digits = 3),
                  round(meta.ar1.delta2$seTE.random, digits = 3),
                  round(meta.ar1.delta3$seTE.random, digits = 3),
                  round(meta.ar1.delta4$seTE.random, digits = 3),
                  round(meta.ar1.delta5$seTE.random, digits = 3),
                  round(meta.ar1.delta9$seTE.random, digits = 3))
mata.ar1     <- paste0(meta.ar1.est," (", meta.ar1.sd, ")")


#### Meta analysis on working user matrix \Sigma1
colnames(param.temp.user.1) <- c("delta1","delta2","delta3","delta4","delta5","delta9")
colnames(var.temp.user.1)   <- c("delta1.sd","delta2.sd","delta3.sd","delta4.sd","delta5.sd","delta9.sd")
df.user <- data.frame(cbind(param.temp.user.1, var.temp.user.1))[which.rank.suff,]
meta.user.delta1 <- metagen(TE = delta1, seTE = delta1.sd, data = df.user, method.tau = "REML" )
meta.user.delta1$TE.random
meta.user.delta1$seTE.random
meta.user.delta1$pval.random

meta.user.delta2 <- metagen(TE = delta2, seTE = delta2.sd, data = df.user, method.tau = "REML" )
meta.user.delta2$TE.random
meta.user.delta2$seTE.random
meta.user.delta2$pval.random

meta.user.delta3 <- metagen(TE = delta3, seTE = delta3.sd, data = df.user, method.tau = "REML" )
meta.user.delta3$TE.random
meta.user.delta3$seTE.random
meta.user.delta3$pval.random

meta.user.delta4 <- metagen(TE = delta4, seTE = delta4.sd, data = df.user, method.tau = "REML" )
meta.user.delta4$TE.random
meta.user.delta4$seTE.random
meta.user.delta4$pval.random

meta.user.delta5 <- metagen(TE = delta5, seTE = delta5.sd, data = df.user, method.tau = "REML" )
meta.user.delta5$TE.random
meta.user.delta5$seTE.random
meta.user.delta5$pval.random

meta.user.delta9 <- metagen(TE = delta9, seTE = delta9.sd, data = df.user, method.tau = "REML" )
meta.user.delta9$TE.random
meta.user.delta9$seTE.random
meta.user.delta9$pval.random

meta.user.est <- c(round(meta.user.delta1$TE.random, digits = 3),
                   round(meta.user.delta2$TE.random, digits = 3),
                   round(meta.user.delta3$TE.random, digits = 3),
                   round(meta.user.delta4$TE.random, digits = 3),
                   round(meta.user.delta5$TE.random, digits = 3),
                   round(meta.user.delta9$TE.random, digits = 3))

meta.user.sd  <- c(round(meta.user.delta1$seTE.random, digits = 3),
                   round(meta.user.delta2$seTE.random, digits = 3),
                   round(meta.user.delta3$seTE.random, digits = 3),
                   round(meta.user.delta4$seTE.random, digits = 3),
                   round(meta.user.delta5$seTE.random, digits = 3),
                   round(meta.user.delta9$seTE.random, digits = 3))
mata.user     <- paste0(meta.user.est," (", meta.user.sd, ")")

#### Meta analysis on working user matrix \Sigma2
colnames(param.temp.user.2) <- c("delta1","delta2","delta3","delta4","delta5","delta9")
colnames(var.temp.user.2)   <- c("delta1.sd","delta2.sd","delta3.sd","delta4.sd","delta5.sd","delta9.sd")
df.user <- data.frame(cbind(param.temp.user.2, var.temp.user.2))[which.rank.suff,]
meta.user.delta1 <- metagen(TE = delta1, seTE = delta1.sd, data = df.user, method.tau = "REML" )
meta.user.delta1$TE.random
meta.user.delta1$seTE.random
meta.user.delta1$pval.random

meta.user.delta2 <- metagen(TE = delta2, seTE = delta2.sd, data = df.user, method.tau = "REML" )
meta.user.delta2$TE.random
meta.user.delta2$seTE.random
meta.user.delta2$pval.random

meta.user.delta3 <- metagen(TE = delta3, seTE = delta3.sd, data = df.user, method.tau = "REML" )
meta.user.delta3$TE.random
meta.user.delta3$seTE.random
meta.user.delta3$pval.random

meta.user.delta4 <- metagen(TE = delta4, seTE = delta4.sd, data = df.user, method.tau = "REML" )
meta.user.delta4$TE.random
meta.user.delta4$seTE.random
meta.user.delta4$pval.random

meta.user.delta5 <- metagen(TE = delta5, seTE = delta5.sd, data = df.user, method.tau = "REML" )
meta.user.delta5$TE.random
meta.user.delta5$seTE.random
meta.user.delta5$pval.random

meta.user.delta9 <- metagen(TE = delta9, seTE = delta9.sd, data = df.user, method.tau = "REML" )
meta.user.delta9$TE.random
meta.user.delta9$seTE.random
meta.user.delta9$pval.random


meta.user2.est <- c(round(meta.user.delta1$TE.random, digits = 3),
                    round(meta.user.delta2$TE.random, digits = 3),
                    round(meta.user.delta3$TE.random, digits = 3),
                    round(meta.user.delta4$TE.random, digits = 3),
                    round(meta.user.delta5$TE.random, digits = 3),
                    round(meta.user.delta9$TE.random, digits = 3))

meta.user2.sd  <- c(round(meta.user.delta1$seTE.random, digits = 3),
                    round(meta.user.delta2$seTE.random, digits = 3),
                    round(meta.user.delta3$seTE.random, digits = 3),
                    round(meta.user.delta4$seTE.random, digits = 3),
                    round(meta.user.delta5$seTE.random, digits = 3),
                    round(meta.user.delta9$seTE.random, digits = 3))
mata.user2     <- paste0(meta.user2.est," (", meta.user2.sd, ")")
#### Meta analysis on working user matrix \Sigma3
colnames(param.temp.user.3) <- c("delta1","delta2","delta3","delta4","delta5","delta9")
colnames(var.temp.user.3)   <- c("delta1.sd","delta2.sd","delta3.sd","delta4.sd","delta5.sd","delta9.sd")
df.user <- data.frame(cbind(param.temp.user.3, var.temp.user.3))[which.rank.suff,]
meta.user.delta1 <- metagen(TE = delta1, seTE = delta1.sd, data = df.user, method.tau = "REML" )
meta.user.delta1$TE.random
meta.user.delta1$seTE.random
meta.user.delta1$pval.random

meta.user.delta2 <- metagen(TE = delta2, seTE = delta2.sd, data = df.user, method.tau = "REML" )
meta.user.delta2$TE.random
meta.user.delta2$seTE.random
meta.user.delta2$pval.random

meta.user.delta3 <- metagen(TE = delta3, seTE = delta3.sd, data = df.user, method.tau = "REML" )
meta.user.delta3$TE.random
meta.user.delta3$seTE.random
meta.user.delta3$pval.random

meta.user.delta4 <- metagen(TE = delta4, seTE = delta4.sd, data = df.user, method.tau = "REML" )
meta.user.delta4$TE.random
meta.user.delta4$seTE.random
meta.user.delta4$pval.random

meta.user.delta5 <- metagen(TE = delta5, seTE = delta5.sd, data = df.user, method.tau = "REML" )
meta.user.delta5$TE.random
meta.user.delta5$seTE.random
meta.user.delta5$pval.random

meta.user.delta9 <- metagen(TE = delta9, seTE = delta9.sd, data = df.user, method.tau = "REML" )
meta.user.delta9$TE.random
meta.user.delta9$seTE.random
meta.user.delta9$pval.random


meta.user3.est <- c(round(meta.user.delta1$TE.random, digits = 3),
                    round(meta.user.delta2$TE.random, digits = 3),
                    round(meta.user.delta3$TE.random, digits = 3),
                    round(meta.user.delta4$TE.random, digits = 3),
                    round(meta.user.delta5$TE.random, digits = 3),
                    round(meta.user.delta9$TE.random, digits = 3))

meta.user3.sd  <- c(round(meta.user.delta1$seTE.random, digits = 3),
                    round(meta.user.delta2$seTE.random, digits = 3),
                    round(meta.user.delta3$seTE.random, digits = 3),
                    round(meta.user.delta4$seTE.random, digits = 3),
                    round(meta.user.delta5$seTE.random, digits = 3),
                    round(meta.user.delta9$seTE.random, digits = 3))
mata.user3     <- paste0(meta.user3.est," (", meta.user3.sd, ")")
#### Meta analysis on working user matrix \Sigma4
colnames(param.temp.user.4) <- c("delta1","delta2","delta3","delta4","delta5","delta9")
colnames(var.temp.user.4)   <- c("delta1.sd","delta2.sd","delta3.sd","delta4.sd","delta5.sd","delta9.sd")
df.user <- data.frame(cbind(param.temp.user.4, var.temp.user.4))[which.rank.suff,]
meta.user.delta1 <- metagen(TE = delta1, seTE = delta1.sd, data = df.user, method.tau = "REML" )
meta.user.delta1$TE.random
meta.user.delta1$seTE.random
meta.user.delta1$pval.random

meta.user.delta2 <- metagen(TE = delta2, seTE = delta2.sd, data = df.user, method.tau = "REML" )
meta.user.delta2$TE.random
meta.user.delta2$seTE.random
meta.user.delta2$pval.random

meta.user.delta3 <- metagen(TE = delta3, seTE = delta3.sd, data = df.user, method.tau = "REML" )
meta.user.delta3$TE.random
meta.user.delta3$seTE.random
meta.user.delta3$pval.random

meta.user.delta4 <- metagen(TE = delta4, seTE = delta4.sd, data = df.user, method.tau = "REML" )
meta.user.delta4$TE.random
meta.user.delta4$seTE.random
meta.user.delta4$pval.random

meta.user.delta5 <- metagen(TE = delta5, seTE = delta5.sd, data = df.user, method.tau = "REML" )
meta.user.delta5$TE.random
meta.user.delta5$seTE.random
meta.user.delta5$pval.random

meta.user.delta9 <- metagen(TE = delta9, seTE = delta9.sd, data = df.user, method.tau = "REML" )
meta.user.delta9$TE.random
meta.user.delta9$seTE.random
meta.user.delta9$pval.random


meta.user4.est <- c(round(meta.user.delta1$TE.random, digits = 3),
                    round(meta.user.delta2$TE.random, digits = 3),
                    round(meta.user.delta3$TE.random, digits = 3),
                    round(meta.user.delta4$TE.random, digits = 3),
                    round(meta.user.delta5$TE.random, digits = 3),
                    round(meta.user.delta9$TE.random, digits = 3))

meta.user4.sd  <- c(round(meta.user.delta1$seTE.random, digits = 3),
                    round(meta.user.delta2$seTE.random, digits = 3),
                    round(meta.user.delta3$seTE.random, digits = 3),
                    round(meta.user.delta4$seTE.random, digits = 3),
                    round(meta.user.delta5$seTE.random, digits = 3),
                    round(meta.user.delta9$seTE.random, digits = 3))
mata.user4     <- paste0(meta.user4.est," (", meta.user4.sd, ")")

#### unstructured
colnames(param.temp.unstruc) <- c("delta1","delta2","delta3","delta4","delta5","delta9")
colnames(var.temp.unstruc)   <- c("delta1.sd","delta2.sd","delta3.sd","delta4.sd","delta5.sd","delta9.sd")
df.unstruc <- data.frame(cbind(param.temp.unstruc, var.temp.unstruc))[which.rank.suff,]

meta.unstruc.delta1 <- metagen(TE = delta1, seTE = delta1.sd, data = df.unstruc, method.tau = "REML" )
meta.unstruc.delta1$TE.random
meta.unstruc.delta1$seTE.random
meta.unstruc.delta1$pval.random

meta.unstruc.delta2 <- metagen(TE = delta2, seTE = delta2.sd, data = df.unstruc, method.tau = "REML" )
meta.unstruc.delta2$TE.random
meta.unstruc.delta2$seTE.random
meta.unstruc.delta2$pval.random

meta.unstruc.delta3 <- metagen(TE = delta3, seTE = delta3.sd, data = df.unstruc, method.tau = "REML" )
meta.unstruc.delta3$TE.random
meta.unstruc.delta3$seTE.random
meta.unstruc.delta3$pval.random

meta.unstruc.delta4 <- metagen(TE = delta4, seTE = delta4.sd, data = df.unstruc, method.tau = "REML" )
meta.unstruc.delta4$TE.random
meta.unstruc.delta4$seTE.random
meta.unstruc.delta4$pval.random

meta.unstruc.delta5 <- metagen(TE = delta5, seTE = delta5.sd, data = df.unstruc, method.tau = "REML" )
meta.unstruc.delta5$TE.random
meta.unstruc.delta5$seTE.random
meta.unstruc.delta5$pval.random

meta.unstruc.delta9 <- metagen(TE = delta9, seTE = delta9.sd, data = df.unstruc, method.tau = "REML" )
meta.unstruc.delta9$TE.random
meta.unstruc.delta9$seTE.random
meta.unstruc.delta9$pval.random

meta.un.est <- c(round(meta.unstruc.delta1$TE.random, digits = 3),
                 round(meta.unstruc.delta2$TE.random, digits = 3),
                 round(meta.unstruc.delta3$TE.random, digits = 3),
                 round(meta.unstruc.delta4$TE.random, digits = 3),
                 round(meta.unstruc.delta5$TE.random, digits = 3),
                 round(meta.unstruc.delta9$TE.random, digits = 3))

meta.un.sd  <- c(round(meta.unstruc.delta1$seTE.random, digits = 3),
                 round(meta.unstruc.delta2$seTE.random, digits = 3),
                 round(meta.unstruc.delta3$seTE.random, digits = 3),
                 round(meta.unstruc.delta4$seTE.random, digits = 3),
                 round(meta.unstruc.delta5$seTE.random, digits = 3),
                 round(meta.unstruc.delta9$seTE.random, digits = 3))
mata.un     <- paste0(meta.un.est," (", meta.un.sd, ")")

df <- rbind(msm.ind.cov.4000re, 
            msm.exch.cov.4000re,
            msm.ar1.cov.4000re,
            msm.user.cov.4000re,
            msm.user2.cov.4000re,
            msm.user3.cov.4000re,
            msm.user4.cov.4000re,
            msm.un.cov.4000re,
            mata.exch,mata.ar1,
            mata.user, mata.user2,
            mata.user3, mata.user4,
            mata.un)
name <- c("Independent (individual)",
          "Exchangeable (individual)",
          "AR1 (individual)",
          "\\Sigma_1 (individual)",
          "\\Sigma_2 (individual)",
          "\\Sigma_3 (individual)",
          "\\Sigma_4 (individual)",
          "Unstructured (individual)",
          "Exchangeable Meta",
          "AR1 Meta",
          "\\Sigma_1 Meta",
          "\\Sigma_2 Meta",
          "\\Sigma_3 Meta",
          "\\Sigma_4 Meta",
          "Unstructured Meta")
df <- cbind(name, df)
write.csv(df, file="mixed.csv")
