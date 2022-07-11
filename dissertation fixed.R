## Please install these 
#install.packages("doSNOW")
#install.packages("doParallel")

### COuld you please run thsi first to update all packages?
#update.packages()
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(faux)
library(multilevel)
library(dplyr)
library(lattice) 
library(car)
library(lme4)
library(texreg)
library(MASS)
library(multilevel)
library(misty)
library(data.table)
library(lavaan)
library(lme4)
library(MplusAutomation)
library(lavaan)
library(haven)
library(lmerTest)
library(SurvDisc)
library(arm)
library(foreach)
library(doSNOW)
library(doParallel)
library(semTools)
library(lsr)
library(effectsize)
library(heplots)
options(scipen=999)


## Data Simulation begins here: 
ismail_data_sim.v2.2 <- function (varx, vary, icc_x, icc_y,corl1,corl2,nobs,NJ) 
{
  #number of observations per school
  nobs <-nobs
  #Number of schools
  J <- NJ
  ## Defining Simulation Parameters
  ## Addind Noise for Level-1 sample sizes
  MIN = nobs* 80/100
  MAX = nobs*120/100
  newN <- round(runif(J, min = MIN, max = MAX)) ## This is a list of the new samples size for each sample
  datasize <- sum(newN) ## Total number of observations
  
  var_x <- varx
  var_y <-vary
  ICCx <- icc_x
  ICCy <- icc_y
  
  cor_level_1 <- corl1
  cor_level_2 <- corl2
  nij <- datasize
  nUoj <- J
  varx_l2 <- ICCx*(var_x)
  varx_l1 <- (1-ICCx)*(var_x)
  vary_l2 <- ICCx*(var_y)
  vary_l1 <- (1-ICCx)*(var_y)
  # Setting the R Matrix and simulating eijs
  desiredcov_level_1 = (cor_level_1)*sqrt(varx_l1*vary_l1)
  
  R <- matrix(c(varx_l1,desiredcov_level_1,desiredcov_level_1,vary_l1),nrow=2,ncol=2)
  mu <- c(X=0,Y=0)
  eijs <- mvtnorm::rmvnorm(nij,mean=mu,sigma = R)
  eijs <- as.data.frame(eijs)
  ## Setting G Matrix and simulating Uojs( for both X and Y)
  desiredcov_level_2 = (cor_level_2)*sqrt(varx_l2*vary_l2)
  G <- matrix(c(varx_l2,desiredcov_level_2,desiredcov_level_2,vary_l2),nrow=2,ncol=2)
  mu <- c(X=0,Y=0)
  UOjs<- mvtnorm::rmvnorm(J,mean=mu,sigma = G)
  UOjs <- as.data.frame(UOjs)
  newN <-  as.data.frame(newN)
  index <- rownames(UOjs )
  UOjs <- cbind(index=index, UOjs )
  UOjs  <- cbind(UOjs,newN)
  UOjs  <- data.frame(lapply(UOjs , rep, UOjs$newN))
  
  
  UOjs  <-UOjs  %>% 
    arrange(X,Y)
  
  fixed.intercepts <- c(0)
  fixed.intercepts <- as.data.frame(fixed.intercepts)
  row = 1
  times = J
  fixed.intercepts <- fixed.intercepts[rep(row, times),]
  fixed.intercepts <- as.data.frame(fixed.intercepts)
  index <- rownames(fixed.intercepts)
  fixed.intercepts <- cbind(index=index, fixed.intercepts)
  fixed.intercepts <- cbind(fixed.intercepts,newN)
  fixed.intercepts <- data.frame(lapply(fixed.intercepts, rep, fixed.intercepts$newN))
  index <- rownames(fixed.intercepts)
  fixed.intercepts.y <- c(0)
  fixed.intercepts.y  <- as.data.frame(fixed.intercepts.y )
  row = 1
  times = J
  fixed.intercepts.y  <- fixed.intercepts.y [rep(row, times),]
  fixed.intercepts.y  <- as.data.frame(fixed.intercepts.y )
  index <- rownames(fixed.intercepts.y )
  fixed.intercepts.y  <- cbind(index=index, fixed.intercepts.y )
  fixed.intercepts.y  <- cbind(fixed.intercepts.y ,newN)
  fixed.intercepts.y  <- data.frame(lapply(fixed.intercepts.y , rep, fixed.intercepts.y $newN))
  fixed.intercepts = subset(fixed.intercepts, select = -c(newN) )
  fixed.intercepts.y = subset(fixed.intercepts.y, select = -c(newN) )
  ##Creating a new data set with fixed intercepts
  sim.data <- cbind(fixed.intercepts,fixed.intercepts.y)
  sim.data = subset(sim.data, select = -c(index) )
  sim.data.2 <- as.data.frame(sim.data)
  ##Creating BOjs(for X and Y)
  sim.data.2$Bojx <- sim.data.2$fixed.intercepts+UOjs$X
  sim.data.2$Bojy <- sim.data.2$fixed.intercepts.y+UOjs$Y
  sim.data.2$Uojx <- UOjs$X
  sim.data.2$Uojy <- UOjs$Y
  ##Creating overall Xijs and Yijs
  sim.data.2$Xij <- sim.data.2$Bojx+eijs$X
  sim.data.2$Yij <- sim.data.2$Bojy+eijs$Y
  ## Adding eijs to the data frame
  sim.data.2$eijx <- eijs$X
  sim.data.2$eijy <-eijs$Y
  
  sim.data.3 <- sim.data.2 %>% 
    group_by(index) %>% 
    mutate(X.cm = mean(Xij),
           X.cwc = Xij-X.cm,
           Y.cm = mean(Yij)) %>%
    ungroup %>%
    return(sim.data.3)
  sim.data.3 <- as.data.frame(sim.data.3)
}


## Model-1 Multivariate Manifest Covariate Approach Function working with data name only 
run_ri <- function(df) {
  # Only requires input of a data frame
  ##Model-1
  candid.1 <- lmer(Yij~ X.cwc+X.cm + (1|index),data=df,REML= TRUE)
  
  
  Beta.mod.1 <- candid.1@beta
  pvalues<- coef(summary(candid.1))[, 5]
  B.bet.mod1 <- Beta.mod.1[3]
  B.bet.mod1  <- as.data.frame(B.bet.mod1)
  rownames(B.bet.mod1) <- NULL
  B.with.mod1 <- Beta.mod.1[2]
  beta.model1.with <- as.data.frame(B.with.mod1)
  rownames(B.with.mod1) <- NULL
  
  fixed.vars <- se.fixef(candid.1)
  fixed.vars
  se.bet.mod1=fixed.vars[3]
  se.bet.mod1 <- as.data.frame(se.bet.mod1)
  rownames(se.bet.mod1) <- NULL
  se.with.mod1=fixed.vars[2]
  se.with.mod1 <- as.data.frame(se.with.mod1)
  rownames(se.with.mod1) <- NULL
  # get approximate p-values
  pvalues.mod1<- coef(summary(candid.1))[, 5]
  pvalues.mod1 <- as.data.frame(pvalues.mod1)
  p.with.mod1 <- pvalues.mod1[2,]
  p.with.mod1 <- as.data.frame(p.with.mod1)
  rownames(p.with.mod1) <- NULL
  p.bet.mod1 <- pvalues.mod1[3,]
  p.bet.mod1 <- as.data.frame(p.bet.mod1)
  p.bet.mod1 <- if (p.bet.mod1<0.05) {1} else {
    0}
  p.bet.mod1 <- as.data.frame(p.bet.mod1)
  rownames(p.bet.mod1) <- NULL
  re_dat = as.data.frame(VarCorr(candid.1))
  int.var.mod1 = re_dat[1,'vcov']
  int.var.mod1 <- as.data.frame(int.var.mod1)
  rownames(int.var.mod1) <- NULL
  resid.var.mod1 = re_dat[2,'vcov']
  resid.var.mod1<- as.data.frame(resid.var.mod1)
  rownames(resid.var.mod1) <- NULL
  
  ##Model-2
  model <- '
level: 1

Yij ~ Xij
Yij ~~ Yij
Xij ~~ Xij
level: 2

Bojdep =~ Yij
Bojindep=~ Xij
Bojdep ~ Bojindep
'
  candid.2 <- sem(model,
                  data = df,std.lv=FALSE,
                  cluster = "index")
  
  estimates.mod.2 <- parameterEstimates(candid.2)
  B.bet.mod2= estimates.mod.2$est[8]
  se.bet.mod2=estimates.mod.2$se[8]
  p.bet.mod2=estimates.mod.2$pvalue[8]
  B.bet.mod2 <- as.data.frame(B.bet.mod2)
  se.bet.mod2 <- as.data.frame(se.bet.mod2)
  p.bet.mod2 <- as.data.frame(p.bet.mod2)
  p.bet.mod2 <- if (p.bet.mod2<0.05) {1} else {
    0}
  p.bet.mod2 <- as.data.frame(p.bet.mod2)
  rownames(se.bet.mod2) <- NULL
  B.with.mod2= estimates.mod.2$est[1]
  se.with.mod2=estimates.mod.2$se[1]
  p.with.mod2=estimates.mod.2$pvalue[1]
  B.with.mod2 <- as.data.frame(B.with.mod2)
  se.with.mod2 <- as.data.frame(se.with.mod2)
  p.with.mod2 <- as.data.frame(p.with.mod2)
  int.var.mod2 = estimates.mod.2$est[12]
  se.intvar.mod2=estimates.mod.2$se[12]
  p.intvar.mod2=estimates.mod.2$pvalue[12]
  int.var.mod2 <- as.data.frame(int.var.mod2)
  se.intvar.mod2<- as.data.frame(se.intvar.mod2)
  p.intvar.mod2<- as.data.frame(p.intvar.mod2)
  residuals(candid.2)
  
  ##Model-3
  ##candid.3.1 <- lme(Yij~X.cwc+X.cm,random=~1|index,data=df)
  icc.modx<-aov(Xij~as.factor(index),data=df)
 
  level2rel <- ICC2(icc.modx)
  ##summary(level2rel)
  level2meanrelmod3 = mean(level2rel) # Getting Level-2 Mean reliability
  Correctionfactor = (1-level2meanrelmod3)*var(df$X.cm)# Creating multiplying factor to fix Observed level-2 predictor error
  # In the model X.cm ~~ Correctionfactor*X.cm this line fixes the error variance with the (1-meanrel)*variance(X.cm)
  
  # In the model X.cm ~~ Correctionfactor*X.cm this line fixes the error variance with the (1-meanrel)*variance(X.cm)
  model.3 <- '
level: 1
Yij ~ X.cwc

level: 2
new.lat =~ X.cm
Bojdep =~ Yij
Bojdep ~ new.lat
new.lat ~~ Correctionfactor*(new.lat)

'
  candid.3 <- sem(model.3,
                  data = df,
                  cluster = "index",fixed.x=TRUE,estimator="ML")
  estimates.mod.3 <- parameterEstimates(candid.3)
  B.bet.mod3= estimates.mod.3$est[8]
  se.bet.mod3=estimates.mod.3$se[8]
  p.bet.mod3=estimates.mod.3$pvalue[8]
  B.bet.mod3 <- as.data.frame(B.bet.mod3)
  se.bet.mod3 <- as.data.frame(se.bet.mod3)
   p.bet.mod3 <- if (p.bet.mod3<0.05) {1} else {
    0}
   p.bet.mod3 <- as.data.frame(p.bet.mod3)
  rownames(se.bet.mod2) <- NULL
  B.with.mod3= estimates.mod.3$est[1]
  se.with.mod3=estimates.mod.3$se[1]
  p.with.mod3=estimates.mod.3$pvalue[1]
  B.with.mod3 <- as.data.frame(B.with.mod3)
  se.with.mod3 <- as.data.frame(se.with.mod3)
  p.with.mod3 <- as.data.frame(p.with.mod3)
  int.var.mod3 = abs(estimates.mod.3$est[13])
  se.intvar.mod3=estimates.mod.3$se[13]
  p.intvar.mod3=estimates.mod.3$pvalue[13]
  int.var.mod3 <- as.data.frame(int.var.mod3)
  se.intvar.mod3<- as.data.frame(se.intvar.mod3)
  p.intvar.mod3<- as.data.frame(p.intvar.mod3)
  newlat.var.mod3 = abs(estimates.mod.3$est[9])
  se.newlatvar.mod3=estimates.mod.3$se[9]
  p.newlatvar.mod3=estimates.mod.3$pvalue[9]
  newlat.var.mod3 <- as.data.frame(newlat.var.mod3)
  se.newlatvar.mod3<- as.data.frame(se.newlatvar.mod3)
  p.newlatvar.mod3<- as.data.frame(p.newlatvar.mod3)
  
  
  ##Model-4
  ## finding Intraclass correlation for Xij
  icc.modx<-aov(Xij~as.factor(index),data=df)
  iccx<- ICC1(icc.modx)
  ## Finding Average sample size for the formula 
  samplesize <- as.numeric(nrow(df))
  numberofgroups <- as.numeric(max(df$index))
  naverg= samplesize/numberofgroups
  ##Getting Beta weight associated with Between Person and ICCx to create their new Beta weight
  B.bet.mod4 <- (B.bet.mod1)/(sqrt(abs((naverg*iccx)/(1+((naverg-1)*iccx)))))
  B.bet.mod4 <- as.data.frame(B.bet.mod4)
  colnames(B.bet.mod4) <- c("B.bet.mod4")
  
  
  
  
  ## model 5:
  
  ## Finding Average sample size for the formula 
  ##samplesize <- as.numeric(nrow(df))
  ##numberofgroups <- as.numeric(max(df$index))
  ##naverg= samplesize/numberofgroups
  ##Getting Beta weight associated with Between Person and ICCx to create their new Beta weight
  ##icc.mod5 <- level2rel$ICC
  ##icc.mod5<- ICC2(icc.modx)
  B.bet.mod5 <- (B.bet.mod1)/sqrt(level2meanrelmod3)
  B.bet.mod5 <- as.data.frame(B.bet.mod5)
  colnames(B.bet.mod5) <- c("B.bet.mod5")
  
  
  
  ## Model 6 begins here
  ## Extracting Standard errors by Group
  model6.summary <-  df %>%
    group_by(index) %>%
    summarise(
      sd = sd(Xij),
      n = n(),
      se = sd / sqrt(n)
    )%>% 
    arrange(index)
  model6.summary <- data.frame(lapply(model6.summary, rep, model6.summary$n))
  model6.summary$se
  model6.summary = subset(model6.summary, select = -c(sd,n) )
  ## Calculating the weights using inverse squared standard errors
  denom <- (model6.summary$se)^2
  model6.summary$weight<- 1/denom
  model6.summary$index <- as.numeric(model6.summary$index) 
  model6.summary <- model6.summary %>% 
    arrange(index)
  df <- cbind(df,model6.summary)
  df = subset(df, select = -c(index) )
  ## Fitting the model using the inverse squared standard errors as weights
  candid.6 <-lmer(Yij~X.cwc+X.cm+(1|index),data=df,weights = weight)
  summary(candid.6)
  Beta.mod.6 <- candid.6@beta
  pvalues<- coef(summary(candid.6))[, 5]
  B.bet.mod6 <- Beta.mod.6[3]
  B.bet.mod6  <- as.data.frame(B.bet.mod6)
  rownames(B.bet.mod1) <- NULL
  B.with.mod6 <- Beta.mod.6[2]
  beta.model1.with <- as.data.frame(B.with.mod1)
  rownames(B.with.mod6) <- NULL
  
  fixed.vars.mod6 <- se.fixef(candid.6)
  fixed.vars.mod6
  se.bet.mod6=fixed.vars.mod6[3]
  se.bet.mod6 <- as.data.frame(se.bet.mod6)
  rownames(se.bet.mod6) <- NULL
  se.with.mod6=fixed.vars.mod6[2]
  se.with.mod6 <- as.data.frame(se.with.mod6)
  rownames(se.with.mod6) <- NULL
  # get approximate p-values
  pvalues.mod6<- coef(summary(candid.6))[, 5]
  pvalues.mod6 <- as.data.frame(pvalues.mod6)
  p.with.mod6 <- pvalues.mod6[2,]
  p.with.mod6 <- as.data.frame(p.with.mod6)
  rownames(p.with.mod6) <- NULL
  p.bet.mod6 <- pvalues.mod6[3,]
  p.bet.mod6<- if (p.bet.mod6 <0.05) {1} else {
    0}
  p.bet.mod6 <- as.data.frame(p.bet.mod6)
  rownames(p.bet.mod1) <- NULL
  re_dat.mod6 = as.data.frame(VarCorr(candid.6))
  int.var.mod6 = re_dat.mod6[1,'vcov']
  int.var.mod6 <- as.data.frame(int.var.mod6)
  rownames(int.var.mod6) <- NULL
  resid.var.mod6 = re_dat.mod6[2,'vcov']
  resid.var.mod6<- as.data.frame(resid.var.mod6)
  rownames(resid.var.mod6) <- NULL
  
  
  ##Model-7:
  
  fm2=nlme::lme.formula(fixed = Yij~X.cwc+X.cm,
                        data = df, random = ~1| index,
                        control=nlme::lmeControl(returnObject=TRUE))
  s1 = simexlme(model=fm2, model.model=df[,c("index","Yij","X.cm","X.cwc")],
                SIMEXvariable="X.cm",respvar="Yij",grpvar="index",corform="~1| index",
                measurement.error=res.sd,measurement.error.resp=res.sd,
                lambda = c(0.5,2),B = 100, fitting.method = "linear",
                jackknife.estimation = FALSE)
  fm2$coefficients$fixed
  Beta.mod.7 <- as.vector(s1$coefficients)
  B.bet.mod7 <- Beta.mod.7[3]
  B.with.mod7 <- Beta.mod.7[2]
  
  data.frame(B.bet.mod1,se.bet.mod1,p.bet.mod1,B.with.mod1,se.with.mod1,
             p.with.mod1,int.var.mod1,resid.var.mod1,B.bet.mod2,se.bet.mod2,p.bet.mod2,B.with.mod2,se.with.mod2,p.with.mod2,int.var.mod2,se.intvar.mod2,p.intvar.mod2,B.bet.mod3,se.bet.mod3,p.bet.mod3,B.with.mod3,se.with.mod3,
             p.with.mod3,int.var.mod3,se.intvar.mod3,p.intvar.mod3,newlat.var.mod3,
             se.newlatvar.mod3,p.newlatvar.mod3,B.bet.mod4,B.bet.mod5,B.bet.mod6,se.bet.mod6,p.bet.mod6,B.with.mod6,se.with.mod6,
             p.with.mod6,int.var.mod6,resid.var.mod6,B.bet.mod7,B.with.mod7,row.names=NULL)
}
## single file to check if the function is running 
##realdata <- ismail_data_sim.v2.2(varx =1,vary =1,icc_x =0.10,icc_y = 0.1,corl1 = -0.2,corl2 =0.5,nobs=40,NJ=700)
##run_ri(realdata)

### Entering Simulation conditions here 
icc_x.1= c(0.05,0.10,0.20)
icc_y.1= 0.10
varx.1= 1
vary.1=1
corl1.1=c(-0.20,0,0.20)
corl2.1=c(-0.5,0,0.5)
nobs.1= c(4,8,12,16,20)
NJ.1=25
set.seed(1234567)
nrep <-1000# number of replication is supposed to be about 1000 but for testing only 10
## Creating a vector which contains conditions to be able to match general results with replication
cond <- expand.grid(varx=varx.1, vary=vary.1, icc_x=icc_x.1, icc_y=icc_y.1,corl1=corl1.1,corl2=corl2.1,nobs=nobs.1,NJ=NJ.1) # all conditions
pb <- txtProgressBar(0, nrep, style = 3)
progress <- function(nrep) setTxtProgressBar(pb, nrep)
opts <- list(progress=progress)
## calling data generation function
dat_gen <-  function(varx, vary, icc_x, icc_y,corl1,corl2,nobs,NJ){ # function calling on the SMD2 function
  ismail_data_sim.v2.2 (varx, vary, icc_x, icc_y,corl1,corl2,nobs,NJ) }

registerDoSNOW(makeCluster(16))
#registerDoParallel()

## Creating one data set which matches results of functions with conditions
crossdat <- foreach(i = 1:nrow(cond)) %:% # you need to pass all needed packages to this function
  foreach(j = seq_len(nrep), .packages = 
            c("ggplot2","dplyr","tidyr","tidyverse","faux","multilevel",
              "dplyr","lattice","car","lme4","texreg","MASS","multilevel"
              ,"misty","data.table","lavaan","lme4","lavaan","haven","lmerTest","SurvDisc","arm","foreach"),
          .combine = "rbind",
          .errorhandling = 'remove',
          progress) %dopar% {
            mydata <-    dat_gen(varx= cond[i, "varx"], vary = cond[i, "vary"], 
                                 icc_x =  cond[i, "icc_x"],icc_y =  cond[i, "icc_y"],
                                 corl1 =  cond[i, "corl1"],corl2 =  cond[i, "corl2"],
                                 nobs = cond[i, "nobs"],NJ=  cond[i, "NJ"])
            #------------------------------------------------------------------
            # Computin all models and extracting the results to one final data set for each replication
            m0 <- run_ri(mydata) 
            
            ##m2 <- run_rimod2(mydata)# this is the function tested above 
            
            #------------------------------------------------------------------
            # creating final data set for each replication with each condition
            
            dat_all  <- cbind(cond[i, ],m0)
            dat_all
            
          } # end 

results.dat.2 <- foreach(i = 1:nrow(cond)) %:% # you need to pass all needed packages to this function
  foreach(j = seq_len(1), 
          .combine = "rbind",
          .errorhandling = 'remove',
          progress) %dopar% {
            
            rep.results <-   colMeans(crossdat[[i]],na.rm=TRUE,dims = 1L)
            rep.results <-t(rep.results)
            rep.results <- as.data.frame(rep.results)
            rep.results$B.bet.mod1.sd <- sd(crossdat[[i]]$B.bet.mod1,na.rm =TRUE)
            rep.results$B.bet.mod2.sd <- sd(crossdat[[i]]$B.bet.mod2,na.rm =TRUE)
            rep.results$B.bet.mod3.sd <- sd(crossdat[[i]]$B.bet.mod3,na.rm =TRUE)
            rep.results$B.bet.mod4.sd <- sd(crossdat[[i]]$B.bet.mod4,na.rm =TRUE)
            rep.results$B.bet.mod6.sd <- sd(crossdat[[i]]$B.bet.mod6,na.rm =TRUE)
            rep.results$B.bet.mod5.sd <- sd(crossdat[[i]]$B.bet.mod5,na.rm =TRUE)
            rep.results$B.bet.mod7.sd <- sd(crossdat[[i]]$B.bet.mod7,na.rm =TRUE)
            
            rep.results$sum.1 <- sum(((crossdat[[i]]$B.bet.mod1-crossdat[[i]]$corl2)^2),na.rm =TRUE)
            rep.results$sq.sum1 <-sqrt(rep.results$sum.1) 
            rep.results$RMSE.mod1 <- (rep.results$sq.sum1)/sqrt(nrep)
            rep.results$sum.2 <- sum(((crossdat[[i]]$B.bet.mod2-crossdat[[i]]$corl2)^2),na.rm =TRUE)
            rep.results$sq.sum2 <-sqrt(rep.results$sum.2)
            rep.results$RMSE.mod2 <- (rep.results$sq.sum2)/sqrt(nrep) 
            rep.results$sum.3 <- sum(((crossdat[[i]]$B.bet.mod3-crossdat[[i]]$corl2)^2),na.rm =TRUE)
            rep.results$sq.sum3 <-sqrt(rep.results$sum.3) 
            rep.results$RMSE.mod3 <- (rep.results$sq.sum3)/sqrt(nrep) 
            rep.results$sum.4 <- sum(((crossdat[[i]]$B.bet.mod4-crossdat[[i]]$corl2)^2),na.rm =TRUE)
            rep.results$sq.sum4 <-sqrt(rep.results$sum.4) 
            rep.results$RMSE.mod4 <- (rep.results$sq.sum4)/sqrt(nrep) 
            rep.results$sum.5 <- sum(((crossdat[[i]]$B.bet.mod5-crossdat[[i]]$corl2)^2),na.rm =TRUE)
            rep.results$sq.sum5 <-sqrt(rep.results$sum.5) 
            rep.results$RMSE.mod5 <- (rep.results$sq.sum5)/sqrt(nrep) 
            rep.results$sum.6 <- sum(((crossdat[[i]]$B.bet.mod6-crossdat[[i]]$corl2)^2),na.rm =TRUE)
            rep.results$sq.sum6 <-sqrt(rep.results$sum.6) 
            rep.results$RMSE.mod6 <- (rep.results$sq.sum6)/sqrt(nrep)  
            rep.results$sum.7 <- sum(((crossdat[[i]]$B.bet.mod7-crossdat[[i]]$corl2)^2),na.rm =TRUE)
            rep.results$sq.sum7 <-sqrt(rep.results$sum.7) 
            rep.results$RMSE.mod7 <- (rep.results$sq.sum7)/sqrt(nrep)
            rep.results$SEratio.mod1 <- (rep.results$se.bet.mod1)/(rep.results$B.bet.mod1.sd)
            rep.results$SEratio.mod2 <- (rep.results$se.bet.mod2)/(rep.results$B.bet.mod2.sd)
            rep.results$SEratio.mod3 <- (rep.results$se.bet.mod3)/(rep.results$B.bet.mod3.sd)
            rep.results$SEratio.mod6 <- (rep.results$se.bet.mod6)/(rep.results$B.bet.mod6.sd)
            rep.results$bias.mod1 <-  rep.results$B.bet.mod1 -  rep.results$corl2 
            rep.results$bias.mod2 <-  rep.results$B.bet.mod2 -  rep.results$corl2
            rep.results$bias.mod3 <-  rep.results$B.bet.mod3 -  rep.results$corl2 
            rep.results$bias.mod4 <-  rep.results$B.bet.mod4 -  rep.results$corl2 
            rep.results$bias.mod5 <-  rep.results$B.bet.mod5 -  rep.results$corl2 
            rep.results$bias.mod6 <-  rep.results$B.bet.mod6 -  rep.results$corl2 
            rep.results$bias.mod7 <-  rep.results$B.bet.mod7 -  rep.results$corl2 
            
            
            
            rep.results <- subset(rep.results, select = -c(varx,vary,NJ) )
            dat_all  <- cbind(rep.results)
            
          }

results.dat.2<- reduce(results.dat.2, bind_rows)
write.csv(results.dat.2,file="C:/Users/ismae/Desktop/Research and Work/Dissertation/averageresults.csv")


summary(results.dat.2$RMSE.mod1)
sd(results.dat.2$RMSE.mod1)
summary(results.dat.2$RMSE.mod2)
sd(results.dat.2$RMSE.mod2)
summary(results.dat.2$RMSE.mod3)
sd(results.dat.2$RMSE.mod3)
summary(results.dat.2$RMSE.mod4)
sd(results.dat.2$RMSE.mod4)
summary(results.dat.2$RMSE.mod5)
sd(results.dat.2$RMSE.mod5)
summary(results.dat.2$RMSE.mod6)
sd(results.dat.2$RMSE.mod6)
summary(results.dat.2$RMSE.mod7)
sd(results.dat.2$RMSE.mod7)
summary(results.dat.2$SEratio.mod1)
sd(results.dat.2$SEratio.mod1)
summary(results.dat.2$SEratio.mod2)
sd(results.dat.2$SEratio.mod2)
summary(results.dat.2$SEratio.mod3)
sd(results.dat.2$SEratio.mod3)
summary(results.dat.2$SEratio.mod6)
sd(results.dat.2$SEratio.mod6)
summary(results.dat.2$bias.mod1)
sd(results.dat.2$bias.mod1)
summary(results.dat.2$bias.mod2)
sd(results.dat.2$bias.mod2)
summary(results.dat.2$bias.mod3)
sd(results.dat.2$bias.mod3)
summary(results.dat.2$bias.mod4)
sd(results.dat.2$bias.mod4)
summary(results.dat.2$bias.mod5)
sd(results.dat.2$bias.mod5)
summary(results.dat.2$bias.mod6)
sd(results.dat.2$bias.mod6)
summary(results.dat.2$bias.mod7)
sd(results.dat.2$bias.mod7)

bias.dat.mod1 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","bias.mod1"))
bias.dat.mod1$model.no <- "model.1"
bias.dat.mod1 <-bias.dat.mod1 %>% 
  rename(
    bias.mod = bias.mod1  )
bias.dat.mod2 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","bias.mod2"))
bias.dat.mod2$model.no <- "model.2"
bias.dat.mod2 <-bias.dat.mod2 %>% 
  rename(
    bias.mod = bias.mod2 )
bias.dat.mod3 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","bias.mod3"))
bias.dat.mod3$model.no <- "model.3"
bias.dat.mod3 <-bias.dat.mod3 %>% 
  rename(
    bias.mod = bias.mod3  )
bias.dat.mod4 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","bias.mod4"))
bias.dat.mod4$model.no <- "model.4"
bias.dat.mod4 <-bias.dat.mod4 %>% 
  rename(
    bias.mod = bias.mod4  )
bias.dat.mod5 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","bias.mod5"))
bias.dat.mod5$model.no <- "model.5"
bias.dat.mod5 <-bias.dat.mod5 %>% 
  rename(
    bias.mod = bias.mod5  )
bias.dat.mod6 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","bias.mod6"))
bias.dat.mod6$model.no <- "model.6"
bias.dat.mod6 <-bias.dat.mod6 %>% 
  rename(
    bias.mod = bias.mod6  )
bias.dat.mod7 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","bias.mod7"))
bias.dat.mod7$model.no <- "model.7"
bias.dat.mod7 <-bias.dat.mod7 %>% 
  rename(
    bias.mod = bias.mod7  )
bias.dene.2 <- rbind(bias.dat.mod1,bias.dat.mod2,bias.dat.mod3,bias.dat.mod4,bias.dat.mod5,bias.dat.mod6,bias.dat.mod7)
write.csv(bias.dene.2,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/biasdene.csv")
bias.dene.2$model.no <- as.factor(bias.dene.2$model.no)
bias.dene.2$icc_x <- as.factor(bias.dene.2$icc_x)
bias.dene.2$corl1 <- as.factor(bias.dene.2$corl1)
bias.dene.2$corl2 <- as.factor(bias.dene.2$corl2)
bias.dene.2$nobs <- as.factor(bias.dene.2$nobs)

bias.effect <- as.data.frame(etasq(ModelBias.an.mod , anova = TRUE,type=3))
write.csv(bias.effect,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/biaseta.csv")
cor_l1 <- bias.dene.2%>%
  group_by(corl1,model.no) %>%
  summarise(
    sd = sd(bias.mod),
    mean= mean(bias.mod)
  )%>% 
  arrange(model.no)
write.csv(cor_l1,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/corl1bias.csv")
cor_l2.bias <- bias.dene.2%>%
  group_by(corl2,model.no) %>%
  summarise(
    sd = sd(bias.mod),
    mean= mean(bias.mod)
  )%>% 
  arrange(model.no)
write.csv(cor_l2.bias,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/corl2bias.csv")
modeleffect.bias <- bias.dene.2%>%
  group_by(model.no) %>%
  summarise(
    sd = sd(bias.mod),
    mean= mean(bias.mod)
  )%>% 
  arrange(model.no)

nobs <- bias.dene.2%>%
  group_by(nobs,model.no) %>%
  summarise(
    sd = sd(bias.mod),
    mean= mean(bias.mod)
  )%>% 
  arrange(model.no)
write.csv(nobs,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/nobsbias.csv")

icc.bias <- bias.dene.2%>%
  group_by(icc_x,model.no) %>%
  summarise(
    sd = sd(bias.mod),
    mean= mean(bias.mod)
  )%>% 
  arrange(model.no)
write.csv(icc.bias,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/iccbias.csv")





write.csv(modeleffect.bias,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/modeleffectbias.csv")
RMSE.dat.mod1 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","RMSE.mod1"))
RMSE.dat.mod1$model.no <- "model.1"
RMSE.dat.mod1 <-RMSE.dat.mod1 %>% 
  rename(
    RMSE.mod = RMSE.mod1  )
summary(RMSE.dat.mod7)
RMSE.dat.mod2 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","RMSE.mod2"))
RMSE.dat.mod2$model.no <- "model.2"
RMSE.dat.mod2 <-RMSE.dat.mod2 %>% 
  rename(
    RMSE.mod = RMSE.mod2 )
RMSE.dat.mod3 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","RMSE.mod3"))
RMSE.dat.mod3$model.no <- "model.3"
RMSE.dat.mod3 <-RMSE.dat.mod3 %>% 
  rename(
    RMSE.mod = RMSE.mod3  )
RMSE.dat.mod4 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","RMSE.mod4"))
RMSE.dat.mod4$model.no <- "model.4"
RMSE.dat.mod4 <-RMSE.dat.mod4 %>% 
  rename(
    RMSE.mod = RMSE.mod4  )
RMSE.dat.mod5 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","RMSE.mod5"))
RMSE.dat.mod5$model.no <- "model.5"
RMSE.dat.mod5 <-RMSE.dat.mod5 %>% 
  rename(
    RMSE.mod = RMSE.mod5  )
RMSE.dat.mod6 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","RMSE.mod6"))
RMSE.dat.mod6$model.no <- "model.6"
RMSE.dat.mod6 <-RMSE.dat.mod6 %>% 
  rename(
    RMSE.mod = RMSE.mod6  )
RMSE.dat.mod7 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","RMSE.mod7"))
RMSE.dat.mod7$model.no <- "model.7"
RMSE.dat.mod7 <-RMSE.dat.mod7 %>% 
  rename(
    RMSE.mod = RMSE.mod7  )
RMSE.dene.2 <- rbind(RMSE.dat.mod1,RMSE.dat.mod2,RMSE.dat.mod3,RMSE.dat.mod4,RMSE.dat.mod5,RMSE.dat.mod6,RMSE.dat.mod7)
RMSE.dene.2$model.no <- as.factor(RMSE.dene.2$model.no)
RMSE.dene.2$icc_x <- as.factor(RMSE.dene.2$icc_x)
RMSE.dene.2$corl1 <- as.factor(RMSE.dene.2$corl1)
RMSE.dene.2$coll2 <- as.factor(RMSE.dene.2$coll2)
RMSE.dene.2$nobs <- as.factor(RMSE.dene.2$nobs)
options(contrasts = c("contr.sum", "contr.poly"))

icceffect.RMSE <- RMSE.dene.2%>%
  group_by(icc_x, model.no) %>%
  summarise(
    sd = sd(RMSE.mod),
    mean= mean(RMSE.mod)
  )%>% 
  arrange(model.no)
write.csv(icceffect.RMSE ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/ICCeffectRMSE.csv")
write.csv(etasq(ModelRMSE.an.mod , anova = TRUE,type=3),file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/ETaRMSE.csv")
nobseffect.RMSE <- RMSE.dene.2%>%
  group_by(nobs, model.no) %>%
  summarise(
    sd = sd(RMSE.mod),
    mean= mean(RMSE.mod)
  )%>% 
  arrange(model.no)
write.csv(nobseffect.RMSE ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/nobseffectRMSE.csv")
modeleffect.RMSE <- RMSE.dene.2%>%
  group_by(model.no) %>%
  summarise(
    sd = sd(RMSE.mod),
    mean= mean(RMSE.mod)
  )%>% 
  arrange(model.no)
write.csv(modeleffect.RMSE ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/modeleffectRMSE.csv")
SEratio.dat.mod1 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","SEratio.mod1"))
SEratio.dat.mod1$model.no <- "model.1"
SEratio.dat.mod1 <-SEratio.dat.mod1 %>% 
  rename(
    SEratio.mod = SEratio.mod1  )
SEratio.dat.mod2 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","SEratio.mod2"))
SEratio.dat.mod2$model.no <- "model.2"
SEratio.dat.mod2 <-SEratio.dat.mod2 %>% 
  rename(
    SEratio.mod = SEratio.mod2 )
SEratio.dat.mod3 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","SEratio.mod3"))
SEratio.dat.mod3$model.no <- "model.3"
SEratio.dat.mod3 <-SEratio.dat.mod3 %>% 
  rename(
    SEratio.mod = SEratio.mod3  )


SEratio.dat.mod6 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","SEratio.mod6"))
SEratio.dat.mod6$model.no <- "model.6"
SEratio.dat.mod6 <-SEratio.dat.mod6 %>% 
  rename(
    SEratio.mod = SEratio.mod6  )

SEratio.dene.2 <- rbind(SEratio.dat.mod1,SEratio.dat.mod2,SEratio.dat.mod3,SEratio.dat.mod6)
SEratio.dene.2$model.no <- as.factor(SEratio.dene.2$model.no)
SEratio.dene.2$icc_x<- as.factor(SEratio.dene.2$icc_x)
SEratio.dene.2$corl1<- as.factor(SEratio.dene.2$corl1)
SEratio.dene.2$corl2 <- as.factor(SEratio.dene.2$corl2)
SEratio.dene.2$nobs<- as.factor(SEratio.dene.2$nobs)

##write.csv(Modelserati.anova,file="C:/Users/ismae/Desktop/Research and Work/Dissertation/Modelserati.anova.csv")

modeleffect.SEratio <- SEratio.dene.2%>%
  group_by( model.no) %>%
  summarise(
    sd = sd(SEratio.mod),
    mean= mean(SEratio.mod)
  )%>% 
  arrange(model.no)
write.csv(modeleffect.SEratio ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/ICCeffectSEratio.csv")
write.csv(etasq(SEratio.an.mod , anova = TRUE,type=3),file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/ETaSEratio.csv")
nobseffect.SEratio <- SEratio.dene.2%>%
  group_by(nobs, model.no) %>%
  summarise(
    sd = sd(SEratio.mod),
    mean= mean(SEratio.mod)
  )%>% 
  arrange(model.no)
write.csv(nobseffect.SEratio ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/nobseffectSEratio.csv")








##write.csv(Modelserationova,file="C:/Users/ismae/Desktop/Research and Work/Dissertation/Modelserationova.csv")
##summary(results.dat.2$B.with.mod6)
##write.csv(results.dat.2,file="C:/Users/ismae/Desktop/Research and Work/Dissertation/results.csv")
##summary(results.dat.2$B.with.mod6)
p.dat.mod1 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","p.bet.mod1"))
p.dat.mod1$model.no <- "model.1"
p.dat.mod1 <-p.dat.mod1 %>% 
  rename(
    p.bet.mod =p.bet.mod1  )
p.dat.mod2 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","p.bet.mod2"))
p.dat.mod2$model.no <- "model.2"
p.dat.mod2 <-p.dat.mod2 %>% 
  rename(
    p.bet.mod =p.bet.mod2  )
p.dat.mod3 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","p.bet.mod3"))
p.dat.mod3$model.no <- "model.3"
p.dat.mod3 <-p.dat.mod3 %>% 
  rename(
    p.bet.mod =p.bet.mod3  )
p.dat.mod6 <- subset(results.dat.2, select=c("icc_x","corl1","corl2","nobs","p.bet.mod6"))
p.dat.mod6$model.no <- "model.6"
p.dat.mod6 <-p.dat.mod6 %>% 
  rename(
    p.bet.mod =p.bet.mod6  )
p.dat.dene.2 <- rbind(p.dat.mod1,p.dat.mod2,p.dat.mod3,p.dat.mod6)
p.dat.dene.2$model.no <- as.factor(p.dat.dene.2$model.no)
p.dat.dene.2$icc_x <- as.factor(p.dat.dene.2$icc_x)
p.dat.dene.2$corl1 <- as.factor(p.dat.dene.2$corl1)
p.dat.dene.2$corl2 <- as.factor(p.dat.dene.2$corl2)
p.dat.dene.2$nobs <- as.factor(p.dat.dene.2$nobs)
p.dat.power.1 <- filter(p.dat.dene.2, corl1 != "0",corl2!="0")
p.dat.power.2 <- filter(p.dat.dene.2, corl1 == "0",corl2!="0")

p.dat.power <- rbind(p.dat.power.1,p.dat.power.2)
p.dat.power$corl2 <- as.factor(p.dat.power$corl2)
p.dat.type1.1 <- filter(p.dat.dene.2, corl1 != "0", corl2 == "0")
p.dat.type1.2 <- filter(p.dat.dene.2, corl1 == "0", corl2 == "0")
p.dat.type1 <- rbind(p.dat.type1.1,p.dat.type1.2)
summary(p.dat.power)
summary(p.dat.type1)
p.dat.power$model.no <- as.factor(p.dat.power$model.no)
p.dat.power$icc_x <- as.factor(p.dat.power$icc_x)
p.dat.power$corl1 <- as.factor(p.dat.power$corl1)
p.dat.power$corl2 <- as.factor(p.dat.power$corl2)
p.dat.power$nobs <- as.factor(p.dat.power$nobs)
p.dat.power <- p.dat.power%>% group_by(model.no) %>% 
  arrange(model.no) 

p.dat.power <-as.data.frame(p.dat.power)


options(contrasts = c("contr.sum", "contr.poly"))
Modelpower.an.mod <- lm(p.bet.mod~ icc_x+nobs+corl1+corl2+model.no+icc_x*model.no+icc_x*corl1+icc_x*corl2+icc_x*nobs+model.no*corl1+model.no*corl2+model.no*nobs+corl1*corl2+corl1*nobs+corl2*nobs+icc_x*model.no*corl1+icc_x*model.no*corl2+icc_x*model.no*nobs+icc_x*corl1*corl2+icc_x*corl2*nobs+model.no*corl1*corl2+model.no*corl1*nobs+model.no*corl2*nobs+corl1*corl2*nobs+icc_x*corl1*nobs+icc_x*model.no*corl1*corl2+icc_x*corl1*corl2*nobs+model.no*corl1*corl2*nobs+icc_x*model.no*corl2*nobs+icc_x*model.no*corl1*nobs,data=p.dat.power)
Anova(Modelpower.an.mod,type="III")
options(contrasts = c("contr.sum","contr.poly"))
ModeltypeI.an.mod <- lm(I(p.bet.mod * 1e6)~ icc_x+nobs+corl1+corl2+model.no+icc_x*model.no+icc_x*corl1+icc_x*corl2+icc_x*nobs+model.no*corl1+model.no*corl2+model.no*nobs+corl1*corl2+corl1*nobs+corl2*nobs+icc_x*model.no*corl1+icc_x*model.no*corl2+icc_x*model.no*nobs+icc_x*corl1*corl2+icc_x*corl2*nobs+model.no*corl1*corl2+model.no*corl1*nobs+model.no*corl2*nobs+corl1*corl2*nobs+icc_x*corl1*nobs+icc_x*model.no*corl1*corl2+icc_x*corl1*corl2*nobs+model.no*corl1*corl2*nobs+icc_x*model.no*corl2*nobs+icc_x*model.no*corl1*nobs,data=p.dat.type1)
Anova(ModeltypeI.an.mod,type="III")


#### SEmi partial Eta Squared Calculations 
###Bias
options(contrasts = c("contr.sum", "contr.poly"))
ModelBias.an.mod <- lm(bias.mod~ icc_x+nobs+corl1+corl2+model.no+icc_x*model.no+icc_x*corl1+icc_x*corl2+icc_x*nobs+model.no*corl1+model.no*corl2+model.no*nobs+corl1*corl2+corl1*nobs+corl2*nobs+icc_x*model.no*corl1+icc_x*model.no*corl2+icc_x*model.no*nobs+icc_x*corl1*corl2+icc_x*corl2*nobs+model.no*corl1*corl2+model.no*corl1*nobs+model.no*corl2*nobs+corl1*corl2*nobs+icc_x*corl1*nobs+icc_x*model.no*corl1*corl2+icc_x*corl1*corl2*nobs+model.no*corl1*corl2*nobs+icc_x*model.no*corl2*nobs+icc_x*model.no*corl1*nobs,data=bias.dene.2)
Anova(ModelBias.an.mod,type=3 )
ModelBias.Anov <- as.data.frame(Anova(ModelBias.an.mod,type=3 ))
ModelBias.Anov <-ModelBias.Anov[-c(1), ]
ModelBias.SStot <- sum(ModelBias.Anov$`Sum Sq`)
ModelBias.Anov$semipartialetasquared <- (ModelBias.Anov$`Sum Sq`)/ModelBias.SStot


##RMSE
options(contrasts = c("contr.sum", "contr.poly"))
ModelRMSE.an.mod <- lm(RMSE.mod~ icc_x+nobs+corl1+corl2+model.no+icc_x*model.no+icc_x*corl1+icc_x*corl2+icc_x*nobs+model.no*corl1+model.no*corl2+model.no*nobs+corl1*corl2+corl1*nobs+corl2*nobs+icc_x*model.no*corl1+icc_x*model.no*corl2+icc_x*model.no*nobs+icc_x*corl1*corl2+icc_x*corl2*nobs+model.no*corl1*corl2+model.no*corl1*nobs+model.no*corl2*nobs+corl1*corl2*nobs+icc_x*corl1*nobs+icc_x*model.no*corl1*corl2+icc_x*corl1*corl2*nobs+model.no*corl1*corl2*nobs+icc_x*model.no*corl2*nobs+icc_x*model.no*corl1*nobs,data=RMSE.dene.2)
Anova(ModelRMSE.an.mod,type=3 )
ModelRMSE.Anov <- as.data.frame(Anova(ModelRMSE.an.mod,type=3 ))
ModelRMSE.Anov <-ModelRMSE.Anov[-c(1), ]
ModelRMSE.SStot <- sum(ModelRMSE.Anov$`Sum Sq`)
ModelRMSE.Anov$semipartialetasquared <- (ModelRMSE.Anov$`Sum Sq`)/ModelRMSE.SStot

## SEratio
options(contrasts = c("contr.sum", "contr.poly"))
ModelSEratio.an.mod <- lm(SEratio.mod~ icc_x+nobs+corl1+corl2+model.no+icc_x*model.no+icc_x*corl1+icc_x*corl2+icc_x*nobs+model.no*corl1+model.no*corl2+model.no*nobs+corl1*corl2+corl1*nobs+corl2*nobs+icc_x*model.no*corl1+icc_x*model.no*corl2+icc_x*model.no*nobs+icc_x*corl1*corl2+icc_x*corl2*nobs+model.no*corl1*corl2+model.no*corl1*nobs+model.no*corl2*nobs+corl1*corl2*nobs+icc_x*corl1*nobs+icc_x*model.no*corl1*corl2+icc_x*corl1*corl2*nobs+model.no*corl1*corl2*nobs+icc_x*model.no*corl2*nobs+icc_x*model.no*corl1*nobs,data=SEratio.dene.2)
Anova(ModelSEratio.an.mod,type=3 )
ModelSEratio.Anov <- as.data.frame(Anova(ModelSEratio.an.mod,type=3 ))
ModelSEratio.Anov <-ModelSEratio.Anov[-c(1), ]
ModelSEratio.SStot <- sum(ModelSEratio.Anov$`Sum Sq`)
ModelSEratio.Anov$semipartialetasquared <- (ModelSEratio.Anov$`Sum Sq`)/ModelSEratio.SStot

## Type-1 Analysis

options(contrasts = c("contr.treatment", "contr.poly"))
p.dat.type1$corl2 <- as.numeric(p.dat.type1$corl2)
type1.an.mod <- lm(p.bet.mod~ icc_x+nobs+corl1+model.no+icc_x*nobs+icc_x*corl1+icc_x*model.no+ nobs*corl1+nobs*model.no+corl1*model.no+icc_x*nobs*corl1+icc_x*nobs*model.no+nobs*corl1*model.no,data=p.dat.type1)
alias( type1.an.mod  )
Anova(type1.an.mod,type=3 )
Modeltype1.Anov <- as.data.frame(Anova(type1.an.mod ,type=3 ))
Modeltype1.Anov <-Modeltype1.Anov[-c(1), ]
Modeltype1.SStot <- sum(Modeltype1.Anov$`Sum Sq`)
Modeltype1.Anov$semipartialetasquared <- (Modeltype1.Anov$`Sum Sq`)/Modeltype1.SStot 
### Power Analysis 
options(contrasts = c("contr.sum", "contr.poly"))
power.an.mod <- lm(p.bet.mod~ icc_x+nobs+corl1+corl2+model.no+icc_x*model.no+icc_x*corl1+icc_x*corl2+icc_x*nobs+model.no*corl1+model.no*corl2+model.no*nobs+corl1*corl2+corl1*nobs+corl2*nobs+icc_x*model.no*corl1+icc_x*model.no*corl2+icc_x*model.no*nobs+icc_x*corl1*corl2+icc_x*corl2*nobs+model.no*corl1*corl2+model.no*corl1*nobs+model.no*corl2*nobs+corl1*corl2*nobs+icc_x*corl1*nobs+icc_x*model.no*corl1*corl2+icc_x*corl1*corl2*nobs+model.no*corl1*corl2*nobs+icc_x*model.no*corl2*nobs+icc_x*model.no*corl1*nobs,data=p.dat.power)
Modelpower.Anov <- as.data.frame(Anova(power.an.mod ,type=3 ))
Modelpower.Anov <-Modelpower.Anov[-c(1), ]
Modelpower.SStot <- sum(Modelpower.Anov$`Sum Sq`)
Modelpower.Anov$semipartialetasquared <- (Modelpower.Anov$`Sum Sq`)/Modelpower.SStot


write.csv(Modelpower.Anov, file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/ETapower.csv")
write.csv(Modeltype1.Anov, file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/ETatype1.csv")
write.csv(ModelSEratio.Anov,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/ETaSERATIO.csv")
write.csv(ModelRMSE.Anov,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/ETaRMSE.csv")
write.csv(ModelBias.Anov,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/BiasRMSE.csv")

### Interaction Exploration
icceffect.corl1.bias<- bias.dene.2 %>%
  group_by(icc_x, corl1) %>%
  summarise(
    sd = sd(bias.mod),
    mean= mean(bias.mod)
  )%>% 
  arrange(icc_x)
write.csv(icceffect.corl1.bias ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/icceffect.corl1.bias.csv")

icceffect.corl2.bias<- bias.dene.2 %>%
  group_by(icc_x, corl2) %>%
  summarise(
    sd = sd(bias.mod),
    mean= mean(bias.mod)
  )%>% 
  arrange(icc_x)
write.csv(icceffect.corl2.bias ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/icceffect.corl2.bias.csv")

icceffect.corl2.model.bias<- bias.dene.2 %>%
  group_by(icc_x, corl2,model.no) %>%
  summarise(
    sd = sd(bias.mod),
    mean= mean(bias.mod)
  )%>% 
  arrange(model.no)
write.csv(icceffect.corl2.model.bias ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/icceffectcorl1model.csv")

nobs.corl1.model.bias<- bias.dene.2 %>%
  group_by(nobs, corl1,model.no) %>%
  summarise(
    sd = sd(bias.mod),
    mean= mean(bias.mod)
  )%>% 
  arrange(model.no)
write.csv(nobs.corl1.model.bias ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/nobscorl1model.csv")

nobs.iccx.model.corl2.bias<- bias.dene.2 %>%
  group_by(nobs, corl2,icc_x,model.no) %>%
  summarise(
    sd = sd(bias.mod),
    mean= mean(bias.mod)
  )%>% 
  arrange(model.no)
write.csv(nobs.iccx.model.corl2.bias ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/nobs.iccx.model.corl2.bias.csv")

nobscorl2.bias<- bias.dene.2 %>%
  group_by(nobs, corl2) %>%
  summarise(
    sd = sd(bias.mod),
    mean= mean(bias.mod)
  )%>% 
  arrange(nobs)
write.csv(nobscorl2.bias ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/nobscorl2.bias.csv")







iccmodel.RMSE<- RMSE.dene.2 %>%
  group_by(icc_x, model.no) %>%
  summarise(
    sd = sd(RMSE.mod),
    mean= mean(RMSE.mod)
  )%>% 
  arrange(model.no)
write.csv(iccmodel.RMSE ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/iccmodel.RMSE.csv")

nobsmodel.RMSE<- RMSE.dene.2 %>%
  group_by(nobs, model.no) %>%
  summarise(
    sd = sd(RMSE.mod),
    mean= mean(RMSE.mod)
  )%>% 
  arrange(model.no)
write.csv(nobsmodel.RMSE ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/nobsmodel.RMSE.csv")

corl1corl2model.RMSE<- RMSE.dene.2 %>%
  group_by(corl1, corl2) %>%
  summarise(
    sd = sd(RMSE.mod),
    mean= mean(RMSE.mod)
  )%>% 
  arrange(corl1)
write.csv(corl1corl2model.RMSE,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/corl1corl2model.RMSE.csv")

iccxnobsmodelmodel.RMSE<- RMSE.dene.2 %>%
  group_by(icc_x,nobs,model.no) %>%
  summarise(
    sd = sd(RMSE.mod),
    mean= mean(RMSE.mod)
  )%>% 
  arrange(model.no,nobs)
write.csv(iccxnobsmodelmodel.RMSE,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/iccxnobsmodelmodel.RMSE.RMSE.csv")

iccxnobscorl1modelmodel.RMSE<- RMSE.dene.2 %>%
  group_by(icc_x,nobs,model.no,corl1) %>%
  summarise(
    sd = sd(RMSE.mod),
    mean= mean(RMSE.mod)
  )%>% 
  arrange(icc_x,model.no)
write.csv(iccxnobscorl1modelmodel.RMSE,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/iccxnobscorl1modelmodel.RMSE.csv")


### Se Ratio interaction analysis 
iccmodel.SEratio<- SEratio.dene.2 %>%
  group_by(icc_x, model.no) %>%
  summarise(
    sd = sd(SEratio.mod),
    mean= mean(SEratio.mod)
  )%>% 
  arrange(model.no)
write.csv(iccmodel.SEratio ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/iccmodel.SEratio.csv")

iccnobs.SEratio<- SEratio.dene.2 %>%
  group_by(icc_x, nobs) %>%
  summarise(
    sd = sd(SEratio.mod),
    mean= mean(SEratio.mod)
  )%>% 
  arrange(nobs)
write.csv(iccnobs.SEratio ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/iccnobs.SEratio.csv")

modelnobs.SEratio<- SEratio.dene.2 %>%
  group_by(model.no ,nobs) %>%
  summarise(
    sd = sd(SEratio.mod),
    mean= mean(SEratio.mod)
  )%>% 
  arrange(nobs)
write.csv(modelnobs.SEratio ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/modelnobs.SEratio.csv")

iccxmodelnobs.SEratio<- SEratio.dene.2 %>%
  group_by(model.no,nobs,icc_x) %>%
  summarise(
    sd = sd(SEratio.mod),
    mean= mean(SEratio.mod)
  )%>% 
  arrange(nobs,model.no)
write.csv(iccxmodelnobs.SEratio ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/iccxmodelnobs.SEratio.csv")

iccxnobscorl1corl2.SEratio<- SEratio.dene.2 %>%
  group_by(icc_x,nobs,corl2,corl1) %>%
  summarise(
    sd = sd(SEratio.mod),
    mean= mean(SEratio.mod)
  )%>% 
  arrange(icc_x)
write.csv(iccxnobscorl1corl2.SEratio,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/iccxnobscorl1corl2.SEratio.csv")

modelnobscorl1corl2.SEratio<- SEratio.dene.2 %>%
  group_by(model.no,nobs,corl2,corl1) %>%
  summarise(
    sd = sd(SEratio.mod),
    mean= mean(SEratio.mod)
  )%>% 
  arrange(nobs)
write.csv(modelnobscorl1corl2.SEratio,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/modelnobscorl1corl2.SEratio.csv")
##Power
corl1corl2.power<- p.dat.power %>%
  group_by(corl1,corl2) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod)
  )%>% 
  arrange(corl1)
write.csv(corl1corl2.power,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/corl1corl2.power.csv")

icc.power<- p.dat.power %>%
  group_by(icc_x) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod)
  )%>% 
  arrange(icc_x)
write.csv(icc.power ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/icc.power.csv")

nobs.power<- p.dat.power %>%
  group_by(nobs) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod)
  )%>% 
  arrange(nobs)
write.csv(nobs.power ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/nobs.power.csv")
##Type-I
model.type1<- p.dat.type1 %>%
  group_by(model.no) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod)
  )%>% 
  arrange(model.no)
write.csv(model.type1 ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/model.type1.csv")







nobsmodel.type1<- p.dat.type1 %>%
  group_by(nobs,model.no) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod)
  )%>% 
  arrange(model.no)
write.csv(nobsmodel.type1 ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/nobsmodel.type1.csv")

iccnobs.type1<- p.dat.type1 %>%
  group_by(icc_x,nobs) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod)
  )%>% 
  arrange(nobs)
write.csv(iccnobs.type1 ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/iccnobs.type1.csv")

corl1model.type1<- p.dat.type1 %>%
  group_by(corl1,model.no) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod)
  )%>% 
  arrange(model.no)
write.csv(corl1model.type1 ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/corl1model.type1.csv")

iccnobsmodel.type1<- p.dat.type1 %>%
  group_by(icc_x,nobs,model.no) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod)
  )%>% 
  arrange(nobs,model.no)
write.csv(corl1corl2model.type1 ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/corl1corl2model.type1.csv")





write.csv(descriptives ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/descriptives.csv")

corl1.p.type1<- p.dat.type1 %>%
  group_by( corl1,model.no) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod)
  )%>% 
  arrange(model.no)
write.csv(corl1.p.type1 ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/corl1effecttype1.csv")

corl2.p.type1<- p.dat.type1 %>%
  group_by( corl2,model.no) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod)
  )%>% 
  arrange(model.no)
write.csv(corl2.p.type1 ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/corl2effecttype1.csv")


p.type1.summary<- p.dat.type1 %>%
  group_by( model.no) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod),
    minimum= min(p.bet.mod),
    maximum= max(p.bet.mod)
  )%>% 
  arrange(model.no)
write.csv(p.type1.summary ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/p.type1.summary.csv")

p.power.summary<- p.dat.power %>%
  group_by( model.no) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod),
    minimum= min(p.bet.mod),
    maximum= max(p.bet.mod)
  )%>% 
  arrange(model.no)
write.csv(p.power.summary ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/p.power.summary.csv")
icc.p.type1.summary<- p.dat.type1 %>%
  group_by( icc_x) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod),
    minimum= min(p.bet.mod),
    maximum= max(p.bet.mod)
  )%>% 
  arrange(icc_x)
write.csv(icc.p.type1.summary ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/icc.p.type1.summary.csv")
nobs.p.type1.summary<- p.dat.type1 %>%
  group_by( nobs) %>%
  summarise(
    sd = sd(p.bet.mod),
    mean= mean(p.bet.mod),
    minimum= min(p.bet.mod),
    maximum= max(p.bet.mod)
  )%>% 
  arrange(nobs)
write.csv(nobs.p.type1.summary ,file = "C:/Users/ismae/Desktop/Research and Work/Dissertation/nobs.p.type1.summary.csv")
