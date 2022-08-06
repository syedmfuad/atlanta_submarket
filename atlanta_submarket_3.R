
library(tidyverse)
library(spdep)
library(spatialreg)
library(rgdal)
library(rgeos)
library(dplyr)
library(sp)
library(rgeos)
library(geosphere)
library(rstatix)
library(gap)

rm(list=ls())

start.time <- Sys.time()

set.seed(12345)
data <- read.csv("atlanta single family new.csv")
data$sqft <- as.numeric(data$sqft)
data <- subset(data, sqft>0)

data$lon <- -1*data$lon
data$plon <- 0
data$plat <- 0

#fully endogenized

data %>% select(salesprice, calcacres, ENGMeanScaleScore15, median_income, age, totbath, rmtot, rmbed, fixtot, fixhalf, fixaddl, 
                fourthquart, pct_renter_occupied, dcdu, pct_above_f, pct_below_f, MATHMeanScaleScore15, pct_HSdiploma, log_median_income,
                sqft, lon, plon, lat, plat, pct_white, pct_black, pct_collegeDegree, median_house_value, HHsize, tract_cover, 
                total_crime_house) -> data

data$fixaddl <- as.factor(data$fixaddl)
data$dcdu <- as.factor(data$dcdu)
data$fourthquart <- as.factor(data$fourthquart)

data <- data[complete.cases(data), ]

#data$pct_change <- ifelse(data$pct_above_f>0, data$pct_above_f, -1*data$pct_below_f)

Y <- I(data$salesprice)/100000

#house attributes

#I(data$pct_above_f), I(data$pct_below_f),

data$pct_above_f <- ifelse(data$pct_above_f > 0.05, 1, 0)
data$pct_below_f <- ifelse(data$pct_below_f > 0.05, 1, 0)

X <- cbind(1, I(log(data$calcacres)), I(data$age), 
           I(data$age*data$age)/1000, I(data$totbath),
           I(data$sqft)/1000, I(data$tract_cover), I(data$total_crime_house),
           I(data$MATHMeanScaleScore15), I(data$pct_black), I(data$pct_above_f), I(data$pct_below_f),
           I(data$lon-data$plon), I(data$lat-data$plat), I((data$lon-data$plon)^2)/1000, I((data$lat-data$plat)^2)/1000)

#demographic variables/mixing variables 

Z <- cbind(1,I(data$pct_black), I(data$pct_renter_occupied),
           I(data$pct_collegeDegree), I(data$MATHMeanScaleScore15), I(data$pct_HSdiploma), I(data$log_median_income))


ols_agg <- lm(Y~X-1);
summary(ols_agg)

#starting values for the hedonic estimates/betas for each type i.e. for mixing algorithm

types <- 3;

beta_start <- matrix(ols_agg$coef,(types*ncol(X)),1);  

#starting values for the gamma estimates for the demographic variables

gamma_start <- matrix(0.01,(1*ncol(Z)),1);

#starting values for sigma
sigma_start <- matrix(sqrt(mean(ols_agg$residuals^2)),types,1)     

#collecting initializing values

val_start <- c(beta_start,gamma_start,sigma_start);

vals <- val_start;
types <- 3;

#convergence criteria comparing new and old estimates:

Iter_conv <- 0.0001;
j <- types; 

#number of independent variables or beta estimates we need to keep track of - so to use when indexing

niv <- ncol(X); 

#number of demographic variables to use when indexing

gvs <- ncol(Z); 

#row dim of aggregate

n <- nrow(X);
conv_cg = 5000; 
conv_cb = 5000; 

par <- val_start

#FnOne prob density  of observing prices given mean of cross product of house attributes and current 
#iteration of hedonic estimates and sigma

FnOne <- function(par,x,y) 
{
  dnorm(y, mean=x%*%par[-1], sd = par[1], log=FALSE)   
}

#FnTwo max prob densities over type probabilities

FnTwo <- function(par,d,x,y)   
{
  pdy <- matrix(0,n,j) 
  b <- par[1:(niv*j)] 
  s <- par[(niv*j+1):((niv+1)*j)]
  for (i in 1:j)
  {	
    pdy[,i] <- FnOne(c(s[i],b[((i-1)*niv+1):(i*niv)]),X,Y)        
  }
  sum(d*log(pdy))
}

#FnThree logit for gamma estimates 

FnThree <- function(g,z)  
{ 
  L <- exp(z%*%g)  	
}

#FnFour max gamma estimates, type probabilities

FnFour <- function(par,d,z,y)   
{
  L <- matrix(0,n,j) 
  L[,1] <- 1
  L[,3] <- 1
  for (m in 1:(j-2))   
  { 
    L[,(m+1)] <- FnThree(par[((m-1)*gvs+1):(m*gvs)],z)     
  }
  Pi <- L / apply(L,1,sum) 
  sum(apply(d*log(Pi),1,sum))  
}

#mixing algorithm

FMM <- function(par,X,Z,y) 
{
  b <- par[1:(j*niv)]; 
  #g <- par[(j*niv+1):((j*(niv+gvs)-gvs))]; 
  g <- par[(j*niv+1):(length(par)-types)];
  #s <- par[-(1:(j*(niv+gvs)-gvs))];
  s <- par[-(1:(length(par)-types))];
  L <- matrix(0,n,j); 
  f <- L; 
  d <- L;
  b <- matrix(b,niv,j); 
  iter <- 0
  
  while (abs(conv_cg) + abs(conv_cb) > Iter_conv)   {   
    
    #store parameter estimates of preceding iteration of mix through loop
    beta_old <- b; 
    gamma_old <- g; 
    
    #counter for while loop
    iter <- iter+1     
    
    for (i in 1:j) 
    { 
      f[,i] <- FnOne(c(s[i],b[,i]),X,Y)		
    }
    for (i in 1:(j-2)) 
    { 
      L[,1] <- 0
      L[,i+1] <- 0
      L[,(i+1)] <- Z%*%g[((i-1)*gvs+1):(i*gvs)] 
    }
    
    #estimate Pi (P) and individual probabilities of belonging to a certain type (d):
    
    P <- exp(L)/(1+apply(exp(L[,(1:j)]),1,sum))    
    
    for (i in 1:n) 
    {	
      d[i,] <- P[i,]*f[i,]/sum(P[i,]*f[i,]) 
    }
    
    #use individual probs (d) to estimate beta (b), gamma (g)
    
    b1 <- matrix(b,(niv*j),1); par1 <- c(b1,s);
    beta_m <- optim(par1,FnTwo,d=d,x=X,y=Y,control=list(fnscale=-1,maxit=100000))
    b <- matrix(beta_m$par[1:(j*niv)],niv,j) 
    
    s <- beta_m$par[(j*niv+1):(j*(niv+1))]
    
    gam_m <- optim(g,FnFour,d,z=Z,Y,control=list(fnscale=-1,maxit=100000))
    g <- gam_m$par
    
    #setting up convergence check
    
    conv_cg <- sum(abs(g-gamma_old)) 
    conv_cb <- sum(abs(b-beta_old))  
    
    #collecting parameter estimates to use to impute LL
    
    par2 <- matrix(b,(niv*j),1)
    par2 <- c(par2,s)
    LL <- FnTwo(par2,d=d,x=X,y=Y) + FnFour(g,d=d,z=Z,y=Y);
    
    #storing 
    
    bvector <- matrix(b,j*niv,1)
    vals_fin <- c(bvector,g,s)	
    dvector <- d
  }
  #collecting parameters for output
  
  out_pars <- list("vals_fin" = vals_fin, "i_type" = d)
  print(b)
  print(g)
  print(iter)
  
  #return list of estimates - index for subsetting in final updating 
  return(out_pars)
}

#calling:

mix <- FMM(val_start,X=X,Z=Z,y=Y)

end.time <- Sys.time()

start.time-end.time

#final updating:

d <- mix$i_type

b <- mix$vals_fin[1:(j*niv)]; 

#g <- mix$vals_fin[(j*niv+1):((j*(niv+gvs)-gvs))];
g <- mix$vals_fin[(j*niv+1):(length(par)-types)];

#s <- mix$vals_fin[-(1:(j*(niv+gvs)-gvs))];
s <- mix$vals_fin[-(1:(length(par)-types))];

b <- matrix(b,niv,j);

b1 <- matrix(b,(niv*j),1); 
par3 <- c(b1,s);

LL <- FnTwo(par3,d=d,x=X,y=Y) + FnFour(g,d=d,z=Z,y=Y);
AIC <- -2*LL+2*niv

########## WRITE CSV ##########

write.csv(b, file="beta3.csv")
write.csv(d, file="prob3.csv")

########## EXTRA ##########

prob <- read.csv("prob3.csv")
datanew <- cbind(data, prob)
datanew$prob <- colnames(datanew[, c(33:35)])[apply(datanew[, c(33:35)],1,which.max)]
datanew$prob <- ifelse(datanew$prob=="V1", 1, ifelse(datanew$prob=="V2", 2, 3))

data1 <- subset(datanew, prob==1)
data2 <- subset(datanew, prob==2)
data3 <- subset(datanew, prob==3)

data1 %>% get_summary_stats(salesprice, MATHMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, tract_cover, total_crime_house, pct_renter_occupied, type = "mean_sd") -> sumdata1
data2 %>% get_summary_stats(salesprice, MATHMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, tract_cover, total_crime_house, pct_renter_occupied, type = "mean_sd") -> sumdata2
data3 %>% get_summary_stats(salesprice, MATHMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, tract_cover, total_crime_house, pct_renter_occupied, type = "mean_sd") -> sumdata3

########## OLS ##########

Y <- I(data1$salesprice)/100000

X <- cbind(1, I(log(data1$calcacres)), I(data1$age), 
           I(data1$age*data1$age)/1000, I(data1$totbath),
           I(data1$sqft)/1000, I(data1$tract_cover), I(data1$total_crime_house),
           I(data1$MATHMeanScaleScore15), I(data1$pct_black), I(data1$pct_above_f), I(data1$pct_below_f),
           I(data1$lon-data1$plon), I(data1$lat-data1$plat), I((data1$lon-data1$plon)^2)/1000, I((data1$lat-data1$plat)^2)/1000)

ols_agg <- lm(Y~X-1);
summary(ols_agg)

deviance(ols_agg)
pred <- predict(ols_agg, newdata = data.frame(X))
sqrt(mean((Y - pred)^2))*100000
extractAIC(ols_agg)


Y <- I(data2$salesprice)/100000

X <- cbind(1, I(log(data2$calcacres)), I(data2$age), 
           I(data2$age*data2$age)/1000, I(data2$totbath),
           I(data2$sqft)/1000, I(data2$tract_cover), I(data2$total_crime_house),
           I(data2$MATHMeanScaleScore15), I(data2$pct_black), I(data2$pct_above_f), I(data2$pct_below_f),
           I(data2$lon-data2$plon), I(data2$lat-data2$plat), I((data2$lon-data2$plon)^2)/1000, I((data2$lat-data2$plat)^2)/1000)

ols_agg <- lm(Y~X-1);
summary(ols_agg)

deviance(ols_agg)
pred <- predict(ols_agg, newdata = data.frame(X))
sqrt(mean((Y - pred)^2))*100000
extractAIC(ols_agg)


Y <- I(data3$salesprice)/100000

X <- cbind(1, I(log(data3$calcacres)), I(data3$age), 
           I(data3$age*data3$age)/1000, I(data3$totbath),
           I(data3$sqft)/1000, I(data3$tract_cover), I(data3$total_crime_house),
           I(data3$MATHMeanScaleScore15), I(data3$pct_black), I(data3$pct_above_f), I(data3$pct_below_f),
           I(data3$lon-data3$plon), I(data3$lat-data3$plat), I((data3$lon-data3$plon)^2)/1000, I((data3$lat-data3$plat)^2)/1000)

ols_agg <- lm(Y~X-1);
summary(ols_agg)

deviance(ols_agg)
pred <- predict(ols_agg, newdata = data.frame(X))
sqrt(mean((Y - pred)^2))*100000
extractAIC(ols_agg)


########## BETAS ##########

beta <- read.csv("beta3.csv")

Y <- I(data$salesprice)/100000

X <- cbind(1, I(log(data$calcacres)), I(data$age), 
           I(data$age*data$age)/1000, I(data$totbath),
           I(data$sqft)/1000, I(data$tract_cover), I(data$total_crime_house),
           I(data$MATHMeanScaleScore15), I(data$pct_black), I(data$pct_above_f), I(data$pct_below_f),
           I(data$lon-data$plon), I(data$lat-data$plat), I((data$lon-data$plon)^2)/1000, I((data$lat-data$plat)^2)/1000)

pred <- (X%*%beta[,2]*prob[,2] + X%*%beta[,3]*prob[,3] + X%*%beta[,4]*prob[,4])
resid <- Y-pred
sqrt(mean((resid)^2))*100000

########## CHOW TEST ##########

Y1 <- I(data1$salesprice)/100000

X1 <- cbind(1, I(log(data1$calcacres)), I(data1$age), 
            I(data1$age*data1$age)/1000, I(data1$totbath),
            I(data1$sqft)/1000, I(data1$tract_cover), I(data1$total_crime_house),
            I(data1$MATHMeanScaleScore15), I(data1$pct_black), I(data1$pct_above_f), I(data1$pct_below_f),
            I(data1$lon-data1$plon), I(data1$lat-data1$plat), I((data1$lon-data1$plon)^2)/1000, I((data1$lat-data1$plat)^2)/1000)

ols_agg1 <- lm(Y1~X1-1);

Y2 <- I(data2$salesprice)/100000

X2 <- cbind(1, I(log(data2$calcacres)), I(data2$age), 
            I(data2$age*data2$age)/1000, I(data2$totbath),
            I(data2$sqft)/1000, I(data2$tract_cover), I(data2$total_crime_house),
            I(data2$MATHMeanScaleScore15), I(data2$pct_black), I(data2$pct_above_f), I(data2$pct_below_f),
            I(data2$lon-data2$plon), I(data2$lat-data2$plat), I((data2$lon-data2$plon)^2)/1000, I((data2$lat-data2$plat)^2)/1000)

ols_agg2 <- lm(Y2~X2-1);

Y3 <- I(data3$salesprice)/100000

X3 <- cbind(1, I(log(data3$calcacres)), I(data3$age), 
            I(data3$age*data3$age)/1000, I(data3$totbath),
            I(data3$sqft)/1000, I(data3$tract_cover), I(data3$total_crime_house),
            I(data3$MATHMeanScaleScore15), I(data3$pct_black), I(data3$pct_above_f), I(data3$pct_below_f),
            I(data3$lon-data3$plon), I(data3$lat-data3$plat), I((data3$lon-data3$plon)^2)/1000, I((data3$lat-data3$plat)^2)/1000)

ols_agg3 <- lm(Y3~X3-1);

chow.test(Y1,X1[,-1],Y2,X2[,-1])




