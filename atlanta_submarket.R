#https://rural-urban.eu/sites/default/files/05_Spatial%20Autocorrelation%20and%20the%20Spatial%20Durbin%20Model_Eilers.pdf
#https://d-nb.info/1188635956/34
#http://rstudio-pubs-static.s3.amazonaws.com/5027_52298866e7924b18b54e5c9a0a21b450.html
#https://stats.stackexchange.com/questions/149415/how-do-i-interpret-lagsarlm-output-from-rs-spdep
#https://rpubs.com/corey_sparks/109650
#https://rpubs.com/corey_sparks/108130
#https://www.stata.com/training/webinar_series/spatial-autoregressive-models/spatial/resource/spillover.html

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

#example with ethiopia

rm(list=ls())

set.seed(12345)
data <- read.csv("atlanta single family.csv")
data$sqft <- as.numeric(data$sqft)
data <- subset(data, sqft>0)

data$lon <- -1*data$lon
data$plon <- 0 #84.3563
data$plat <- 0 #33.7663

data %>% select(salesprice, calcacres, ENGMeanScaleScore15, median_income, age, totbath, stories, centheat, fourthquart, pct_renter_occupied,
                sqft, lon, plon, lat, plat, pct_white, pct_black, pct_collegeDegree, median_house_value, HHsize) -> data
data$lon <- -1*data$lon
data <- data[complete.cases(data), ]

data_sp <- data

coordinates(data_sp) <- ~lon+lat

class(data_sp) #check if SpatialPointsDataFrame
d <- distm(data_sp)

#d <- gDistance(data_sp, byid=T)
min.d <- apply(d, 1, function(x) order(x, decreasing=F)[2])

newdata <- cbind(data, data[min.d,], apply(d, 1, function(x) sort(x, decreasing=F)[2]))

colnames(newdata) <- c(colnames(data), "salesprice_n", "calcacres_n", "ENGMeanScaleScore15_n", "median_income_n", "age_n", "totbath_n",
                       "stories_n", "centheat_n", "fourthquart_n", "pct_renter_occupied_n", "sqft_n", "lon_n", "plon_n", "lat_n", "plat_n",
                       "pct_white_n", "pct_black_n", "pct_collegeDegree_n", "median_house_value_n", "HHsize_n", "distance")


#fully endogenized

data <- newdata

data %>% select(salesprice, calcacres, ENGMeanScaleScore15, median_income, age, totbath, stories, centheat, fourthquart, pct_renter_occupied,
                sqft, lon, plon, lat, plat, pct_white, pct_black, pct_collegeDegree, median_house_value, HHsize) -> data
data <- data[complete.cases(data), ]

Y <- I(data$salesprice)/100000

#house attributes

X <- cbind(1, I(log(data$calcacres)), I(log((data$ENGMeanScaleScore15)/100)), I(log((data$median_income)/10000)), I(data$age), 
           I(data$age*data$age)/1000, I(data$totbath), I(data$stories), I(data$centheat), I(data$fourthquart), I(data$pct_renter_occupied),
           I(data$sqft)/1000, I(data$lon-data$plon), I(data$lat-data$plat), I((data$lon-data$plon)^2)/1000, I((data$lat-data$plat)^2)/1000)

#demographic variables/mixing variables 

Z <- cbind(1,I(data$ENGMeanScaleScore15)/100, I(data$calcacres), I(data$median_income)/10000, I(data$age), I(data$pct_white), I(data$pct_black),
           I(data$pct_collegeDegree), I(data$median_house_value)/100000, I(data$HHsize))


ols_agg <- lm(Y~X-1);
summary(ols_agg)

k = ncol(X)
n=nrow(X)
beta = matrix(NA,nrow=n,ncol=k)
sigma = c(1,rep(NA,n))
psi = rep(NA,n)

set.seed(12345)

#starting values for the hedonic estimates/betas for each type i.e. for mixing algorithm

beta_start <- matrix(ols_agg$coef,(2*ncol(X)),1);  

#starting values for the gamma estimates for the demographic variables

gamma_start <- matrix(0.01,(1*ncol(Z)),1);

#starting values for sigma
sigma_start <- matrix(sqrt(mean(ols_agg$residuals^2)),2,1)     

#collecting initializing values

val_start <- c(beta_start,gamma_start,sigma_start);

vals <- val_start;
types <- 2;

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
  #par[1] <- ifelse(par[1]>0, par[1], -1*par[1])
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
  #pdy[pdy < 0.0005] <- 0.05
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
  for (m in 1:(j-1))   
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
  g <- par[(j*niv+1):((j*(niv+gvs)-gvs))]; 
  s <- par[-(1:(j*(niv+gvs)-gvs))];
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
    for (i in 1:(j-1)) 
    { 
      L[,1] <- 0
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

#final updating:

d <- mix$i_type

b <- mix$vals_fin[1:(j*niv)]; 

g <- mix$vals_fin[(j*niv+1):((j*(niv+gvs)-gvs))]; 

s <- mix$vals_fin[-(1:(j*(niv+gvs)-gvs))];

b <- matrix(b,niv,j);

b1 <- matrix(b,(niv*j),1); 

par3 <- c(b1,s);

########## EXTRA ##########

datanew <- cbind(newdata, d)
datanew$prob <- ifelse(datanew$`1` > datanew$`2`, 1, 2)

data1 <- subset(datanew, prob==1)
data2 <- subset(datanew, prob==2)

data1 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, type = "mean_sd") -> sumdata1
data2 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, type = "mean_sd") -> sumdata2

LL <- FnTwo(par3,d=d,x=X,y=Y) + FnFour(g,d=d,z=Z,y=Y);
AIC <- -2*LL+2*niv

Ds=d;
beta=b;
#bse=cbind(bse1,bse2);
gamma=cbind(g[1:gvs],g[(gvs+1):(2*gvs)]);
#gse=cbind(gse1);

X <- cbind(1, I(log(data$calcacres)), I(log((data$ENGMeanScaleScore15)/100)), I(log((data$median_income)/10000)), I(data$age), 
           I(data$age*data$age)/1000, I(data$totbath), I(data$stories), I(data$centheat), I(data$fourthquart), I(data$pct_renter_occupied),
           I(data$sqft)/1000, I(data$lon-data$plon), I(data$lat-data$plat), I((data$lon-data$plon)^2)/1000, I((data$lat-data$plat)^2)/1000)

pred <- (X%*%b[,1]*Ds[,1] + X%*%b[,2]*Ds[,2])*100000

resid <- Y-pred
rmse2 <- sqrt(mean((resid)^2))

########## EXTRA OLS ##########

prob <- read.csv("prob2.csv")
datanew <- cbind(newdata, prob)
datanew$prob <- ifelse(datanew$V1 > datanew$V2, 1, 2)

data1 <- subset(datanew, prob==1)
data2 <- subset(datanew, prob==2)

data1 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, type = "mean_sd") -> sumdata1
data2 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, type = "mean_sd") -> sumdata2

LL <- FnTwo(par3,d=d,x=X,y=Y) + FnFour(g,d=d,z=Z,y=Y);
AIC <- -2*LL+2*niv

##### OLS #####

data_ols <- subset(datanew, prob==2) #change prob==2

Y <- I(data_ols$salesprice)/100000

X <- cbind(1, I(log(data_ols$calcacres)), I(log((data_ols$ENGMeanScaleScore15)/100)), I(log((data_ols$median_income)/10000)), I(data_ols$age), 
           I(data_ols$age*data_ols$age)/1000, I(data_ols$totbath), I(data_ols$stories), I(data_ols$centheat), I(data_ols$fourthquart),
           I(data_ols$pct_renter_occupied), I(data_ols$sqft)/1000, I(data_ols$lon-data_ols$plon), I(data_ols$lat-data_ols$plat), 
           I((data_ols$lon-data_ols$plon)^2)/1000, I((data_ols$lat-data_ols$plat)^2)/1000)

X <- cbind(1, I(log(data_ols$calcacres)), I(log((data_ols$ENGMeanScaleScore15)/100)), I(log((data_ols$median_income)/10000)), I(data_ols$age), 
           I(data_ols$age*data_ols$age)/1000, I(data_ols$totbath), I(data_ols$stories), I(data_ols$centheat), I(data_ols$fourthquart),
           I(data_ols$pct_renter_occupied), I(data_ols$sqft)/1000, I(data_ols$lon-data_ols$plon), I(data_ols$lat-data_ols$plat), 
           I((data_ols$lon-data_ols$plon)^2)/1000, I((data_ols$lat-data_ols$plat)^2)/1000, 
           I(log(data_ols$calcacres_n)), I(log((data_ols$ENGMeanScaleScore15_n)/100)), I(log((data_ols$median_income_n)/10000)), I(data_ols$age_n), 
           I(data_ols$age_n*data_ols$age_n)/1000, I(data_ols$totbath_n), I(data_ols$stories_n), I(data_ols$centheat_n), I(data_ols$fourthquart_n),
           I(data_ols$pct_renter_occupied_n), I(data_ols$sqft_n)/1000, I(data_ols$lon_n-data_ols$plon_n), I(data_ols$lat_n-data_ols$plat_n), 
           I((data_ols$lon_n-data_ols$plon_n)^2)/1000, I((data_ols$lat_n-data_ols$plat_n)^2)/1000)


ols_agg <- lm(Y~X-1);
summary(ols_agg)
deviance(ols_agg)

#weighted RMSE

Y <- I(datanew$salesprice)/100000

X <- cbind(1, I(log(datanew$calcacres)), I(log((datanew$ENGMeanScaleScore15)/100)), I(log((datanew$median_income)/10000)), I(datanew$age), 
           I(datanew$age*datanew$age)/1000, I(datanew$totbath), I(datanew$stories), I(datanew$centheat), I(datanew$fourthquart),
           I(datanew$pct_renter_occupied), I(datanew$sqft)/1000, I(datanew$lon-datanew$plon), I(datanew$lat-datanew$plat), 
           I((datanew$lon-datanew$plon)^2)/1000, I((datanew$lat-datanew$plat)^2)/1000)

ols_agg <- lm(Y~X-1);
pred <- predict(ols_agg, newdata = data.frame(X))
sqrt(mean((Y - pred*datanew$V1-pred*datanew$V2)^2))*100000






#standard errors

#this is original

par3 <- c(-2.109626e+05,  4.460316e-01,  1.270029e+01,  5.553325e-01,  3.686364e-03,  2.205276e-02,  8.017481e-01,  1.067674e-01, -4.144621e-03,  1.743981e-01,
          -4.281453e-03,  6.328598e-01, -4.762233e+03,  6.200320e+02, -2.830400e+04, -9.039645e+03, -2.109611e+05,  2.116792e+00,  1.354729e+01,  9.884247e-01,
          -2.360048e-02,  3.985402e-01,  1.470050e+00,  5.869073e-01, -1.061358e-01, -7.041625e-02, -5.112193e-03,  9.428571e-01, -4.762243e+03,  6.199476e+02,
          -2.830388e+04, -9.040922e+03,  1.176798e+00,  3.537552e+00)

#this is for bayesian (this works)

par3 <- c(-2.109596e+05,  4.686432e-01,  1.175717e+01,  5.572309e-01,  1.033144e-03,  4.809619e-05,  7.988188e-01,  1.340362e-02, -4.321235e-04,  1.722443e-01,
          -4.724813e-03,  7.169791e-01, -4.762159e+03,  6.200727e+02, -2.830370e+01, -9.038494e+00, -2.109625e+05,  2.078251e+00,  1.285127e+01,  1.076832e+00,
          -2.017951e-02,  3.609336e-04,  1.472500e+00,  6.867002e-01, -8.974575e-02, -5.156546e-02, -3.996879e-03,  9.197446e-01, -4.762208e+03,  6.201338e+02,
          -2.830390e+01, -9.041939e+00,  1.158014e+00,  3.573728e+00)


beta_opt <- optim(par3,FnTwo,d=d,x=X,y=Y,control=list(fnscale=-1,maxit=100000),hessian=TRUE)

b <- matrix(beta_opt$par[1:(j*niv)],niv,j); 
bse1 <- sqrt(-diag(solve(beta_opt$hessian[1:niv,1:niv])))
bse2 <- sqrt(-diag(solve(beta_opt$hessian[(niv+1):(2*niv),(niv+1):(2*niv)])))

s <- beta_opt$par[(j*niv+1):(j*(niv+1))]

gamma_opt <- optim(g,FnFour,d=d,z=Z,y=Y,control=list(fnscale=-1,maxit=100000),hessian=TRUE)
g <- gamma_opt$par
gse1 <- sqrt(-diag(solve(gamma_opt$hessian[1:gvs,1:gvs])))

par2 <- matrix(b,(niv*j),1);  
par2 <- c(par2,s)

LL <- FnTwo(par2,d=d,x=X,y=Y) + FnFour(g,d=d,z=Z,y=Y);
AIC <- -2*LL+2*niv

Ds=d;
beta=b;
bse=cbind(bse1,bse2);
gamma=cbind(g[1:gvs],g[(gvs+1):(2*gvs)]);
gse=cbind(gse1);

data$yhat <- ((1*beta[1,1]+data$calcacres*beta[2,1]+data$stories*beta[3,1]+
                 data$age*beta[4,1]+data$centheat*beta[5,1]+data$totbath*beta[6,1]+data$sqft*beta[7,1])*Ds[,1] +
                (1*beta[1,2]+data$calcacres*beta[2,2]+data$stories*beta[3,2]+
                   data$age*beta[4,2]+data$centheat*beta[5,2]+data$totbath*beta[6,2]+data$sqft*beta[7,2])*Ds[,2])
data$resid <- data$salesprice-data$yhat
rmse2 <- sqrt(mean((data$resid)^2))








