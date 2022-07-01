rm(list=ls())

set.seed(12345)
data <- read.csv("atlanta single family.csv")
data$sqft <- as.numeric(data$sqft)
data <- subset(data, sqft>0)

data$lon <- -1*data$lon
data$plon <- 0 #84.3563
data$plat <- 0 #33.7663

data %>% select(salesprice, calcacres, ENGMeanScaleScore15, median_income, age, totbath, stories, centheat, fourthquart, pct_renter_occupied,
                sqft, lon, plon, lat, plat, pct_white, pct_black, pct_collegeDegree, median_house_value, HHsize, tract_cover,
                total_crime_house) -> data
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
                       "pct_white_n", "pct_black_n", "pct_collegeDegree_n", "median_house_value_n", "HHsize_n", "tract_cover_n", 
                       "total_crime_house_n", "distance")


#fully endogenized

data <- newdata

data %>% select(salesprice, calcacres, ENGMeanScaleScore15, median_income, age, totbath, stories, centheat, fourthquart, pct_renter_occupied,
                sqft, lon, plon, lat, plat, pct_white, pct_black, pct_collegeDegree, median_house_value, HHsize, tract_cover, 
                total_crime_house) -> data
data <- data[complete.cases(data), ]

Y <- I(data$salesprice)/100000

#house attributes

X <- cbind(1, I(log(data$calcacres)), I(log((data$ENGMeanScaleScore15)/100)), I(log((data$median_income)/10000)), I(data$age), 
           I(data$age*data$age)/1000, I(data$totbath), I(data$stories), I(data$centheat), I(data$fourthquart), I(data$pct_renter_occupied),
           I(data$sqft)/1000, I(data$lon-data$plon), I(data$lat-data$plat), I((data$lon-data$plon)^2)/1000, I((data$lat-data$plat)^2)/1000,
           I(data$tract_cover), I(data$total_crime_house))

#demographic variables/mixing variables 

Z <- cbind(1,I(data$ENGMeanScaleScore15)/100, I(data$calcacres), I(data$median_income)/10000, I(data$pct_black), I(data$pct_renter_occupied),
           I(data$pct_collegeDegree), I(data$median_house_value)/100000, I(data$tract_cover), I(data$total_crime_house))

ols_agg <- lm(Y~X-1);
summary(ols_agg)

#starting values for the hedonic estimates/betas for each type i.e. for mixing algorithm

types <- 4;

beta_start <- matrix(ols_agg$coef,(types*ncol(X)),1);  

#starting values for the gamma estimates for the demographic variables

gamma_start <- matrix(0.01,(1*ncol(Z)),1);

#starting values for sigma
sigma_start <- matrix(sqrt(mean(ols_agg$residuals^2)),types,1)     

#collecting initializing values

val_start <- c(beta_start,gamma_start,sigma_start);

vals <- val_start;
types <- 4;

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
  L[,4] <- 1
  for (m in 1:(j-3))   
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
    for (i in 1:(j-3)) 
    { 
      L[,1] <- 0
      L[,i+1] <- 0
      L[,i+2] <- 0
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

#g <- mix$vals_fin[(j*niv+1):((j*(niv+gvs)-gvs))];
g <- mix$vals_fin[(j*niv+1):(length(par)-types)];

#s <- mix$vals_fin[-(1:(j*(niv+gvs)-gvs))];
s <- mix$vals_fin[-(1:(length(par)-types))];

b <- matrix(b,niv,j);

b1 <- matrix(b,(niv*j),1); 
par3 <- c(b1,s);


########## EXTRA ##########

datanew <- cbind(newdata, d)
datanew$prob <- colnames(datanew[, c(46:49)])[apply(datanew[, c(46:49)],1,which.max)]
datanew$prob <- as.numeric(datanew$prob)

data1 <- subset(datanew, prob==1)
data2 <- subset(datanew, prob==2)
data3 <- subset(datanew, prob==3)
data4 <- subset(datanew, prob==4)

data1 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, tract_cover, total_crime_house, pct_renter_occupied, type = "mean_sd") -> sumdata1
data2 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, tract_cover, total_crime_house, pct_renter_occupied, type = "mean_sd") -> sumdata2
data3 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, tract_cover, total_crime_house, pct_renter_occupied, type = "mean_sd") -> sumdata3
data4 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, tract_cover, total_crime_house, pct_renter_occupied, type = "mean_sd") -> sumdata4

LL <- FnTwo(par3,d=d,x=X,y=Y) + FnFour(g,d=d,z=Z,y=Y);
AIC <- -2*LL+2*niv

Ds=d;
beta=b;
#bse=cbind(bse1,bse2);
gamma=cbind(g[1:gvs],g[(gvs+1):(2*gvs)]);
#gse=cbind(gse1);

X <- cbind(1, I(log(data$calcacres)), I(log((data$ENGMeanScaleScore15)/100)), I(log((data$median_income)/10000)), I(data$age), 
           I(data$age*data$age)/1000, I(data$totbath), I(data$stories), I(data$centheat), I(data$fourthquart), I(data$pct_renter_occupied),
           I(data$sqft)/1000, I(data$lon-data$plon), I(data$lat-data$plat), I((data$lon-data$plon)^2)/1000, I((data$lat-data$plat)^2)/1000,
           I(data$tract_cover), I(data$total_crime_house))

pred <- (X%*%b[,1]*Ds[,1] + X%*%b[,2]*Ds[,2] + X%*%b[,3]*Ds[,3] + X%*%b[,4]*Ds[,4])*100000

resid <- Y-pred
rmse2 <- sqrt(mean((resid)^2))

########## EXTRA OLS ##########

prob <- read.csv("prob4.csv")
datanew <- cbind(newdata, prob)
datanew$prob <- colnames(datanew[, c(46:49)])[apply(datanew[, c(46:49)],1,which.max)]
datanew$prob <- ifelse(datanew$prob=="V1", 1, ifelse(datanew$prob=="V2", 2, ifelse(datanew$prob=="V3", 3, 4)))

data1 <- subset(datanew, prob==1)
data2 <- subset(datanew, prob==2)
data3 <- subset(datanew, prob==3)
data4 <- subset(datanew, prob==4)

data1 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, tract_cover, total_crime_house, pct_renter_occupied, type = "mean_sd") -> sumdata1
data2 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, tract_cover, total_crime_house, pct_renter_occupied, type = "mean_sd") -> sumdata2
data3 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, tract_cover, total_crime_house, pct_renter_occupied, type = "mean_sd") -> sumdata3
data4 %>% get_summary_stats(salesprice, ENGMeanScaleScore15, calcacres, median_income, age, pct_white, pct_black, 
                            pct_collegeDegree, median_house_value, HHsize, tract_cover, total_crime_house, pct_renter_occupied, type = "mean_sd") -> sumdata4

LL <- FnTwo(par3,d=d,x=X,y=Y) + FnFour(g,d=d,z=Z,y=Y);
AIC <- -2*LL+2*niv

##### OLS #####

data_ols <- subset(datanew, prob==4) #change prob==2

Y <- I(data_ols$salesprice)/100000

X <- cbind(1, I(log(data_ols$calcacres)), I(log((data_ols$ENGMeanScaleScore15)/100)), I(log((data_ols$median_income)/10000)), I(data_ols$age), 
           I(data_ols$age*data_ols$age)/1000, I(data_ols$totbath), I(data_ols$stories), I(data_ols$centheat), I(data_ols$fourthquart), I(data_ols$pct_renter_occupied),
           I(data_ols$sqft)/1000, I(data_ols$lon-data_ols$plon), I(data_ols$lat-data_ols$plat), I((data_ols$lon-data_ols$plon)^2)/1000, I((data_ols$lat-data_ols$plat)^2)/1000,
           I(data_ols$tract_cover), I(data_ols$total_crime_house))

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

pred <- predict(ols_agg, newdata = data.frame(X))
sqrt(mean((Y*100000 - pred*100000)^2))

#weighted RMSE

Y <- I(datanew$salesprice)/100000

X <- cbind(1, I(log(datanew$calcacres)), I(log((datanew$ENGMeanScaleScore15)/100)), I(log((datanew$median_income)/10000)), I(datanew$age), 
           I(datanew$age*datanew$age)/1000, I(datanew$totbath), I(datanew$stories), I(datanew$centheat), I(datanew$fourthquart),
           I(datanew$pct_renter_occupied), I(datanew$sqft)/1000, I(datanew$lon-datanew$plon), I(datanew$lat-datanew$plat), 
           I((datanew$lon-datanew$plon)^2)/1000, I((datanew$lat-datanew$plat)^2)/1000, I(datanew$tract_cover), I(datanew$total_crime_house))

ols_agg <- lm(Y~X-1);
pred <- predict(ols_agg, newdata = data.frame(X))
sqrt(mean((Y - pred*datanew$V1-pred*datanew$V2-pred*datanew$V3--pred*datanew$V4)^2))*100000










#standard errors

beta_opt <- optim(par3,FnTwo,d=d,x=X,y=Y,control=list(fnscale=-1,maxit=10000, ndeps=rep(0.00005, 28)),
                  hessian=TRUE, method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))
b <- matrix(beta_opt$par[1:(j*niv)],niv,j); 
bse1 <- sqrt(-diag(solve(beta_opt$hessian[1:niv,1:niv])))
bse2 <- sqrt(-diag(solve(beta_opt$hessian[(niv+1):(2*niv),(niv+1):(2*niv)])))
bse3 <- sqrt(-diag(solve(beta_opt$hessian[(niv*2+1):(3*niv),(niv*2+1):(3*niv)])))

s <- beta_opt$par[(j*niv+1):(j*(niv+1))]

gamma_opt <- optim(g,FnFour,d=d,z=Z,y=Y,control=list(fnscale=-1,maxit=100000),hessian=TRUE)
g <- gamma_opt$par
gse1 <- sqrt(-diag(solve(gamma_opt$hessian[1:gvs,1:gvs])))

par2 <- matrix(b,(niv*j),1);  
par2 <- c(par2,s)

LL <- FnTwo(par2,d=d,x=X,y=Y) + FnFour(g,d=d,z=Z,y=Y);
AIC <- -2*LL+2*niv

Ds=d; #probabilities of belonging in submarkets
beta=b; #beta estimates
bse=cbind(bse1,bse2,bse3); #standard errors
gamma=cbind(g[1:gvs],g[(gvs+1):(2*gvs)]);
gse=cbind(gse1);

row_nombre <- c("Intercept","Lot Size","Stories", "House Age", "Central Heat", "Total Bath")

write.table(b, file= paste0(path,"Beta.csv"),quote = FALSE, row.names= row_nombre, col.names=col_nombre, sep=",")
write.table(bse, file= paste0(path,"Bse.csv"),quote = FALSE, row.names= row_nombre, col.names=col_nombre, sep=",")
write.table(LL, file= paste0(path,"LL.csv"),quote = FALSE, row.names= TRUE, col.names=TRUE, sep=",")
write.table(s, file= paste0(path,"S.csv"),quote = FALSE, row.names= TRUE, col.names=TRUE, sep=",")
write.table(gse1, file= paste0(path,"Gse.csv"),quote = FALSE, row.names= TRUE, col.names=TRUE, sep=",")
write.table(gamma, file= paste0(path,"Gamma.csv"),quote = FALSE, row.names= TRUE, col.names=TRUE, sep=",")
write.table(d, file= paste0(path,"Dhats.csv"),quote = FALSE, row.names= TRUE, col.names=col_nombre, sep=",")