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

#example with ethiopia

set.seed(12345)
data <- read.csv("atlanta single family.csv")
data$sqft <- as.numeric(data$sqft)
data %>% select(salesprice, calcacres, stories, age, centheat, totbath, ENGMeanScaleScore15, median_income, sqft,
                median_house_value, HHsize, lat, lon, pct_black) -> data
data$lon <- -1*data$lon
data <- data[complete.cases(data), ]

data_sp <- data

coordinates(data_sp) <- ~lon+lat

class(data_sp) #check if SpatialPointsDataFrame
d <- distm(data_sp)

#d <- gDistance(data_sp, byid=T)
min.d <- apply(d, 1, function(x) order(x, decreasing=F)[2])

newdata <- cbind(data, data[min.d,], apply(d, 1, function(x) sort(x, decreasing=F)[2]))

colnames(newdata) <- c(colnames(data), "salesprice_n", "calcacres_n", "stories_n", "age_n", "centheat_n", "totbath_n", 
                       "ENGMeanScaleScore15_n", "median_income_n", "sqft_n", "median_house_value_n", "HHsize_n", "lat_n", 
                       "lon_n", "pct_black_n", "distance")








data$c <- paste(data$lat, data$lon, sep="_")
data <- distinct(data, c, .keep_all = TRUE)
data$c <- NULL

spatialreg::set.ZeroPolicyOption(TRUE)
spdep::set.ZeroPolicyOption(TRUE)

#getting neighbors
coor <- cbind(data$lon, data$lat)
k1 <- knn2nb(knearneigh(coor, k=1))
k1dists <- unlist(nbdists(k1, coor, longlat=TRUE))
summary(k1dists) #tells you summary stats for the distance among houses

#cross checking

1
39.97604	-82.916504
0.6977135 -1.447166
10517
39.976151	-82.916496
0.6977155 -1.447166
10441
39.975422	-82.916557
0.6977027 -1.447167

deg2rad <- function(deg) return(deg*pi/180)

gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

#fixed weight matrix
fixed <- dnearneigh(coor, 0, 2, longlat=TRUE) #only houses within 0.5km #generates some households without neighbors

#adaptive weight matrix #my preferred method
adap <- knn2nb(knearneigh(coor, k=25)) #only nearest 25 houses

#inverse distance matrix using adaptive matrix
#if you want to use fixed weight matrix, replace 'adap' with 'fixed'
dist <- nbdists(fixed, coor, longlat=TRUE)

#some houses have the same coordinates, so the distance between them is 0 
#(this will potentially ruin the inverse distance calculation because you cant divide anything by 0)
#so i replaced all 0 distances with 0.005. you can play around with this number
dist2 <- lapply(dist, function(x) ifelse(x==0, 0.0005, x)) 

invd <- lapply(dist2, function(x) (1/(x/100)))
invd.w <- nb2listw(fixed, glist=invd, style="W", zero.policy=TRUE)
#invd.w <- nb2listw(fixed, glist=dist2, style="W", zero.policy=TRUE)
invd.w$weights[1]
invd.w #this is the weight matrix

#converting to matrix for better visibility
invd.m <- listw2mat(invd.w)

t <- apply(invd.m,2,function(x) sum(x > 0))

#model

eq1 <- tranamt ~ onestory + rooms + fullbath + agehouse + buildingsqft + lotsize + 
  pctblack + college + mediany_cbg + science10

eq1 <- salesprice ~ stories + age + centheat + totbath + ENGMeanScaleScore15 + median_income + 
  sqft + median_house_value + HHsize + calcacres

reg1 <- lm(eq1, data=data)
summary(reg1)

#first step

#model selection via lm tests
lmtest <- lm.LMtests(reg1, invd.w, test=c("LMerr", "LMlag", "RLMerr", "RLMlag",
                                          "SARMA"))
lmtest #check results of lm tests

#from lm tests, first compare lmerr and lmlag model. pick the one that is statistically significant 
#by that, i mean if lmerr is significant, proceed to running the spatial error model 
#otherwise, run the spatial lag model
#if both are statistically significant, compare between rlmerr and rlmlag and pick the one that is statistically significant
#if rlmerr is significant, proceed to running the spatial error model 
#otherwise run the spatial lag model

#if none are significant, we go with ols (nothing left to do)

#spatial error model
reg2 <- errorsarlm(eq1, data=data, invd.w, tol.solve=1.0e-30)
summary(reg2) #results of spatial model
#summary(impacts(reg2, listw=invd.w, R=500), zstats = TRUE) #impact 

#spatial lag model
reg3 <- lagsarlm(eq1, data = data, invd.w, tol.solve=1.0e-30)
summary(reg3)
summary(impacts(reg3, listw=invd.w, R=500), zstats = TRUE)

#second step

#run spatial durbin model

#spatial durbin model

reg4 <- lagsarlm(eq1, data = data, invd.w, Durbin=TRUE, tol.solve=1.0e-30)
summary(reg4)
summary(impacts(reg4, listw=invd.w, R=500, Q=5),zstats=TRUE)

#third step
#if you picked spatial error model in step 1, compare spatial durbin with spatial error model

#compare spatial durbin with spatial error
LR.sarlm(reg2, reg4)

#if test is statistically significant, reject null and pick spatial durbin model
#if not, the spatial error model is our final model

#if we pick spatial durbin model, it may make sense to compare with spatially lagged x model

#run spatially lagged x model

#spatially lagged x model
reg5 = lmSLX(eq1, data = data, invd.w)
summary(reg5)
summary(impacts(reg5, listw=invd.w, R=500), zstats = TRUE)

#now compare spatial durbin with spatial x model

#compare spatial durbin with spatial x
LR.sarlm(reg4, reg5)

#if statistically significant, durbin is our final model. if not, spatial x is our final model

#if you picked spatial lag model in step 1, compare spatial durbin with spatial lag model

#compare spatial durbin with spatial lag
LR.sarlm(reg3, reg4)

#if test is statistically significant, reject null and pick spatial durbin model
#if not, the spatial lag model is our final model

#if we pick spatial durbin model, it again may make sense to compare with spatially lagged x model

#compare spatial durbin with spatial x
LR.sarlm(reg4, reg5)

#if statistically significant, durbin is our final model. if not, spatial x is our final model

#mapping

library(readxl)
library(splm)
library(rgdal)
library(spdep)
library(plm)
library(lmtest)
library(fields)
library(tidyverse)
library(sf)
library(leaflet)
library(mapview)
library(ggmap)
library(ggridges)
library(viridis)
library(rayshader)



