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


#fully endogenized

set.seed(12345)
data <- read.csv("atlanta single family.csv")
data$sqft <- as.numeric(data$sqft)
data <- data[complete.cases(data), ]
data$lon <- -1*data$lon
data$plon <- 84.3720
data$plat <- 33.7885

data %>% select(X, salesprice, calcacres, stories, age, centheat, totbath, ENGMeanScaleScore15, median_income, sqft4,
                median_house_value, HHsize, pct_black, pct_white) -> data
data <- data[complete.cases(data), ]

Y <- I(data$salesprice)/100000

#house attributes

X <- cbind(1, I(log(data$calcacres)), I(log((data$ENGMeanScaleScore15)/100)), I(log((data$median_income)/10000)), I(data$age), 
           I(data$age*data$age), I(data$totbath), I(data$stories), I(data$centheat), I(data$fourthquart), I(data$pct_renter_occupied),
           I(data$sqft)/1000, I(data$lon-data$plon), I(data$lat-data$plat), I((data$lat-data$plat)^2), I((data$lon-data$plon)^2))

#demographic variables/mixing variables

Z <- cbind(1,I(data$ENGMeanScaleScore15)/100, I(data$calcacres), I(data$median_income)/10000, I(data$age), I(data$pct_white), I(data$pct_black),
           I(data$pct_collegeDegree))











