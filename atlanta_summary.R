#https://rural-urban.eu/sites/default/files/05_Spatial%20Autocorrelation%20and%20the%20Spatial%20Durbin%20Model_Eilers.pdf
#https://d-nb.info/1188635956/34
#http://rstudio-pubs-static.s3.amazonaws.com/5027_52298866e7924b18b54e5c9a0a21b450.html
#https://stats.stackexchange.com/questions/149415/how-do-i-interpret-lagsarlm-output-from-rs-spdep
#https://rpubs.com/corey_sparks/109650
#https://rpubs.com/corey_sparks/108130
#https://www.stata.com/training/webinar_series/spatial-autoregressive-models/spatial/resource/spillover.html

library(dplyr)
library(sp)
library(rgeos)
library(geosphere)
library(ggplot2)
library(sf)

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

newdata$price_sqft <- newdata$salesprice/newdata$sqft

sub <- subset(newdata, dec<=9)

newdata$dec <- ntile(newdata$salesprice,10)
newdata %>% group_by(dec) %>% summarise(sqft = mean(sqft),
                                        price_sqft = mean(price_sqft))


#crime and canopy cover

#https://www.atlantapd.org/i-want-to/crime-data-downloads
#https://www.mrlc.gov/data/nlcd-2016-usfs-tree-canopy-cover-conus
#https://shiandy.com/post/2020/11/02/mapping-lat-long-to-fips/

library(tigris)
library(leaflet)
library(tidyr)

data <- read.csv("atlanta single family.csv")
data$sqft <- as.numeric(data$sqft)
data <- subset(data, sqft>0)

atl_roads <- roads("GA", "Fulton")

ggplot(atl_roads) + 
  geom_sf() + 
  theme_void()


crime <- read.csv("COBRA-2015-16.csv") #number of incidents

crime <- crime[!is.na(crime$Latitude),]
crime <- crime[!is.na(crime$Longitude),]

coord <- dplyr::select(crime, Latitude, Longitude)
coord <- coord[complete.cases(coord), ]

tarrant <- tracts("GA", "Fulton", cb = TRUE)

leaflet(tarrant) %>%
  addTiles() %>%
  addPolygons(popup = ~NAME)

coord_sf <- coord %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(tarrant))

intersected <- st_intersects(coord_sf, tarrant)

latlong_final <- coord_sf %>%
  mutate(intersection = as.integer(intersected),
         fips = if_else(is.na(intersection), "",
                        tarrant$TRACTCE[intersection]))

crime <- cbind(crime, latlong_final)
crime$fips2 <- as.numeric(crime$fips)

crime$crime <- ifelse(crime$UCR.Literal == "BURGLARY-RESIDENCE" | crime$UCR.Literal == "LARCENY-FROM VEHICLE" | 
                      crime$UCR.Literal == "LARCENY-NON VEHICLE" | crime$UCR.Literal == "ROBBERY-PEDESTRIAN" | 
                      crime$UCR.Literal == "AUTO THEFT" | crime$UCR.Literal == "ROBBERY-RESIDENCE" | 
                      crime$UCR.Literal == "ROBBERY-COMMERCIAL" | crime$UCR.Literal == "BURGLARY-NONRES", "PROP", "VIOLENT")

crime <- separate(crime, col=Possible.Date, into=c('Month', 'Date', 'Year'), sep='/')

crime$Month <- as.numeric(crime$Month)
crime$Date <- as.numeric(crime$Date)
crime$Year <- as.numeric(crime$Year)

crime <- crime[!is.na(crime$fips2),]

crime_cat = crime %>% group_by(Year, fips2)  %>%
  summarise(total_prop = sum(crime=="PROP"),
            total_violent = sum(crime=="VIOLENT"))

crime_cat <- subset(crime_cat, Year>=2015)

crime_cat$GEOID <- paste0(crime_cat$Year, crime_cat$fips2)

data$tract <- data$tract2010*100
data$GEOID <- paste0(data$saledtyearyyyy, data$tract)

#tree canopy by individual houses

library(raster)

raster1 <- raster("nlcd_2016_treecanopy_2019_08_31.img")

coord <- dplyr::select(data, lat, lon)

str(coord)
coordinates(coord)  <-  c("lon",  "lat")
proj4string(coord)  <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
sites_transformed<-spTransform(coord, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

indiv_treecover <- raster::extract(raster1, sites_transformed)

data$indiv_treecover <- indiv_treecover

#the reason I'm merging it here instead of after 'data$GEOID <- paste0(data$saledtyearyyyy, data$tract)' is bc 
#merging before removes 2 obs rendering it impossible to attach the indiv_treecover
jointdataset <- merge(data, crime_cat, by = 'GEOID') 

jointdataset$total_crime <- jointdataset$total_prop+jointdataset$total_violent

#tree canopy by individual houses

tarrant <- tracts("GA", "Fulton", cb = TRUE)

tarrant$geometry %>% st_set_crs(4326) %>% st_transform(crs=32736) -> new_geometry

exactextractr::exact_extract(raster1, new_geometry) -> geom

lapply(geom,rowMeans) -> tbl

lapply(tbl,mean) -> tbl2

tbl3 <- unlist(tbl2)

tarrant$tract_cover <- tbl3
tarrant$tract <- as.numeric(tarrant$TRACTCE)
tarrant2 <- tarrant
tarrant2$geometry <- NULL

tarrant2 <- dplyr::select(tarrant2, tract_cover, tract)

jointdataset2 <- merge(jointdataset, tarrant2, by = 'tract')

#getting tractwise population and houses

library(tidycensus)
library(tidyverse)
library(maps)
library(sf)
library(data.table)

vt1 <- get_acs(geography = "tract", 
               variables = c(hh_number="DP02_0001", tract_pop="B01003_001"), 
               state="GA", county="Fulton",
               year = 2010)

vt1_dt <- data.table::as.data.table(vt1)
data_wide_vt1 <- dcast(vt1_dt, GEOID ~ variable, fun.aggregate = mean,
                       value.var=c("estimate"))

data_wide_vt1$tract = substring(data_wide_vt1$GEOID, 6)
data_wide_vt1$tract <- as.numeric(data_wide_vt1$tract)
data_wide_vt1$GEOID <- NULL

jointdataset3 <- merge(jointdataset2, data_wide_vt1, by = 'tract')

jointdataset3$total_prop_house <- jointdataset3$total_prop/jointdataset3$hh_number
jointdataset3$total_violent_house <- jointdataset3$total_violent/jointdataset3$hh_number
jointdataset3$total_crime_house <- jointdataset3$total_crime/jointdataset3$hh_number

jointdataset3$total_prop_pop <- jointdataset3$total_prop/jointdataset3$tract_pop
jointdataset3$total_violent_pop <- jointdataset3$total_violent/jointdataset3$tract_pop
jointdataset3$total_crime_pop <- jointdataset3$total_crime/jointdataset3$tract_pop

#regression testing

jointdataset3$lon <- -1*jointdataset3$lon
jointdataset3$plon <- 0 #84.3563
jointdataset3$plat <- 0 #33.7663

jointdataset3 %>% dplyr::select(salesprice, calcacres, ENGMeanScaleScore15, median_income, age, totbath, stories, centheat, fourthquart, pct_renter_occupied,
                sqft, lon, plon, lat, plat, pct_white, pct_black, pct_collegeDegree, median_house_value, HHsize, indiv_treecover, tract_cover,
                total_crime_house) -> data

data <- data[complete.cases(data), ]

Y <- I(data$salesprice)/100000

#house attributes

X <- cbind(1, I(log(data$calcacres)), I(log((data$ENGMeanScaleScore15)/100)), I(log((data$median_income)/10000)), I(data$age), 
           I(data$age*data$age)/1000, I(data$totbath), I(data$stories), I(data$centheat), I(data$fourthquart), I(data$pct_renter_occupied),
           I(data$sqft)/1000, I(data$lon-data$plon), I(data$lat-data$plat), I((data$lon-data$plon)^2)/1000, I((data$lat-data$plat)^2)/1000, 
           I(data$tract_cover), I(data$total_crime_house))

#demographic variables/mixing variables 

Z <- cbind(1,I(data$ENGMeanScaleScore15)/100, I(data$calcacres), I(data$median_income)/10000, I(data$age), I(data$pct_white), I(data$pct_black),
           I(data$pct_collegeDegree), I(data$median_house_value)/100000, I(data$HHsize))

ols_agg <- lm(Y~X-1);
summary(ols_agg)






















coord <- dplyr::select(tarrant, geometry)

str(coord)
coordinates(coord)  <-  c("lon",  "lat")
proj4string(coord)  <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
sites_transformed<-spTransform(coord, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

raster::extract(raster1, sites_transformed)




#xy <- cbind(c(-83.25, -82.75), c(39.8, 40.2))
#project(xy, "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#raster1 <- raster("C:/Users/syedm/Downloads/usfs_carto_CONUS_2016/usfs_2016_treecanopy_cartographic_12-14-2018.img")
xy <- cbind(c(min(data$lon), max(data$lon)), c(min(data$lat), max(data$lat)))
new_coord <- proj4::project(xy, "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
e <- extent(new_coord[1,1], new_coord[2,1], new_coord[1,2], new_coord[2,2])
#e <- extent(1045954, 1145839, 1919124, 2006382)
rc <- raster::crop(raster1, e)
plot(rc)
points(sites_transformed, col=data$price)

test_df <- as.data.frame(rc, xy=TRUE)
colnames(test_df) <- c("x", "y", "value")
test_df$value <- ifelse(test_df$value>=50, 50, test_df$value)

r <- raster(nrow=1673, ncol=1194, ext=extent(1078695, 1114515, 1937025, 1987215), crs="+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
values(r) <- test_df$value
plot(r)








