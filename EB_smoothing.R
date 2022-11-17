##EB Smoothing for Sero-prevalence estimates ###
library(dplyr)
library(readxl)
library(writexl)
library(janitor)
library(tidyverse)
library(expss)
library(survey)
library(sp)
library(sf)
library(spdep) #EB smoothing

#Read in dataset
geoSeroprev <- st_read("your_spatial_file.gpkg")  #merged seroprevalence data by LGA in QGIS v. 3.26
class(geoSeroprev) #sf

options(scipen= 999) #prevent numbers from coming up as scientific notation

#Queen neighbors
queen_nb <- poly2nb(geoSeroprev, queen= TRUE)
summary(queen_nb)

### Visualizing the neighbors
# Create a matrix of the x,y coordinates for each LGA centroid
NG_cent <- st_centroid(st_geometry(geoSeroprev))

# Plot the outlines of LGA with light grey boundaries
plot(st_geometry(geoSeroprev), border = 'grey')
# Add the plot of the Queen contiguity connections
plot.nb(queen_nb, NG_cent, points = F, add = T)

#### Estimate spatial (local) EB under the Queen contiguity neighbor definition ####
eb_Pm <- EBlocal(geoSeroprev$n_PmMSP1pos, geoSeroprev$n, nb = queen_nb)
eb_Po <- EBlocal(geoSeroprev$n_PoMSP1pos, geoSeroprev$n, nb = queen_nb)
eb_Pv <- EBlocal(geoSeroprev$n_PvMSP1pos, geoSeroprev$n, nb = queen_nb)
eb_Pf <- EBlocal(geoSeroprev$n_PfMSP1pos, geoSeroprev$n, nb = queen_nb)

# The output from EBlocal() is a 2 column data frame. The second column is the EB estimate. 
geoSeroprev$EB_Pm <- eb_Pm[,2]
geoSeroprev$EB_Po <- eb_Po[,2]
geoSeroprev$EB_Pv <- eb_Pv[,2]
geoSeroprev$EB_Pf <- eb_Pf[,2]

setwd("yourworkingdirectory")
#Save the gpkg file
st_write(geoSeroprev, "EB_seroprev_by_LGA.gpkg", driver= "GPKG")
