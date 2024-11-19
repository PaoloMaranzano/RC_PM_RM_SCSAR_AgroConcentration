#########################################################################################################################################
##########                              Cerqueti, R., Maranzano, P. & Mattera, R. (2024+)                                      ##########
########## "Spatially-clustered spatial autoregressive models with application to agricultural market concentration in Europe" ##########
##########                                  Code for the application (Sections 4 and 5)                                        ##########
#########################################################################################################################################

#################
##### Setup #####
#################
library(SparseM)
library(MASS)
library(spatialreg)
library(spdep)
library(sp)
library(dplyr)
library(spdep)
library(spatialreg)
library(sf)
library(tidyverse)
library(viridis)

rm(list=ls())



#########################################
########## Auxiliary functions ##########
#########################################
setwd("C:/Users/paulm/OneDrive/Documenti/GitHub/RC_PM_RM_SCSAR_AgroConcentration/Data and code")


#################################################################################################
########## Dataset for the analysis of the agricultural market concentration in Europe ##########
#################################################################################################

##### Dataset from Eurostat
load("RC_PM_RM_SCSAR_Data.RData")

##### Select aggregate data for 2010 and 2020
Data <- Data %>%
  dplyr::filter(time %in% c(2010,2020),
                farmtype_lab %in% c("Total"), uaarea_lab %in% c("Total"),
                organic_lab %in% c("Total"),so_eur %in% c("KE0")) %>%
  mutate(
    ##### Output distribution (Gini on standard output)
    Gini_SO = Gini_SO,
    ##### Wealth
    # Per capita GDP: overall regional wealth
    GDPPC_PPS2020 = GDPPC_PPS2020,
    ##### Labor market
    # Share of employment in agriculture: relevance of agricultural industry on the regional labor market
    Share_AgroEmp = EmpPersThs_NACE_A/EmpPersThs_NACE_Total*100,
    # Hours worked per agro-employed: agricultural labor market intensity
    HoursWorked_AgroEmp = WorkedHoursThs_NACE_A/EmpPersThs_NACE_A,
    ##### Economy
    # Gross value added per agro-employed: agricultural productivity intensity
    GVA_AgroEmp = GVA_MillEuro_NACE_A/EmpPersThs_NACE_A,
    # Investment per agro-employed: propensity to invest according to the economic size
    GFCF_AgroEmp = GFCF_MillEuro_NACE_A/EmpPersThs_NACE_A,
    # Share of agricultural GVA on total GVA: relevance of agricultural industry on the regional economy
    Share_AgroGVA = GVA_MillEuro_NACE_A/GVA_MillEuro_NACE_Total*100,
    ##### Landscape and environment
    # Share of agricultural land: relevance of agricultural industry on the regional activities
    Share_AgroLand = LandUseKM2_Agro/LandUseKM2_Total*100,
    # Average altitude: geography and landscape
    Alt_mean = Alt_mean,
    # Heating degree days (HDD): proxy of temperature and weather conditions
    HDD = HDD) %>%
  select(Year = time, geo_lab, geo, Gini_SO, GDPPC_PPS2020, Share_AgroEmp,
         HoursWorked_AgroEmp, GVA_AgroEmp, GFCF_AgroEmp, Share_AgroGVA, Share_AgroLand,
         Alt_mean, HDD)

##### Select complete rows
Data <- Data[!st_is_empty(Data$geometry), ]
Data <- na.omit(Data)

##### Split 2010 and 2020
Data2010 <- Data %>%
  filter(Year == 2010)
Data2020 <- Data %>%
  filter(Year == 2020)

### Sp and W matrix
sp_object <- as(Data2020, "Spatial")
nb <- poly2nb(sp_object)
listW <- nb2listw(nb, style = "W", zero.policy=TRUE)

##### Export data
save(Data,Data2010,Data2020,Countries,listW,file = "Data_RC_PM_RM_JABES2024.RData")
save(Data2010,Data2020,listW,file = "Data_RC_PM_RM_JABES2024.rda")
