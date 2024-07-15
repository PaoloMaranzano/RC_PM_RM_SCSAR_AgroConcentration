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
library(ggplot2)
library(tidyverse)
library(texreg)
library(xtable)
library(ggpubr)

rm(list=ls())

'%notin%' <- Negate('%in%')
ToyExample <- FALSE
setwd("C:/Users/paulm/OneDrive/Documenti/GitHub/RC_PM_RM_SCSAR_AgroConcentration/Data and code")


#########################################
########## Auxiliary functions ##########
#########################################
source("RC_PM_RM_SCSAR_AuxFuns", encoding = 'UTF-8')


#################################
########## Toy example ##########
#################################

if (ToyExample == TRUE) {
  # Example data
  Data <- rio::import("Matrice_SLL_96_05_19.xlsx", sheet=3) # 1 for 1996, 3 for 2019
  Data$vaa <- log(Data$vaa)
  Data$vap <- log(Data$vap)
  Data$vas <- log(Data$vas)
  Data$den <- log(Data$den)
  Data$dis <- log(Data$dis)
  
  Data[,3:10] <- scale(Data[,3:10])
  
  Y0 <- cbind(Data$vaa)
  X0 <- Data[,6:ncol(Data)]
  X0 <- X0[,-6] # rimuovo "ind" perchè crea clusters dove ind è sempre 0
  # attenzione con dummy variables, queste cose possono capitare!!
  
  shpprov <- readOGR("SLL_2011_2018.shp") # 2019
  rswm_q <- poly2nb(shpprov, queen = TRUE)
  set.ZeroPolicyOption(TRUE)
  set.ZeroPolicyOption(TRUE)
  Bmat <- nb2listw(rswm_q, style = "B", zero.policy = T)
  Wmat2 <- listw2mat(Bmat)
  Sp <- coordinates(shpprov)
  #Wmat1 <- as.matrix(dist(Sp))
  
  # Do methods work?
  
  
  
  #############################
  
  reg0 <- Cl_spatialReg(Y=Y0, X=X0, Sp=Sp, W=Wmat2, type="lm")
  summary(reg0$Results[[1]])
  summary(reg0$Results[[2]])
  
  reg0 <- Cl_spatialReg(Y=Y0, X=X0, G=2, Sp=Sp, W=Wmat2, type="lagsarlm")
  summary(reg0$Results[[1]])
  summary(reg0$Results[[2]])
  reg0$group
  table(reg0$group)
  reg0$ML
  
  a=G_select(Y=Y0, X=X0, Sp=Sp, W=Wmat2, type="lagsarlm")
  
  reg0 <- Cl_spatialReg(Y=Y0, X=X0, Sp=Sp, W=Wmat2, type="errorsarlm")
  summary(reg0$Results[[1]])
  summary(reg0$Results[[2]])
  reg0$group
  table(reg0$group)
  reg0$ML
  
  reg0 <- Cl_spatialReg(Y=Y0, X=X0, Sp=Sp, W=Wmat2, type="lmSLX")
  summary(reg0$Results[[1]])
  summary(reg0$Results[[2]])
  reg0$group
  table(reg0$group)
  reg0$ML
}



##############################################################################################
########## Application: Analysis of the agricultural market concentration in Europe ##########
##############################################################################################

##### Dataset from Eurostat
load("RC_PM_RM_SCSAR_Data.RData")

##### Graphical representation
Data$RoadsKM_KM2 <- Data$RoadsKM/Data$LandKM2_Total
Y_name <- "RoadsKM_KM2"
Data %>%
  select(time,Y = .data[[Y_name]]) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Y)) + 
  geom_sf(data = Countries, linewidth = 1.05, show.legend=FALSE, alpha=0,color="#000000") + 
  facet_wrap(~time) + 
  scale_fill_gradient2(mid = "#FFFFFF",low = "#00FF00",high = "#FF0000",
                       midpoint = mean(Data[[Y_name]],na.rm=T)) + 
  labs(x = "", y = "", title = Y_name) + 
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 12),strip.text = element_text(size = 16),
        axis.title = element_text(size = 14))





######################################################
##### Parameters tuning with likelihood criteria #####
######################################################

# 2010
G_selection_2010_1 <- G_select(Y=dta2010b_non_empty$Gini_SO, X = Xmat_2010,Sp = Sp, W = Wmat, G.set=2:4, type="lagsarlm",
                               maxitr = 50,Phi = 1)
G_selection_2010_75 <- G_select(Y=dta2010b_non_empty$Gini_SO, X = Xmat_2010,Sp = Sp, W = Wmat, G.set=2:4, type="lagsarlm",
                                maxitr = 50,Phi = 0.75)
G_selection_2010_50 <- G_select(Y=dta2010b_non_empty$Gini_SO, X = Xmat_2010,Sp = Sp, W = Wmat, G.set=2:4, type="lagsarlm",
                                maxitr = 50,Phi = 0.50)

# 2020
G_selection_2020_1 <- G_select(Y=dta2020b_non_empty$Gini_SO, X = Xmat_2020,Sp = Sp, W = Wmat, G.set=2:4, type="lagsarlm",
                               maxitr = 50,Phi = 1)
G_selection_2020_75 <- G_select(Y=dta2020b_non_empty$Gini_SO, X = Xmat_2020,Sp = Sp, W = Wmat, G.set=2:4, type="lagsarlm",
                                maxitr = 50,Phi = 0.75)
G_selection_2020_50 <- G_select(Y=dta2020b_non_empty$Gini_SO, X = Xmat_2020,Sp = Sp, W = Wmat, G.set=2:4, type="lagsarlm",
                                maxitr = 50,Phi = 0.50)

TuningBIC <- data.frame(rbind(cbind(Year = 2010, phi = 1, G = 2:4, BIC = G_selection_2010_1$BIC),
                              cbind(Year = 2010, phi = 0.75, G = 2:4, BIC = G_selection_2010_75$BIC),
                              cbind(Year = 2010, phi = 0.50, G = 2:4, BIC = G_selection_2010_50$BIC),
                              cbind(Year = 2020, phi = 1, G = 2:4, BIC = G_selection_2020_1$BIC),
                              cbind(Year = 2020, phi = 0.75, G = 2:4, BIC = G_selection_2020_75$BIC),
                              cbind(Year = 2020, phi = 0.50, G = 2:4, BIC = G_selection_2020_50$BIC)))
# save(TuningBIC,file = "TuningBIC.RData")
p <- TuningBIC %>%
  mutate(phi = as.factor(phi)) %>%
  ggplot(mapping = aes(x = as.factor(G), y = BIC, group = phi)) + 
  geom_line(mapping = aes(col = phi), size = 1.5) + 
  geom_point(mapping = aes(col = phi), size = 4) + 
  facet_wrap(~ Year) + 
  theme_bw() + 
  # scale_y_continuous(breaks = seq(1400,1800,500)) +
  labs(x = "Number of clusters (K)",
       title = latex2exp::TeX("Tuning of the spatial penalty parameter ($\\phi$) and number of clusters ($\\K$)"))
ggexport(p,width = 1800, height = 1200, res = 150, filename = paste0("ParsTuning.png"))





################################
##### SCSAR model for 2010 #####
################################
dta2010 <- Data %>% 
  dplyr::filter(time == 2010,
                farmtype_lab %in% c("Total"), uaarea_lab %in% c("Total"),
                organic_lab %in% c("Total"),so_eur %in% c("KE0"))
dta2010b <- na.omit(dta2010)

### Sp and W matrix
non_empty_geometries <- !st_is_empty(dta2010b$geometry)
dta2010b_non_empty <- dta2010b[non_empty_geometries, ]
sp_object <- as(dta2010b_non_empty, "Spatial")
nb <- poly2nb(sp_object)
listw <- nb2listw(nb, style = "B", zero.policy=TRUE)
Wmat <- listw2mat(listw)
Sp <- coordinates(sp_object)

### Covariates
Xmat_2010 <- cbind(
  ##### Wealth
  # Per capita GDP: overall regional wealth
  dta2010b_non_empty$GDPPC_PPS2020,
  ##### Labor market
  # Share of employment in agriculture: relevance of agricultural industry on the regional labor market
  dta2010b_non_empty$EmpPersThs_NACE_A/dta2010b_non_empty$EmpPersThs_NACE_Total*100,
  # Hours worked per agro-employed: agricultural labor market intensity
  dta2010b_non_empty$WorkedHoursThs_NACE_A/dta2010b_non_empty$EmpPersThs_NACE_A,
  ##### Economy
  # Gross value added per agro-employed: agricultural productivity intensity
  dta2010b_non_empty$GVA_MillEuro_NACE_A/dta2010b_non_empty$EmpPersThs_NACE_A,
  # Investment per agro-employed: propensity to invest according to the economic size
  dta2010b_non_empty$GFCF_MillEuro_NACE_A/dta2010b_non_empty$EmpPersThs_NACE_A,
  # Share of agricultural GVA on total GVA: relevance of agricultural industry on the regional economy
  dta2010b_non_empty$GVA_MillEuro_NACE_A/dta2010b_non_empty$GVA_MillEuro_NACE_Total*100,
  ##### Landscape and environment
  # Share of agricultural land: relevance of agricultural industry on the regional activities
  dta2010b_non_empty$LandUseKM2_Agro/dta2010b_non_empty$LandUseKM2_Total*100,
  # Average altitude: geography and landscape
  dta2010b_non_empty$Alt_mean,
  # Heating degree days (HDD): proxy of temperature and weather conditions
  dta2010b_non_empty$HDD
  )

### Estimate the models
K_opt <- 3
reg0 <- Cl_spatialReg(Y=dta2010b_non_empty$Gini_SO, X=Xmat_2010, Sp=Sp, G=K_opt, W=Wmat, type="lagsarlm", Phi = 0.50)
regnull <- spatialreg::lagsarlm(dta2010b_non_empty$Gini_SO ~ Xmat_2010, listw=listw, zero.policy = T)
# Change names to coefficients
for (m in 1:K_opt) {
  names(reg0$Results[[m]]$coefficients) <- c("Intercept", "GDP pc",
                                             "Share of agro employment","Worked hours per agro worker",
                                             "GVA per agro worker","GFCF per agro worker","Share of agro GVA",
                                             "Share of agro land","Average altitude","HDD")
}
names(regnull$coefficients) <- c("Intercept", "GDP pc",
                                 "Share of agro employment","Worked hours per agro worker",
                                 "GVA per agro worker","GFCF per agro worker","Share of agro GVA",
                                 "Share of agro land","Average altitude","HDD")
### Store in latex
reglist <- c(list(regnull),reg0$Results[1:K_opt])
texreg_table <- texreg(l = reglist,
                       file = paste0("ResClusterAgro_K",K_opt,"_2010.tex"),
                       caption = paste0("Estimated pooled and cluster-wise coefficients for 2010 with K = ",K_opt),
                       custom.model.names = c("Pooled",paste0("G=",1:K_opt)))
texreg(l = reglist)
table(reg0$group)
dta2010_clust <- dta2010b_non_empty
dta2010_clust$Cl_hat <- as.factor(reg0$group)
colnames(dta2010_clust)[which(colnames(dta2010_clust) == "Cl_hat")] <- paste0("clust_K", K_opt)

### Mapping
p1 <- dta2010_clust %>% 
  select(time,'Clustering' = .data[[paste0("clust_K", K_opt)]]) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Clustering)) + 
  geom_sf(data = Countries, linewidth = 1.05, show.legend=FALSE, alpha=0,color="#000000") + 
  labs(x = "", y = "",  title = "Estimated cluster structure") +
  theme_bw() + 
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),strip.text = element_text(size = 16),
        axis.title = element_text(size = 14))
p2 <- dta2010_clust %>% 
  select(time,Gini = Gini_SO) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Gini)) + 
  geom_sf(data = Countries, linewidth = 1.05, show.legend=FALSE, alpha=0,color="#000000") + 
  scale_fill_gradient2(mid = "#FFFFFF",low = "#00FF00",high = "#FF0000", midpoint = 50) +
  labs(x = "", y = "", title = "Observed Gini index for standard output") +
  theme_bw() + 
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),strip.text = element_text(size = 16),
        axis.title = element_text(size = 14))
p <- ggpubr::ggarrange(p1,p2,ncol = 2,align = "v")
p <- ggpubr::annotate_figure(p,
                             top = text_grob(paste0("Results for 2010 with K = ",K_opt," clusters"),face = "bold",size = 24))
p
ggexport(p,width = 1800, height = 1200, res = 150, filename = paste0("ResClusterAgro_K",K_opt,"_2010.png"))





################################
##### SCSAR model for 2020 #####
################################
dta2020 <- Data %>% 
  dplyr::filter(time == 2020,
                farmtype_lab %in% c("Total"), uaarea_lab %in% c("Total"),
                organic_lab %in% c("Total"),so_eur %in% c("KE0"))
dta2020b <- na.omit(dta2020)

### Sp and W matrix
non_empty_geometries <- !st_is_empty(dta2020b$geometry)
dta2020b_non_empty <- dta2020b[non_empty_geometries, ]
sp_object <- as(dta2020b_non_empty, "Spatial")
nb <- poly2nb(sp_object)
listw <- nb2listw(nb, style = "B", zero.policy=TRUE)
Wmat <- listw2mat(listw)
Sp <- coordinates(sp_object)

### Covariates
Xmat_2020 <- cbind(
  ##### Wealth
  # Per capita GDP: overall regional wealth
  dta2020b_non_empty$GDPPC_PPS2020,
  ##### Labor market
  # Share of employment in agriculture: relevance of agricultural industry on the regional labor market
  dta2020b_non_empty$EmpPersThs_NACE_A/dta2020b_non_empty$EmpPersThs_NACE_Total*100,
  # Hours worked per agro-employed: agricultural labor market intensity
  dta2020b_non_empty$WorkedHoursThs_NACE_A/dta2020b_non_empty$EmpPersThs_NACE_A,
  ##### Economy
  # Gross value added per agro-employed: agricultural productivity intensity
  dta2020b_non_empty$GVA_MillEuro_NACE_A/dta2020b_non_empty$EmpPersThs_NACE_A,
  # Investment per agro-employed: propensity to invest according to the economic size
  dta2020b_non_empty$GFCF_MillEuro_NACE_A/dta2020b_non_empty$EmpPersThs_NACE_A,
  # Share of agricultural GVA on total GVA: relevance of agricultural industry on the regional economy
  dta2020b_non_empty$GVA_MillEuro_NACE_A/dta2020b_non_empty$GVA_MillEuro_NACE_Total*100,
  ##### Landscape and environment
  # Share of agricultural land: relevance of agricultural industry on the regional activities
  dta2020b_non_empty$LandUseKM2_Agro/dta2020b_non_empty$LandUseKM2_Total*100,
  # Average altitude: geography and landscape
  dta2020b_non_empty$Alt_mean,
  # Heating degree days (HDD): proxy of temperature and weather conditions
  dta2020b_non_empty$HDD
)

### Estimate the models
reg0 <- Cl_spatialReg(Y=dta2020b_non_empty$Gini_SO, X=Xmat_2020, Sp=Sp, G=K_opt, W=Wmat, type="lagsarlm", Phi = 0.50)
regnull <- spatialreg::lagsarlm(dta2020b_non_empty$Gini_SO ~ Xmat_2020, listw=listw, zero.policy = T)
# Change names to coefficients
for (m in 1:K_opt) {
  names(reg0$Results[[m]]$coefficients) <- c("Intercept", "GDP pc",
                                             "Share of agro employment","Worked hours per agro worker",
                                             "GVA per agro worker","GFCF per agro worker","Share of agro GVA",
                                             "Share of agro land","Average altitude","HDD")
}
names(regnull$coefficients) <- c("Intercept", "GDP pc",
                                 "Share of agro employment","Worked hours per agro worker",
                                 "GVA per agro worker","GFCF per agro worker","Share of agro GVA",
                                 "Share of agro land","Average altitude","HDD")
### Store in latex
reglist <- c(list(regnull),reg0$Results[1:K_opt])
texreg_table <- texreg(l = reglist,
                       file = paste0("ResClusterAgro_K",K_opt,"_2020.tex"),
                       caption = paste0("Estimated pooled and cluster-wise coefficients for 2020 with K = ",K_opt),
                       custom.model.names = c("Pooled",paste0("G=",1:K_opt)))
texreg(l = reglist)
table(reg0$group)
dta2020_clust <- dta2020b_non_empty
dta2020_clust$Cl_hat <- as.factor(reg0$group)
colnames(dta2020_clust)[which(colnames(dta2020_clust) == "Cl_hat")] <- paste0("clust_K", K_opt)

### Mapping
p1 <- dta2020_clust %>% 
  select(time,'Clustering' = .data[[paste0("clust_K", K_opt)]]) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Clustering)) + 
  geom_sf(data = Countries, linewidth = 1.05, show.legend=FALSE, alpha=0,color="#000000") + 
  labs(x = "", y = "",  title = "Estimated cluster structure") +
  theme_bw() + 
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),strip.text = element_text(size = 16),
        axis.title = element_text(size = 14))
p2 <- dta2020_clust %>% 
  select(time,Gini = Gini_SO) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Gini)) + 
  geom_sf(data = Countries, linewidth = 1.05, show.legend=FALSE, alpha=0,color="#000000") + 
  scale_fill_gradient2(mid = "#FFFFFF",low = "#00FF00",high = "#FF0000", midpoint = 50) +
  labs(x = "", y = "", title = "Observed Gini index for standard output") +
  theme_bw() + 
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),strip.text = element_text(size = 16),
        axis.title = element_text(size = 14))
p <- ggpubr::ggarrange(p1,p2,ncol = 2,align = "v")
p <- ggpubr::annotate_figure(p,
                             top = text_grob(paste0("Results for 2020 with K = ",K_opt," clusters"),face = "bold",size = 24))
ggexport(p,width = 1800, height = 1200, res = 150, filename = paste0("ResClusterAgro_K",K_opt,"_2020.png"))



