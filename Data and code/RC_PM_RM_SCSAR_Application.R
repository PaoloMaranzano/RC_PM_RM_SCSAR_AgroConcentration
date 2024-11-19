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
library(viridis)

### Our library
library(SCDA)

rm(list=ls())


###########################
########## Setup ##########
###########################
setwd("C:/Users/paulm/OneDrive/Documenti/GitHub/RC_PM_RM_SCSAR_AgroConcentration/Data and code")
'%notin%' <- Negate('%in%')
cols <- c("k1" = "#F8766D",
          "k2" = "#00BA38",
          "k3" = "#619CFF",
          "k4" = "orange")
Tuning <- FALSE
Correlation <- TRUE
WRowStd <- FALSE





#########################################
########## Auxiliary functions ##########
#########################################
# 
# source("H:/Il mio Drive/SpatialClustering/SpatReg_SCR_MatteraMaranzano/SCSR/R/SCSR_InfoCrit.R")
# source("H:/Il mio Drive/SpatialClustering/SpatReg_SCR_MatteraMaranzano/SCSR/R/SCSR_Estim.R")
# source("H:/Il mio Drive/SpatialClustering/SpatReg_SCR_MatteraMaranzano/SCSR/R/Elbow_finder.R")
# source("H:/Il mio Drive/SpatialClustering/SpatReg_SCR_MatteraMaranzano/SCSR/R/SpatReg_Extract.R")
# source("H:/Il mio Drive/SpatialClustering/SpatReg_SCR_MatteraMaranzano/SCSR/R/SpatReg_GoF.R")
# source("H:/Il mio Drive/SpatialClustering/SpatReg_SCR_MatteraMaranzano/SCSR/R/SpatReg_Perf.R")
# source("H:/Il mio Drive/SpatialClustering/SpatReg_SCR_MatteraMaranzano/SCSR/R/SpatReg_PseudoR2.R")
# source("H:/Il mio Drive/SpatialClustering/SpatReg_SCR_MatteraMaranzano/SCSR/R/sarlogLik_i.R")
# source("H:/Il mio Drive/SpatialClustering/SpatReg_SCR_MatteraMaranzano/SCSR/R/semlogLik_i.R")





##############################################################################################
########## Application: Analysis of the agricultural market concentration in Europe ##########
##############################################################################################

##### Dataset from Eurostat
load("Data_RC_PM_RM_JABES2024.RData")

if (Correlation == TRUE) {
  ##### Linear correlation analysis
  XXnames <- c("GDPPC_PPS2020","Share_AgroEmp",
               "HoursWorked_AgroEmp","GVA_AgroEmp","GFCF_AgroEmp","Share_AgroGVA","Share_AgroLand",
               "Alt_mean", "HDD")
  XX2010 <- Data2010 %>%
    select(Gini_SO, GDPPC_PPS2020, Share_AgroEmp,
           HoursWorked_AgroEmp, GVA_AgroEmp, GFCF_AgroEmp, Share_AgroGVA, Share_AgroLand,
           Alt_mean, HDD)
  st_geometry(XX2010) <- NULL
  colnames(XX2010) <- c("Gini index standard output", "GDP per capita",
                        "Share of agro employment","Worked hours per agro worker",
                        "GVA per agro worker","GFCF per agro worker","Share of agro GVA",
                        "Share of agro land","Average altitude","HDD")
  corr2010 <- XX2010 %>%
    cor(method = "pearson") %>%
    reshape2::melt(na.rm = TRUE) %>%
    rename(Reg1 = Var1, Reg2 = Var2, Pearson = value)
  corrmat2010 <- corr2010 %>%
    pivot_longer(cols = 3:last_col()) %>%
    # Standardize for fill/color scaling in tile
    group_by(name) %>%
    mutate(value_std = scale(x = value, center = T, scale = T)) %>%
    ungroup() %>%
    arrange(desc(Reg1),desc(Reg2))
  p_corr2010 <- corrmat2010 %>%
    ggplot(mapping = aes(x = Reg2, y = Reg1, fill = value_std)) +
    geom_tile(color = "white") +
    scale_fill_viridis(option="magma", alpha = 0.6, name = "Measure", begin = 1, end = 0.60) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
    coord_fixed() +
    geom_text(aes(x = Reg2, y = Reg1, label = round(value,2)), color = "black", size = 3) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "",
      plot.title = element_text(size = 20, face = "bold",hjust = +0.85)
    ) +
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5)) +
    labs(title = "Pearson's linear correlation index for agro-related variables in 2010")
  ggexport(p_corr2010,width = 2400, height = 2000, res = 250, filename = "Correlation2010.png")

  XX2020 <- Data2010 %>%
    select(Gini_SO, GDPPC_PPS2020, Share_AgroEmp,
           HoursWorked_AgroEmp, GVA_AgroEmp, GFCF_AgroEmp, Share_AgroGVA, Share_AgroLand,
           Alt_mean, HDD)
  st_geometry(XX2020) <- NULL
  colnames(XX2020) <- c("Gini index standard output", "GDP per capita",
                        "Share of agro employment","Worked hours per agro worker",
                        "GVA per agro worker","GFCF per agro worker","Share of agro GVA",
                        "Share of agro land","Average altitude","HDD")
  corr2020 <- XX2020 %>%
    cor(method = "pearson") %>%
    reshape2::melt(na.rm = TRUE) %>%
    rename(Reg1 = Var1, Reg2 = Var2, Pearson = value)
  corrmat2020 <- corr2020 %>%
    pivot_longer(cols = 3:last_col()) %>%
    # Standardize for fill/color scaling in tile
    group_by(name) %>%
    mutate(value_std = scale(x = value, center = T, scale = T)) %>%
    ungroup() %>%
    arrange(desc(Reg1),desc(Reg2))
  p_corr2020 <- corrmat2020 %>%
    ggplot(mapping = aes(x = Reg2, y = Reg1, fill = value_std)) +
    geom_tile(color = "white") +
    scale_fill_viridis(option="magma", alpha = 0.6, name = "Measure", begin = 1, end = 0.60) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
    coord_fixed() +
    geom_text(aes(x = Reg2, y = Reg1, label = round(value,2)), color = "black", size = 3) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "",
      plot.title = element_text(size = 20, face = "bold",hjust = +0.85)
    ) +
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5)) +
    labs(title = "Pearson's linear correlation index for agro-related variables in 2020")
  ggexport(p_corr2020,width = 2400, height = 2000, res = 250, filename = "Correlation2020.png")
}





######################################################
##### Parameters tuning with likelihood criteria #####
######################################################

##### Setup
sp_object <- as(Data2020, "Spatial")
nb <- poly2nb(sp_object)
listW <- nb2listw(nb, style = ifelse(WRowStd == TRUE,"W","B"), zero.policy=TRUE)
if (Tuning == TRUE) {
  # 2010
  SCSAR_IC_2010 <- SCSR_InfoCrit(Formula = paste0("Gini_SO ~ ",paste(XXnames,collapse = " + ")),
                                 Data_sf = Data2010,listW = listW, Type="SCSAR",
                                 CenterVars = FALSE, ScaleVars = FALSE,
                                 Maxitr = 100, Phi.set = c(0,0.25,0.50,0.75,1,1.25), G.set=c(1,2,3,4,5,6))
  # 2020
  SCSAR_IC_2020 <- SCSR_InfoCrit(Formula = paste0("Gini_SO ~ ",paste(XXnames,collapse = " + ")),
                                 Data_sf = Data2020,listW = listW, Type="SCSAR",
                                 CenterVars = FALSE, ScaleVars = FALSE,
                                 Maxitr = 100, Phi.set = c(0,0.25,0.50,0.75,1,1.25), G.set=c(1,2,3,4,5,6))

  TuningIC <- data.frame(
    rbind(cbind(Year = 2010,SCSAR_IC_2010$IC),
          cbind(Year = 2020,SCSAR_IC_2020$IC))
  )
  save(TuningIC,file = paste0(ifelse(WRowStd == TRUE,"W_","B_"),"TuningIC.RData"))
  p <- TuningIC %>%
    pivot_longer(cols = c("BIC","AIC"), names_to = "IC", values_to = "Value") %>%
    mutate(Phi = as.factor(Phi)) %>%
    ggplot(mapping = aes(x = as.factor(G), y = Value, group = Phi)) +
    geom_line(mapping = aes(col = Phi), size = 1.5) +
    geom_point(mapping = aes(col = Phi), size = 4) +
    facet_grid(rows = vars(IC), cols = vars(Year)) +
    theme_bw() +
    # scale_y_continuous(breaks = seq(1000,1700,50)) +
    labs(x = "Number of clusters (K)", y = "",
         title = latex2exp::TeX("Tuning of the spatial penalty parameter ($\\phi$) and number of clusters ($\\K$)"),
         subtitle = paste0("Spatial weighting matrix for the autocorrelation term is ",
                           ifelse(WRowStd==TRUE,"row-standardized","binary"),
                           ", while the spatial penalty matrix is binary")) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          plot.title = element_text(face = "bold", size = 24),
          plot.subtitle = element_text(size = 14),
          strip.text = element_text(size = 14))
  p
  ggexport(p,width = 1800, height = 1200, res = 150, filename = paste0(ifelse(WRowStd == TRUE,"W_","B_"),
                                                                       "ParsTuning.png"))
}






################################
##### SCSAR model for 2010 #####
################################

CombFinal <- expand.grid(
  G = c(3,4),
  Phi = c(0.25,0.50,1)
)


for (i in 1:dim(CombFinal)[1]) {
  ### Model estimation
  G_opt <- CombFinal$G[i]
  Phi_opt <- CombFinal$Phi[i]

  SCSAR_2010 <- SCSR_Estim(Formula = paste0("Gini_SO ~ ",paste(XXnames,collapse = " + ")),
                           Data_sf = Data2010, listW = listW, G=G_opt, Type="SCSAR", Phi = Phi_opt,
                           CenterVars = FALSE, ScaleVars = FALSE)
  pooledSAR_2010 <- SCSR_Estim(Formula = paste0("Gini_SO ~ ",paste(XXnames,collapse = " + ")),
                               Data_sf = Data2010, listW = listW,G=1, Type="SCSAR", Phi = Phi_opt,
                               CenterVars = FALSE, ScaleVars = FALSE)
  emptySAR_2010 <- SCSR_Estim(Formula = "Gini_SO ~ 1",
                              Data_sf = Data2010, listW = listW,G=1, Type="SCSAR", Phi = Phi_opt,
                              CenterVars = FALSE, ScaleVars = FALSE)
  SCLM_2010 <- SCSR_Estim(Formula = paste0("Gini_SO ~ ",paste(XXnames,collapse = " + ")),
                          Data_sf = Data2010, listW = listW, G=G_opt, Type="SCLM", Phi = Phi_opt,
                          CenterVars = FALSE, ScaleVars = FALSE)
  pooledLM_2010 <- SCSR_Estim(Formula = paste0("Gini_SO ~ ",paste(XXnames,collapse = " + ")),
                              Data_sf = Data2010, listW = listW,G=1, Type="SCLM", Phi = Phi_opt,
                              CenterVars = FALSE, ScaleVars = FALSE)
  emptyLM_2010 <- SCSR_Estim(Formula = "Gini_SO ~ 1",
                             Data_sf = Data2010, listW = listW,G=1, Type="SCLM", Phi = Phi_opt,
                             CenterVars = FALSE, ScaleVars = FALSE)

  ##### Change names to coefficients
  for (m in 1:G_opt) {
    names(SCSAR_2010$ClusterFitModels[[m]]$coefficients) <-
      names(SCLM_2010$ClusterFitModels[[m]]$coefficients) <- c("Intercept", "GDP pc",
                                                               "Share of agro employment","Worked hours per agro worker",
                                                               "GVA per agro worker","GFCF per agro worker","Share of agro GVA",
                                                               "Share of agro land","Average altitude","HDD")
  }
  names(pooledSAR_2010$ClusterFitModels[[1]]$coefficients) <-
    names(pooledLM_2010$ClusterFitModels[[1]]$coefficients) <- c("Intercept", "GDP pc",
                                                                 "Share of agro employment","Worked hours per agro worker",
                                                                 "GVA per agro worker","GFCF per agro worker","Share of agro GVA",
                                                                 "Share of agro land","Average altitude","HDD")
  names(emptySAR_2010$ClusterFitModels[[1]]$coefficients) <-
    names(emptyLM_2010$ClusterFitModels[[1]]$coefficients) <- c("Intercept")

  ### Add clustering results to data
  Data2010_clust <- Data2010
  Data2010_clust$Cl_hat <- as.factor(SCSAR_2010$Group)
  colnames(Data2010_clust)[which(colnames(Data2010_clust) == "Cl_hat")] <- paste0("clust_SCSAR_G", G_opt)
  Data2010_clust$Cl_hat <- as.factor(SCLM_2010$Group)
  colnames(Data2010_clust)[which(colnames(Data2010_clust) == "Cl_hat")] <- paste0("clust_SCLM_G", G_opt)

  ### Mapping
  # WRowStd == TRUE and G = 3 <- 2-1-3
  # WRowStd == TRUE and G = 3 <- 2-3-1-4
  p1 <- Data2010_clust %>%
    select(Year,'Clustering' = .data[[paste0("clust_SCSAR_G", G_opt)]]) %>% {
      if (G_opt == 3) {
        mutate(., Clustering = case_when(Clustering == 1 ~ "k2",
                                      Clustering == 2 ~ "k1",
                                      Clustering == 3 ~ "k3"))
      } else if (G_opt == 4) {
        mutate(., Clustering = case_when(.data$Clustering == 1 ~ "k2",
                                      .data$Clustering == 2 ~ "k3",
                                      .data$Clustering == 3 ~ "k1",
                                      .data$Clustering == 4 ~ "k4"))
      }
    } %>%
    ggplot() +
    geom_sf(mapping = aes(fill = Clustering)) +
    geom_sf(data = Countries, linewidth = 1.05, show.legend=FALSE, alpha=0,color="#000000") +
    labs(x = "", y = "",  title = "Estimated cluster structure") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.text = element_text(size = 12),strip.text = element_text(size = 16),
          axis.title = element_text(size = 14)) +
    scale_fill_manual("Clusters",values = cols)
  p2 <- Data2010_clust %>%
    select(Year,Gini = Gini_SO) %>%
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
                               top = text_grob(latex2exp::TeX(paste0("SCSAR results for 2010 with K = ",G_opt," clusters and $\\phi$ = ",Phi_opt)),
                                               face = "bold",size = 24,vjust = +1, hjust = +0.50),
                               bottom = text_grob(paste0("Note: the spatial weighting matrix for the autocorrelation term is ",
                                                         ifelse(WRowStd==TRUE,"row-standardized","binary"),
                                                         ", while the spatial penalty matrix is binary"),size=12,
                                                  vjust = -1,hjust=+0.50))
  p
  ggexport(p,width = 1800, height = 1200, res = 150, filename = paste0(ifelse(WRowStd == TRUE,"W_","B_"),
                                                                       "SCSAR_ClusterAgro_G",G_opt,"_phi",tm::removePunctuation(as.character(Phi_opt)),"_2010.png"))


  p3 <- Data2010_clust %>%
    select(Year,'Clustering' = .data[[paste0("clust_SCLM_G", G_opt)]]) %>% {
      if (G_opt == 3) {
        mutate(., Clustering = case_when(Clustering == 1 ~ "k2",
                                         Clustering == 2 ~ "k1",
                                         Clustering == 3 ~ "k3"))
      } else if (G_opt == 4) {
        mutate(., Clustering = case_when(.data$Clustering == 1 ~ "k2",
                                         .data$Clustering == 2 ~ "k3",
                                         .data$Clustering == 3 ~ "k1",
                                         .data$Clustering == 4 ~ "k4"))
      }
    } %>%
    ggplot() +
    geom_sf(mapping = aes(fill = Clustering)) +
    geom_sf(data = Countries, linewidth = 1.05, show.legend=FALSE, alpha=0,color="#000000") +
    labs(x = "", y = "",  title = "Estimated cluster structure") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),strip.text = element_text(size = 16),
          axis.title = element_text(size = 14)) +
    scale_fill_manual("Clusters",values = cols)
  p4 <- Data2010_clust %>%
    select(Year,Gini = Gini_SO) %>%
    ggplot() +
    geom_sf(mapping = aes(fill = Gini)) +
    geom_sf(data = Countries, linewidth = 1.05, show.legend=FALSE, alpha=0,color="#000000") +
    scale_fill_gradient2(mid = "#FFFFFF",low = "#00FF00",high = "#FF0000", midpoint = 50) +
    labs(x = "", y = "", title = "Observed Gini index for standard output") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),strip.text = element_text(size = 16),
          axis.title = element_text(size = 14))
  p <- ggpubr::ggarrange(p3,p4,ncol = 2,align = "v")
  p <- ggpubr::annotate_figure(p,
                               top = text_grob(latex2exp::TeX(paste0("SCLM results for 2010 with K = ",G_opt," clusters and $\\phi$ = ",Phi_opt)),
                                               face = "bold",size = 24,vjust = +1,hjust=+0.50),
                               bottom = text_grob(paste0("Note: the spatial weighting matrix for the autocorrelation term is ",
                                                         ifelse(WRowStd==TRUE,"row-standardized","binary"),
                                                         ", while the spatial penalty matrix is binary"),size=12,
                                                  vjust = -1,hjust=+0.50))
  p
  ggexport(p,width = 1800, height = 1200, res = 150, filename = paste0(ifelse(WRowStd == TRUE,"W_","B_"),
                                                                       "SCLM_ClusterAgro_G",G_opt,"_phi",tm::removePunctuation(as.character(Phi_opt)),"_2010.png"))



}





##########################################
##### SCSAR and SCLM models for 2020 #####
##########################################

for (i in 1:dim(CombFinal)[1]) {
  ### Model estimation
  G_opt <- CombFinal$G[i]
  Phi_opt <- CombFinal$Phi[i]

  SCSAR_2020 <- SCSR_Estim(Formula = paste0("Gini_SO ~ ",paste(XXnames,collapse = " + ")),
                           Data_sf = Data2020, listW = listW, G=G_opt, Type="SCSAR", Phi = Phi_opt,
                           CenterVars = FALSE, ScaleVars = FALSE)
  pooledSAR_2020 <- SCSR_Estim(Formula = paste0("Gini_SO ~ ",paste(XXnames,collapse = " + ")),
                               Data_sf = Data2020, listW = listW,G=1, Type="SCSAR", Phi = Phi_opt,
                               CenterVars = FALSE, ScaleVars = FALSE)
  emptySAR_2020 <- SCSR_Estim(Formula = "Gini_SO ~ 1",
                              Data_sf = Data2020, listW = listW,G=1, Type="SCSAR", Phi = Phi_opt,
                              CenterVars = FALSE, ScaleVars = FALSE)
  SCLM_2020 <- SCSR_Estim(Formula = paste0("Gini_SO ~ ",paste(XXnames,collapse = " + ")),
                          Data_sf = Data2020, listW = listW, G=G_opt, Type="SCLM", Phi = Phi_opt,
                          CenterVars = FALSE, ScaleVars = FALSE)
  pooledLM_2020 <- SCSR_Estim(Formula = paste0("Gini_SO ~ ",paste(XXnames,collapse = " + ")),
                              Data_sf = Data2020, listW = listW,G=1, Type="SCLM", Phi = Phi_opt,
                              CenterVars = FALSE, ScaleVars = FALSE)
  emptyLM_2020 <- SCSR_Estim(Formula = "Gini_SO ~ 1",
                             Data_sf = Data2020, listW = listW,G=1, Type="SCSAR", Phi = Phi_opt,
                             CenterVars = FALSE, ScaleVars = FALSE)

  ##### Change names to coefficients
  for (m in 1:G_opt) {
    names(SCSAR_2020$ClusterFitModels[[m]]$coefficients) <-
      names(SCLM_2020$ClusterFitModels[[m]]$coefficients) <- c("Intercept", "GDP pc",
                                                               "Share of agro employment","Worked hours per agro worker",
                                                               "GVA per agro worker","GFCF per agro worker","Share of agro GVA",
                                                               "Share of agro land","Average altitude","HDD")
  }
  names(pooledSAR_2020$ClusterFitModels[[1]]$coefficients) <-
    names(pooledLM_2020$ClusterFitModels[[1]]$coefficients) <- c("Intercept", "GDP pc",
                                                                 "Share of agro employment","Worked hours per agro worker",
                                                                 "GVA per agro worker","GFCF per agro worker","Share of agro GVA",
                                                                 "Share of agro land","Average altitude","HDD")
  names(emptySAR_2020$ClusterFitModels[[1]]$coefficients) <-
    names(emptyLM_2020$ClusterFitModels[[1]]$coefficients) <- c("Intercept")

  ### Add clustering results to data
  Data2020_clust <- Data2020
  Data2020_clust$Cl_hat <- as.factor(SCSAR_2020$Group)
  colnames(Data2020_clust)[which(colnames(Data2020_clust) == "Cl_hat")] <- paste0("clust_SCSAR_G", G_opt)
  Data2020_clust$Cl_hat <- as.factor(SCLM_2020$Group)
  colnames(Data2020_clust)[which(colnames(Data2020_clust) == "Cl_hat")] <- paste0("clust_SCLM_G", G_opt)

  ### Mapping
  # WRowStd == TRUE and G = 3 <- 2-1-3
  # WRowStd == TRUE and G = 3 <- 2-3-1-4
  p1 <- Data2020_clust %>%
    select(Year,'Clustering' = .data[[paste0("clust_SCSAR_G", G_opt)]]) %>% {
      if (G_opt == 3) {
        mutate(., Clustering = case_when(Clustering == 1 ~ "k2",
                                         Clustering == 2 ~ "k1",
                                         Clustering == 3 ~ "k3"))
      } else if (G_opt == 4) {
        mutate(., Clustering = case_when(.data$Clustering == 1 ~ "k2",
                                         .data$Clustering == 2 ~ "k3",
                                         .data$Clustering == 3 ~ "k1",
                                         .data$Clustering == 4 ~ "k4"))
      }
    } %>%
    ggplot() +
    geom_sf(mapping = aes(fill = Clustering)) +
    geom_sf(data = Countries, linewidth = 1.05, show.legend=FALSE, alpha=0,color="#000000") +
    labs(x = "", y = "",  title = "Estimated cluster structure") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.text = element_text(size = 12),strip.text = element_text(size = 16),
          axis.title = element_text(size = 14)) +
    scale_fill_manual("Clusters",values = cols)
  p2 <- Data2020_clust %>%
    select(Year,Gini = Gini_SO) %>%
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
                               top = text_grob(latex2exp::TeX(paste0("SCSAR results for 2020 with K = ",G_opt," clusters and $\\phi$ = ",Phi_opt)),
                                               face = "bold",size = 24,vjust = +1, hjust = +0.50),
                               bottom = text_grob(paste0("Note: the spatial weighting matrix for the autocorrelation term is ",
                                                         ifelse(WRowStd==TRUE,"row-standardized","binary"),
                                                         ", while the spatial penalty matrix is binary"),size=12,
                                                  vjust = -1,hjust=+0.50))
  p
  ggexport(p,width = 1800, height = 1200, res = 150, filename = paste0(ifelse(WRowStd == TRUE,"W_","B_"),
                                                                       "SCSAR_ClusterAgro_G",G_opt,"_phi",tm::removePunctuation(as.character(Phi_opt)),"_2020.png"))


  p3 <- Data2020_clust %>%
    select(Year,'Clustering' = .data[[paste0("clust_SCLM_G", G_opt)]]) %>% {
      if (G_opt == 3) {
        mutate(., Clustering = case_when(Clustering == 1 ~ "k2",
                                         Clustering == 2 ~ "k1",
                                         Clustering == 3 ~ "k3"))
      } else if (G_opt == 4) {
        mutate(., Clustering = case_when(.data$Clustering == 1 ~ "k2",
                                         .data$Clustering == 2 ~ "k3",
                                         .data$Clustering == 3 ~ "k1",
                                         .data$Clustering == 4 ~ "k4"))
      }
    } %>%
    ggplot() +
    geom_sf(mapping = aes(fill = Clustering)) +
    geom_sf(data = Countries, linewidth = 1.05, show.legend=FALSE, alpha=0,color="#000000") +
    labs(x = "", y = "",  title = "Estimated cluster structure") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),strip.text = element_text(size = 16),
          axis.title = element_text(size = 14)) +
    scale_fill_manual("Clusters",values = cols)
  p4 <- Data2020_clust %>%
    select(Year,Gini = Gini_SO) %>%
    ggplot() +
    geom_sf(mapping = aes(fill = Gini)) +
    geom_sf(data = Countries, linewidth = 1.05, show.legend=FALSE, alpha=0,color="#000000") +
    scale_fill_gradient2(mid = "#FFFFFF",low = "#00FF00",high = "#FF0000", midpoint = 50) +
    labs(x = "", y = "", title = "Observed Gini index for standard output") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12),strip.text = element_text(size = 16),
          axis.title = element_text(size = 14))
  p <- ggpubr::ggarrange(p3,p4,ncol = 2,align = "v")
  p <- ggpubr::annotate_figure(p,
                               top = text_grob(latex2exp::TeX(paste0("SCLM results for 2020 with K = ",G_opt," clusters and $\\phi$ = ",Phi_opt)),
                                               face = "bold",size = 24,vjust = +1,hjust=+0.50),
                               bottom = text_grob(paste0("Note: the spatial weighting matrix for the autocorrelation term is ",
                                                         ifelse(WRowStd==TRUE,"row-standardized","binary"),
                                                         ", while the spatial penalty matrix is binary"),size=12,
                                                  vjust = -1,hjust=+0.50))
  p
  ggexport(p,width = 1800, height = 1200, res = 150, filename = paste0(ifelse(WRowStd == TRUE,"W_","B_"),
                                                                       "SCLM_ClusterAgro_G",G_opt,"_phi",tm::removePunctuation(as.character(Phi_opt)),"_2020.png"))

}





##############################
##### Store Latex tables #####
##############################

for (i in 1:dim(CombFinal)[1]) {
  G_opt <- CombFinal$G[i]
  Phi_opt <- CombFinal$Phi[i]

  ### SCSAR table
  if (G_opt == 3) {
    reglist <- c(
      emptySAR_2010$ClusterFitModels[1],pooledSAR_2010$ClusterFitModels[1],
      SCSAR_2010$ClusterFitModels[2],
      SCSAR_2010$ClusterFitModels[1],
      SCSAR_2010$ClusterFitModels[3],
      emptySAR_2020$ClusterFitModels[1],pooledSAR_2020$ClusterFitModels[1],
      SCSAR_2020$ClusterFitModels[2],
      SCSAR_2020$ClusterFitModels[1],
      SCSAR_2020$ClusterFitModels[3]
    )
    Wlist <- c(
      emptySAR_2010$listW_g,pooledSAR_2010$listW_g,
      SCSAR_2010$listW_g[2],
      SCSAR_2010$listW_g[1],
      SCSAR_2010$listW_g[3],
      emptySAR_2020$listW_g,pooledSAR_2020$listW_g,
      SCSAR_2020$listW_g[2],
      SCSAR_2020$listW_g[1],
      SCSAR_2020$listW_g[3]
    )
  } else if (G_opt == 4) {
    reglist <- c(
      emptySAR_2010$ClusterFitModels[1],pooledSAR_2010$ClusterFitModels[1],
      SCSAR_2010$ClusterFitModels[2],
      SCSAR_2010$ClusterFitModels[3],
      SCSAR_2010$ClusterFitModels[1],
      SCSAR_2010$ClusterFitModels[4],
      emptySAR_2020$ClusterFitModels[1],pooledSAR_2020$ClusterFitModels[1],
      SCSAR_2020$ClusterFitModels[2],
      SCSAR_2020$ClusterFitModels[3],
      SCSAR_2020$ClusterFitModels[1],
      SCSAR_2020$ClusterFitModels[4]
    )
    Wlist <- c(
      emptySAR_2010$listW_g,pooledSAR_2010$listW_g,
      SCSAR_2010$listW_g[2],
      SCSAR_2010$listW_g[3],
      SCSAR_2010$listW_g[1],
      SCSAR_2010$listW_g[4],
      emptySAR_2020$listW_g,pooledSAR_2020$listW_g,
      SCSAR_2020$listW_g[2],
      SCSAR_2020$listW_g[3],
      SCSAR_2020$listW_g[1],
      SCSAR_2020$listW_g[4]
    )
  }
  Mods_gof <- SpatReg_GoF(SRModel_list = reglist, SRModel_W_list = Wlist)$SRModel_GoF
  Mods_gof <- round(Mods_gof,3)
  texreg_table <- texreg(l = reglist,
                         file = paste0(ifelse(WRowStd == TRUE,"W_","B_"),"ResClusterAgro_SCSAR_G",G_opt,"_phi",tm::removePunctuation(as.character(Phi_opt)),".tex"),
                         caption = paste0("Estimated empty (no covariates), pooled and cluster-wise regression coefficients of SCSAR model for 2010 (columns 2 to 6) and 2020 (columns 7 to 11) using $G=$",
                                          G_opt," clusters, spatial penalty $phi=$",Phi_opt, ", and ",ifelse(WRowStd==TRUE,"row-standardized","binary")," weighting matrix for the autoregressive term."),
                         label = paste0("Tab:SCSAR_Estimates_G",G_opt,"_Phi025"),
                         caption.above = TRUE,
                         custom.gof.rows = list(
                           "Resid.SD" = Mods_gof[1,],
                           "Num. \\ obs" = Mods_gof[2,],
                           "Num. \\ pars." = Mods_gof[3,],
                           "Log Likelihood" = Mods_gof[4,],
                           "AIC" = Mods_gof[5,],
                           "BIC" = Mods_gof[6,],
                           "Pseudo R$^2$" = Mods_gof[7,],
                           "LR test: statistic" = Mods_gof[8,],
                           "LR test: p-value" = Mods_gof[9,],
                           "W test: statistic" = Mods_gof[10,],
                           "W test: p-value" = Mods_gof[11,],
                           "Moran's I on Y: statistic" = Mods_gof[12,],
                           "Moran's I on Y: p-value" = Mods_gof[13,],
                           "Moran's I on residuals: statistic" = Mods_gof[14,],
                           "Moran's I on residuals: p-value" = Mods_gof[15,]
                         ),
                         custom.model.names = c("Empty","Pooled",paste0("G=",1:G_opt),
                                                "Empty","Pooled",paste0("G=",1:G_opt))
  )
  texreg(l = reglist)

  ### LM table
  if (G_opt == 3) {
    reglist_LM <- c(
      emptyLM_2010$ClusterFitModels[1],pooledLM_2010$ClusterFitModels[1],
      SCLM_2010$ClusterFitModels[2],
      SCLM_2010$ClusterFitModels[1],
      SCLM_2010$ClusterFitModels[3],
      emptyLM_2020$ClusterFitModels[1],pooledLM_2020$ClusterFitModels[1],
      SCLM_2020$ClusterFitModels[2],
      SCLM_2020$ClusterFitModels[1],
      SCLM_2020$ClusterFitModels[3]
    )
    Wlist_LM <- c(
      emptyLM_2010$listW_g,pooledLM_2010$listW_g,
      SCLM_2010$listW_g[2],
      SCLM_2010$listW_g[1],
      SCLM_2010$listW_g[3],
      emptyLM_2020$listW_g,pooledLM_2020$listW_g,
      SCLM_2020$listW_g[2],
      SCLM_2020$listW_g[1],
      SCLM_2020$listW_g[3]
    )
  } else if (G_opt == 4) {
    reglist_LM <- c(
      emptyLM_2010$ClusterFitModels[1],pooledLM_2010$ClusterFitModels[1],
      SCLM_2010$ClusterFitModels[2],
      SCLM_2010$ClusterFitModels[3],
      SCLM_2010$ClusterFitModels[1],
      SCLM_2010$ClusterFitModels[4],
      emptyLM_2020$ClusterFitModels[1],pooledLM_2020$ClusterFitModels[1],
      SCLM_2020$ClusterFitModels[2],
      SCLM_2020$ClusterFitModels[3],
      SCLM_2020$ClusterFitModels[1],
      SCLM_2020$ClusterFitModels[4]
    )
    Wlist_LM <- c(
      emptyLM_2010$listW_g,pooledLM_2010$listW_g,
      SCLM_2010$listW_g[2],
      SCLM_2010$listW_g[3],
      SCLM_2010$listW_g[1],
      SCLM_2010$listW_g[4],
      emptyLM_2020$listW_g,pooledLM_2020$listW_g,
      SCLM_2020$listW_g[2],
      SCLM_2020$listW_g[3],
      SCLM_2020$listW_g[1],
      SCLM_2020$listW_g[4]
    )
  }
  Mods_gof_LM <- SpatReg_GoF(SRModel_list = reglist_LM, SRModel_W_list = Wlist_LM)$SRModel_GoF
  Mods_gof_LM <- round(Mods_gof_LM,3)
  texreg_table <- texreg(l = reglist_LM,
                         file = paste0(ifelse(WRowStd == TRUE,"W_","B_"),"ResClusterAgro_SCLM_G",G_opt,"_phi",tm::removePunctuation(as.character(Phi_opt)),".tex"),
                         caption = paste0("Estimated empty (no covariates), pooled and cluster-wise regression coefficients of SCLM model for 2010 (columns 2 to 6) and 2020 (columns 7 to 11) using $G=$",
                                          G_opt," clusters, spatial penalty $phi=$",Phi_opt, ", and ",ifelse(WRowStd==TRUE,"row-standardized","binary")," weighting matrix for the autoregressive term."),
                         label = paste0("Tab:SCLM_Estimates_G",G_opt,"_Phi025"),
                         caption.above = TRUE,
                         custom.gof.rows = list(
                           "Resid.SD" = Mods_gof_LM[1,],
                           "Num. \\ obs" = Mods_gof_LM[2,],
                           "Num. \\ pars." = Mods_gof_LM[3,],
                           "Log Likelihood" = Mods_gof_LM[4,],
                           "AIC" = Mods_gof_LM[5,],
                           "BIC" = Mods_gof_LM[6,],
                           "Pseudo R$^2$" = Mods_gof_LM[7,],
                           "LR test: statistic" = Mods_gof_LM[8,],
                           "LR test: p-value" = Mods_gof_LM[9,],
                           "W test: statistic" = Mods_gof_LM[10,],
                           "W test: p-value" = Mods_gof_LM[11,],
                           "Moran's I on Y: statistic" = Mods_gof_LM[12,],
                           "Moran's I on Y: p-value" = Mods_gof_LM[13,],
                           "Moran's I on residuals: statistic" = Mods_gof_LM[14,],
                           "Moran's I on residuals: p-value" = Mods_gof_LM[15,]
                         ),
                         custom.model.names = c("Empty","Pooled",paste0("G=",1:G_opt),
                                                "Empty","Pooled",paste0("G=",1:G_opt))
  )


}
