#########################################################################################################################################
##########                              Cerqueti, R., Maranzano, P. & Mattera, R. (2024+)                                      ##########
########## "Spatially-clustered spatial autoregressive models with application to agricultural market concentration in Europe" ##########
##########                                  Auxiliary functions (Sections 2 and 3)                                             ##########
#########################################################################################################################################

########## Cl_spatialReg
##### Aims: Estimate the SC-spatial autoregressive models
##### Arguments:
# Y = n x 1 vector of response variable
# X = n x K matrix of explanatory/covariate variables
# W = n x n matrix with spatial weights (list2mat object)
# Wstand = binary value stating of the W matrix has to be row-standardized (TRUE) or not (FALSE, default)
# Sp = sp object with spatial weights and coordinates
# G = number of clusters
# Phi = Spatial penalty parameter
# type = model specification to be estimated
#      = "lm" for linear model without spatial effects
#      = "lagsarlm" for spatial autoregressive (SAR) model
#      = "errorsarlm" for linear model with spatial autoregressive error term (SEM, Durbin)
#      = "lmSLX" for linear model with spatially-lagged response variable and covariates (SLX)
# Verbose = ...
# zero.policy = TRUE
SCSR_Estim <- function(Y, X, W, Sp, Wstand=FALSE, G=2, Phi=1, maxitr=100, type=c("lm","lagsarlm","errorsarlm","lmSLX"),
                          zero.policy = TRUE, Verbose = TRUE, seed = 123456789){

  ##### Note generali da risolvere:
  # 1. "Wstand" influisce solo nella definizione dei df del lmSLX
  # 2. Errore in
  # lmSLX(Y[Ind==g] ~ 1 , listw=listWg, zero.policy=TRUE)
  # Errore in lmSLX(Y[Ind == g] ~ 1, listw = listWg, zero.policy = TRUE) : intercept-only model, Durbin invalid
  # 3. Controllo sul numero minimo di unità in ogni gruppo per evitare n<p o n=p

  ########## Setup
  if (Verbose == TRUE) {
    TypeName <- switch(EXPR = type,
                       "lm" = "linear model without spatial effects (lm)",
                       "lagsarlm" = "Spatial autoregressive model (SAR)",
                       "errorsarlm" = "linear model with spatial autoregressive error term (SEM)",
                       "lmSLX" = "linear model with spatially-lagged response variable and covariates (SLX)"
    )
    cat(paste0("Spatially-clustered spatial autoregressive model started at ",round(Sys.time()),"\n"))
    cat(paste0("Selected model: ",TypeName," with g=",G," groups and spatial penalty parameter phi=",Phi,"\n"))
  }

  ### convergence criterion
  ep <- 10^(-5)

  ### Seed for replication
  set.seed(seed)

  ### Check matrix format
  if (!is.null(X)){
    if (is.null(colnames(X))) {
      Xnames <- paste0("X",1:ncol(X))
    } else {
      Xnames <- colnames(X)
    }
    X <- as.matrix(X)
  } else {
    Xnames <- c("Intercept")
    X <- as.matrix(rep(1,length(Y)))
  }

  ### Augment data with constant column for predictions
  if (ncol(X) == 1) {
    XX <- as.matrix(rep(1,length(Y)))
    Xnames <- c("Intercept")
  } else {
    XX <- cbind(Intercept = rep(1,length(Y)),X)
    Xnames <- c("Intercept",Xnames)
  }
  ### W as a sparse matrix
  W <- as(W, "sparseMatrix")
  ########## Attenzione!!
  ### Binary W matrix (for penalty term)
  if (TRUE) {
    Wbin <- (W != 0) + 0
    Wbin <- as(Wbin, "sparseMatrix")
  } else {
    Wbin <- W
  }



  ########## Dimensionality
  ### n: Number of observations
  n <- length(Y)
  ### p: Number of regression coefficients (excluding residual variance)
  if (type=="lm"){
    # covariates (including intercept)
    p <- ncol(XX)
    XXnames <- c(Xnames)
  }
  if (type=="lagsarlm"){
    # covariates (including intercept) + rho
    p <- ncol(XX) + 1
    XXnames <- c("rho",Xnames)
  }
  if (type=="errorsarlm"){
    # covariates (including intercept) + lambda
    p <- ncol(XX) + 1
    XXnames <- c("lambda",Xnames)
  }
  if (type=="lmSLX"){
    if (Wstand==FALSE){
      # covariates (including intercept) + lagged covariates + lambda + rho
      # 2*ncol(X) + 2
      p <- 2*(dim(X)[2]+1)
    } else {
      p <- (2*dim(X))+1
    }
  }



  ########## Clustering algorithm
  ## Initial values
  # the number of initial values of k-means
  M <- 20
  WSS <- c()
  CL <- list()
  for(k in 1:M){
    CL[[k]] <- kmeans(Sp, G)
    WSS[k] <- CL[[k]]$tot.withinss
  }
  Ind <- CL[[which.min(WSS)]]$cluster
  Pen <- rep(0, G)
  Beta <- matrix(0, p, G)
  colnames(Beta) <- paste0("G=",1:G)
  rownames(Beta) <- XXnames
  Sig <- rep(1, G)

  ### Iterative algorithm
  val <- 0
  mval <- 0
  fit <- list()

  for(k in 1:maxitr){

    cat("*")
    cval <- val

    ### Penalty term
    Ind.mat <- matrix(0, n, G)
    for(g in 1:G){
      Ind.mat[Ind==g, g] <- 1
    }
    Ind.mat <- as(Ind.mat, "sparseMatrix")
    # penalty term
    # Pen <- W%*%Ind.mat
    Pen <- Wbin%*%Ind.mat

    ### Model parameters (clustered case)
    for(g in 1:G){
      if(length(Ind[Ind==g]) > p+1){
        if(type=="lm"){
          if (ncol(X) == 1) {
            fit[[g]] <- lm(Y[Ind==g] ~ X[Ind==g,] - 1)
          } else {
            fit[[g]] <- lm(Y[Ind==g] ~ X[Ind==g,])
          }
          Beta[,g] <- as.vector( coef(fit[[g]]) )
          resid <- Y - as.vector(XX%*%Beta[,g])
          Sig[g] <- sqrt(mean(resid[Ind==g]^2))
          Sig[g] <- max(Sig[g], 0.1)
        }
        if (type=="lagsarlm"){
          XXX <- cbind(as.vector(W%*%Y),XX)
          listWg <- mat2listw(W[Ind==g,Ind==g])
          if (ncol(X) == 1) {
            fit[[g]] <- lagsarlm(Y[Ind==g] ~ X[Ind==g,] - 1, listw=listWg, zero.policy=zero.policy)
          } else {
            fit[[g]] <- lagsarlm(Y[Ind==g] ~ X[Ind==g,], listw=listWg, zero.policy=zero.policy)
          }
          Beta[,g] <- as.vector( coef(fit[[g]]) )
          Sig[g] <- fit[[g]]$s2
          Sig[g] <- max(Sig[g], 0.1)
        }
        if (type=="errorsarlm"){
          listWg <- mat2listw(W[Ind==g,Ind==g])
          if (ncol(X) == 1) {
            fit[[g]] <- errorsarlm(Y[Ind==g] ~ X[Ind==g,] - 1, listw=listWg, zero.policy=zero.policy)
          } else {
            fit[[g]] <- errorsarlm(Y[Ind==g] ~ X[Ind==g,], listw=listWg, zero.policy=zero.policy)
          }
          Beta[,g] <- as.vector( coef(fit[[g]]) )
          Sig[g] <- fit[[g]]$s2
          Sig[g] <- max(Sig[g], 0.1)
        }
        if (type=="lmSLX"){
          if (Wstand==FALSE){
            WI <- lag.listw(mat2listw(W), XX[,1], zero.policy = zero.policy)
            WX <- spatialreg:::create_WX(XX[,-1],mat2listw(W), zero.policy=zero.policy)
            XXX <- cbind(XX,WI,WX)
          } else {
            WX <- spatialreg:::create_WX(XX[,-1],mat2listw(W), zero.policy=zero.policy)
            XXX <- cbind(XX,WX)
          }
          listWg <- mat2listw(W[Ind==g,Ind==g])
          if (ncol(X) == 1) {
            fit[[g]] <- lmSLX(Y[Ind==g] ~ as.matrix(X[Ind==g,]) - 1, listw=listWg, zero.policy=zero.policy)
          } else {
            fit[[g]] <- lmSLX(Y[Ind==g] ~ as.matrix(X[Ind==g,]), listw=listWg, zero.policy=zero.policy)
          }
          Beta[,g] <- as.vector( coef(fit[[g]]) )
          resid <-  Y-as.vector(XXX%*%Beta[,g])
          Sig[g] <- sqrt(mean(resid[Ind==g]^2))
          Sig[g] <- max(Sig[g], 0.1)
        }
      } # else print error/warning "N_k < P, estimation not possible!"
    }

    ## Grouping (clustered case)
    if(type=="lm"){
      Mu <- XX%*%Beta      # (n,G)-matrix
      ESig <- t(matrix(rep(Sig,n), G, n))    # (n,G)-matrix
      Q <- dnorm(Y, Mu, ESig, log=T) + Phi*Pen     # penalized likelihood
    }
    if (type=="lagsarlm"){
      LL <- matrix(NA, n, G)
      for (g in 1:G) {
        LL[,g] <- as.vector(sarlogLik_i(Y,XX,W,Beta[-1,g],Beta[1,g],Sig[g]))
      }
      Q <- LL + Phi*Pen   # penalized likelihood
    }
    if (type=="errorsarlm"){
      LL <- matrix(NA, n, G)
      for (g in 1:G) {
        LL[,g] <- as.vector(semlogLik_i(Y,XX,W,Beta[-1,g],Beta[1,g],Sig[g]))
      }
      Q <- LL + Phi*Pen   # penalized likelihood
    }
    if (type=="lmSLX"){ ## come lm ma con X lagged in the space
      Mu <- XXX%*%Beta     # (n,G)-matrix
      ESig <- t(matrix(rep(Sig,n), G, n))    # (n,G)-matrix
      Q <- dnorm(Y, Mu, ESig, log=T) + Phi*Pen
    }

    Ind <- apply(Q, 1, which.max)

    ### Value of objective function
    # val = Sum of the G maximum log-likelihood
    if (type == "lm") {
      val <- sum(unlist(lapply(fit,function(x){as.numeric(logLik(fit[[1]]))})))
    } else {
      val <- sum(unlist(lapply(fit,function(x){as.numeric(x$LL)})))
    }

    ### Check the exit conditions
    dd <- abs(cval-val)/abs(val)
    mval <- max(mval, cval)
    if( dd<ep | abs(mval-val)<ep ){ break }
  }
  cat("\n")

  ##### Varying parameters
  if (ncol(X) == 1 & type == "lm") {
    # Location-wise regression coefficients
    sBeta <- as.matrix(Beta[,Ind])
    colnames(sBeta) <- "Intercept"
    # Location-wise error variance
    sSig <- Sig[Ind]
  } else {
    # Location-wise regression coefficients
    sBeta <- t(Beta[,Ind])
    # Location-wise error variance
    sSig <- Sig[Ind]
  }

  ##### Maximum likelihood
  if(type=="lm"){
    hmu <- apply(XX*sBeta, 1, sum)
    ML <- sum(dnorm(Y, hmu, sSig, log=T))
  }
  if(type %in% c("lagsarlm","errorsarlm","lmSLX")){
    ML <- val
  }



  ########## Final regression results (re-estimation)
  W_g <- listW_g <- vector(mode = "list", length = G)
  if (type=="lm"){
    for (g in 1:G) {
      W_g[[g]] <- W[Ind==g,Ind==g]
      listW_g[[g]] <- mat2listw(W_g[[g]],zero.policy = TRUE)
      if (ncol(X) == 1) {
        fit[[g]] <- lm(Y[Ind==g] ~ X[Ind==g,] - 1)
      } else {
        fit[[g]] <- lm(Y[Ind==g] ~ X[Ind==g,])
      }
    }
  }
  if (type=="lagsarlm"){
    for (g in 1:G) {
      W_g[[g]] <- W[Ind==g,Ind==g]
      listW_g[[g]] <- mat2listw(W_g[[g]],zero.policy = TRUE)
      if (ncol(X) == 1) {
        fit[[g]] <- lagsarlm(Y[Ind==g] ~ X[Ind==g,] - 1, listw=listW_g[[g]], zero.policy=zero.policy)
      } else {
        fit[[g]] <- lagsarlm(Y[Ind==g] ~ X[Ind==g,], listw=listW_g[[g]], zero.policy=zero.policy)
      }
    }
  }
  if (type=="errorsarlm"){
    for (g in 1:G) {
      W_g[[g]] <- W[Ind==g,Ind==g]
      listW_g[[g]] <- mat2listw(W_g[[g]],zero.policy = TRUE)
      if (ncol(X) == 1) {
        fit[[g]] <- errorsarlm(Y[Ind==g] ~ X[Ind==g,] - 1 , listw=listW_g[[g]], zero.policy=zero.policy)
      } else {
        fit[[g]] <- errorsarlm(Y[Ind==g] ~ X[Ind==g,], listw=listW_g[[g]], zero.policy=zero.policy)
      }
    }
  }
  if (type=="lmSLX"){
    for (g in 1:G) {
      W_g[[g]] <- W[Ind==g,Ind==g]
      listW_g[[g]] <- mat2listw(W_g[[g]],zero.policy = TRUE)
      if (ncol(X) == 1) {
        fit[[g]] <- lmSLX(Y[Ind==g] ~ X[Ind==g,] - 1 , listw=listW_g[[g]], zero.policy=zero.policy)
      } else {
        fit[[g]] <- lmSLX(Y[Ind==g] ~ X[Ind==g,], listw=listW_g[[g]], zero.policy=zero.policy)
      }
    }
  }

  ##### Change names to coefficients
  for (g in 1:G) {
    names(fit[[g]]$coefficients) <- Xnames
  }

  ##### Exit message
  if (Verbose == TRUE) {
    cat(paste0("Spatially-clustered spatial autoregressive model ended at ",round(Sys.time()),"\n"))
  }

  ########## Results
  result <- list(Results=fit, Beta=Beta, Sig=Sig, group=Ind, sBeta=sBeta, sSig=sSig, ML=ML, itr=k,
                 W_g = W_g, listW_g = listW_g)
  return(result)

}



########## SCSR_InfoCrit
##### Aims: Automatically select the optimal number of clusters based on BIC
##### Arguments:
# Y = n x 1 vector of response variable
# X = n x K matrix of explanatory/covariate variables
# W = n x n matrix with spatial weights (list2mat object)
# Wstand = binary value stating of the W matrix has to be row-standardized (TRUE) or not (FALSE, default)
# Sp = sp object with spatial weights and coordinates
# G = number of clusters
# Phi = Spatial penalty parameter
# type = model specification to be estimated
#      = "lm" for linear model without spatial effects
#      = "lagsarlm" for spatial autoregressive (SAR) model
#      = "errorsarlm" for linear model with spatial autoregressive error term (SEM, Durbin)
#      = "lm" for linear model with spatially-lagged response variable and covariates (SLX)
SCSR_InfoCrit <- function(Y, X, W, Sp, Wstand=FALSE, Phi.set=1, G.set=c(2,3,4), maxitr=100, type=c("lm","lagsarlm","errorsarlm","lmSLX"),
                          Verbose = TRUE){

  ########## Intro message
  if (Verbose == TRUE) {
    TypeName <- switch(EXPR = type,
                       "lm" = "linear model without spatial effects (lm)",
                       "lagsarlm" = "Spatial autoregressive model (SAR)",
                       "errorsarlm" = "linear model with spatial autoregressive error term (SEM)",
                       "lmSLX" = "linear model with spatially-lagged response variable and covariates (SLX)"
    )
    cat(paste0("Computing information criteria for the spatially-clustered spatial autoregressive model started at ",round(Sys.time()),"\n"))
    cat(paste0("Selected model: ",TypeName,"\n"))
  }

  ########## Setup
  Comb <- expand.grid(G = G.set,Phi = Phi.set)

  ##### Computing information criteria for each combination
  for(l in 1:dim(Comb)[1]) {
    try({
      cat(paste0("Computing combination ",l," of ",dim(Comb)[1],": g=",Comb$G[l]," and phi=",Comb$Phi[l]," started at ",round(Sys.time()),"\n"))
      fit <- SCSR_Estim(Y, X, W, Sp, G=Comb$G[l], Phi=Comb$Phi[l], maxitr=maxitr, type=type, Verbose = FALSE)
      ## Number of regression parameters (covariates + intercept + spatial parameter rho + residual variance sigma squared)
      pp <- summary(fit$Results[[1]])$parameters
      # pp <- dim(fit$Beta)[1]
      ## Total number of estimated parameters (covariates + group-wise variances)*groups
      k <- Comb$G[l]*pp
      ## Total number of observations/units
      n <- dim(fit$sBeta)[1]
      ## Clustering BIC and AIC
      # "Knee Point Detection in BIC for Detecting the Number of Clusters" (Qinpei Zhao, Ville Hautamaki & Pasi Fränti)
      Comb$BIC[l] <- -2*fit$ML + k*log(n)
      Comb$AIC[l] <- -2*fit$ML + k*2
    },
    silent=TRUE)
  }

  ##### Selection
  g_star_BIC <- Comb$G[which.min(Comb$BIC)]
  Phi_star_BIC <- Comb$Phi[which.min(Comb$BIC)]
  g_star_AIC <- Comb$G[which.min(Comb$AIC)]
  Phi_star_AIC <- Comb$Phi[which.min(Comb$AIC)]

  ##### Exit message
  if (Verbose == TRUE) {
    cat(paste0("Optimal hyperparameters (minimum BIC): g* = ",g_star_BIC," and phi* = ",Phi_star_BIC,"\n"))
    cat(paste0("Optimal hyperparameters (minimum AIC): g* = ",g_star_AIC," and phi* = ",Phi_star_AIC,"\n"))
    cat(paste0("Computing information criteria for the spatially-clustered spatial autoregressive model ended at ",round(Sys.time()),"\n"))
  }

  ##### Results
  Result <- list(IC=Comb,
                 g_star_BIC = g_star_BIC, Phi_star_BIC = Phi_star_BIC,
                 g_star_AIC = g_star_AIC, Phi_star_AIC = Phi_star_AIC)
  return(Result)
}



########## sarlogLik_i
##### Aims: log-likelihood function for the SAR model
sarlogLik_i <- function(y, X, W, beta, rho, sigma2) {
  N <- length(y)  # Number of observations
  I <- diag(N)    # Identity matrix

  # Compute the transformed y
  y_tilde <- y - rho * W %*% y

  # Compute the log determinant term
  log_det_term <- log(det(I - rho * W))

  # Compute the residuals
  residuals <- y_tilde - X %*% beta

  # Compute the log-likelihood
  log_likelihood <- - (1 / 2) * log(2 * pi * sigma2) + (log_det_term/N) -
    (1 / (2 * sigma2)) * residuals^2

  return(log_likelihood)
}



########## sarlogLik_i
##### Aims: log-likelihood function for the SEM model
semlogLik_i <- function(y, X, W, beta, lamb, sigma2) {
  N <- length(y)  # Number of observations
  I <- diag(N)    # Identity matrix
  A <- I - lamb * W

  # Compute the transformed y
  y_tilde <- A%*%y

  # Compute the log determinant term
  log_det_term <- log(det(A))

  # Compute the residuals
  residuals <- y_tilde - X %*% beta

  # Compute the log-likelihood
  log_likelihood <- - (1 / 2) * log(2 * pi * sigma2) + (log_det_term/N) -
    (1 / (2 * sigma2)) * residuals^2

  return(log_likelihood)
}



########## elbow_finder
##### Aims: Automatically select the optimal number of clusters based on elbow criterion
##### Arguments:
# x_values = m x 1 vector of integer values (usually, the number of clusters from 2 to G)
# y_values = m x 1 vector of values (usually, the criterion values) associated with the number of groups
Elbow_finder <- function(x_values, y_values) {
  # Max values to create line
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))

  # Creating straight line between the max values
  fit <- lm(max_df$y ~ max_df$x)

  # Distance from point to line
  distances <- c()
  for(i in 1:length(x_values)) {
    distances <- c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
  }

  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]

  #y_max_dist0 <- ifelse(y_max_dist==2, y_max_dist, y_max_dist-1)

  return(c(x_max_dist, y_max_dist))
}




##### Spatial Pseudo R2
SpatReg_PseudoR2 <- function(SRmodel,y) {
  SSE <- sum(residuals(model)^2)
  1-(SSE/(var(y)*(length(y)-1)))
}

##### Performance metrics
SpatReg_Perf <- function(model,y) {
  cbind(performance::performance(model = model, verbose = F),
        PseudoR2 = SpatReg_PseudoR2(model = model,y = y))
}

##### Extraction from spatial autoregression model
SpatReg_Extract <- function(model) {
  ##### Spatial parameters
  rho <- ifelse(is.null(model$rho),NA,model$rho)
  lambda <- ifelse(is.null(model$lambda),NA,model$lambda)
  pars <- c(rho,lambda)
  names(pars) <- c("rho","lambda")
  ##### Output
  return(list(SpatPars = pars))
}

MoranIStat <- function(Data,VarName,ContMat_listw,Labels,cols = c("grey25","grey60","grey85","grey45"),
                       Plot_titles = rep(NA,5), Plot_subtitles = rep(NA,5)) {

  # Data: sf object
  # VarName <- "RAPPORTOM"
  # ContMat_listw <- contnb_q_listw
  # Labels <- "LAU_NAME"
  # Labels <- NULL

  Data_moran <- Data %>%
    select(Y = any_of(VarName),
           Lab = any_of(Labels)) %>%
    mutate(Y_scaled = as.vector(scale(x = Y,center = T, scale = T)),
           Y_scaled_lag = as.vector(lag.listw(x = ContMat_listw, var = Y_scaled)),
           Quarter = case_when(Y_scaled >= 0 & Y_scaled_lag >= 0 ~ "HH",
                               Y_scaled < 0 & Y_scaled_lag >= 0 ~ "LH",
                               Y_scaled < 0 & Y_scaled_lag < 0 ~ "LL",
                               Y_scaled >= 0 & Y_scaled_lag < 0 ~ "HL"))


  ##### Global Moran's I statistic
  # https://mgimond.github.io/simple_moransI_example/?
  # Global Moran's statistic
  MoranI <- moran(x = Data_moran$Y, listw = ContMat_listw, n = length(Data_moran$Y), S0 = Szero(ContMat_listw))
  # Testing significance of Moran's I (Gaussian distribution)
  # "the excess mortality rates are randomly distributed across municipalities following a completely random process"
  MoranI_test_rand <- moran.test(x = Data_moran$Y, listw = ContMat_listw)
  MoranI_test <- moran.test(x = Data_moran$Y, listw = ContMat_listw, randomisation = FALSE)
  # Testing significance of Moran's I (Monte Carlo empirical distribution)
  MoranI_MC <- moran.mc(x = Data_moran$Y, listw = ContMat_listw, nsim = 1000)
  P1 <- plot(MoranI_MC)

  ##### Lagged variable and Moran's plot (Global I)
  # Y_lag <- lag.listw(x = contnb_q_listw, var = Data_moran$Y)
  # plot(Y_lag ~ Data$RAPPORTOMA, pch=16, asp=1)
  # M1 <- lm(Y_lag ~ Data$RAPPORTOMA)
  # abline(M1, col="blue")
  # summary(M1)
  P2 <- moran.plot(x = Data_moran$Y, listw = ContMat_listw, labels = Data_moran$Lab, pch = 19, quiet = F,
                   main = "Moran's scatterplot", xlab = "Variable", ylab = "Lagged variable")

  ##### Map plot (Global I)
  brks <- c("HH","LH","LL","HL")
  P3 <- Data_moran %>%
    ggplot()+
    theme_void() +
    geom_sf(colour = "white") +
    geom_sf(aes(fill= Quarter), size = 0.2, col="black") +
    scale_fill_manual(breaks = brks, values = cols) +
    labs(title = "Moran's quarters") +
    theme(title = element_text(size = 25, face = "bold"))

  ##### Local Moran's I
  alpha <- 0.05
  LocMoranI <- localmoran(x = Data_moran$Y, listw = ContMat_listw, zero.policy = TRUE, na.action = na.omit)
  LocMoranI_MC <- localmoran_perm(x = Data_moran$Y, listw = ContMat_listw, nsim = 1000, zero.policy = TRUE, na.action = na.omit)
  Data_moran <- Data_moran %>%
    mutate(LISA_Moran = as.vector(LocMoranI[,1]),
           LISA_Moran_pv = as.vector(LocMoranI[,5]),
           LISA_Moran_MC = as.vector(LocMoranI_MC[,1]),
           LISA_Moran_MC_pv = as.vector(LocMoranI_MC[,5]),
           Quarter_LISA = case_when(LISA_Moran_pv > alpha ~ "Not significant",
                                    TRUE ~ Quarter),
           Quarter_LISA_MC = case_when(LISA_Moran_MC_pv > alpha ~ "Not significant",
                                       TRUE ~ Quarter))

  ##### Map plot (Local I with significance)
  P4 <- Data_moran %>%
    ggplot()+
    theme_void() +
    geom_sf(colour = "white") +
    geom_sf(aes(fill= Quarter_LISA), size = 0.2, col="black") +
    scale_fill_manual("LISA quarters",breaks = c("Not significant","HH","LH","LL","HL"), values = c("white",cols)) +
    labs(
      #title = ifelse(is.na(Plot_titles[4]),"Local Moran's quarters",Plot_titles[4]),
      subtitle = ifelse(is.na(Plot_subtitles[4]),"Monte Carlo permutations",Plot_subtitles[4])) +
    theme(title = element_text(size = 25, face = "bold"))

  P5 <- Data_moran %>%
    ggplot()+
    theme_void() +
    geom_sf(colour = "white") +
    geom_sf(aes(fill= Quarter_LISA), size = 0.2, col="black") +
    scale_fill_manual("LISA quarters",breaks = c("Not significant","HH","LH","LL","HL"), values = c("white",cols)) +
    labs(title = ifelse(is.na(Plot_titles[5]),"Local Moran's quarters",Plot_titles[5]),
         subtitle = ifelse(is.na(Plot_subtitles[5]),"Monte Carlo permutations",Plot_subtitles[5])) +
    theme(title = element_text(size = 25, face = "bold"))

  ##### Moran’s I spatial correlogram
  contnb_q <- spdep::poly2nb(Data_moran, queen = TRUE, row.names = Data_moran$Lab)
  MoranCorr <- sp.correlogram(contnb_q, Data_moran$Y, order = 10, method = "I", style = "B",zero.policy = TRUE)
  # plot(MoranCorr,main = "Moran's spatial correlogram", xlab = "Neighbors lag", ylab = "Variable")
  P6_data <- MoranCorr$res
  P6_data <- data.frame(Lag = 1:dim(P6_data)[1], P6_data)
  colnames(P6_data) <- c("Lag","Estimate","Expectation","Variance")
  P6 <- P6_data %>%
    ggplot(aes(x=Lag, y=Estimate), col = "black") +
    geom_point()+
    geom_errorbar(aes(ymin=Estimate-2*sqrt(Variance), ymax=Estimate+2*sqrt(Variance)),
                  width=.2, position=position_dodge(0.05)) +
    geom_hline(yintercept = 0, col = "red", size = 1) +
    labs(title = "Moran's spatial correlogram", x = "Lags", y = "Moran's I estimate") +
    scale_x_discrete(limits = 1:10) +
    scale_y_continuous(breaks = round(seq(from = min(P6_data$Estimate)-0.05, to = max(P6_data$Estimate)+0.05, by = 0.05),2)) +
    theme_bw()

  return(list = list(Data_moran = Data_moran,
                     MoranI = MoranI,
                     MoranI_test_rand = MoranI_test_rand,
                     MoranI_test = MoranI_test,
                     MoranI_MC = MoranI_MC,
                     LocMoranI = LocMoranI,
                     LocMoranI_MC = LocMoranI_MC,
                     MoranCorr = MoranCorr,
                     P1 = P1,P2 = P2, P3 = P3,P4 = P4, P5 = P5, P6 = P6))
}

