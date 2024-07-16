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
#      = "lm" for linear model with spatially-lagged response variable and covariates (SLX)
Cl_spatialReg <- function(Y, X, W, Sp, Wstand=FALSE, G=2, Phi=1, maxitr=100, type=c("lm","lagsarlm","errorsarlm","lmSLX")){
  
  ## Preparations
  ep <- 10^(-5)      # convergence criterion 
  if (!is.null(X)){  X <- as.matrix(X) }
  
  n <- length(Y) # number of samples
  
  if (type=="lm"){
    if (is.null(X)) { p <- 1
    } else {
      p <- dim(X)[2]+1    # number of regression coefficients
    }
  }
  if (type=="lagsarlm"){
    if (is.null(X)) { p <- 2
    } else {
      p <- dim(X)[2]+2
    }
  }
  if (type=="errorsarlm"){
    if (is.null(X)) { p <- 2
    } else {
      p <- dim(X)[2]+2
    }
  }
  if (type=="lmSLX"){
    if (Wstand==FALSE){ p <- 2*(dim(X)[2]+1) } else { p <- (2*dim(X))+1}
  }
  
  if (is.null(X)) { 
    XX <- rep(1,n)
  } else {
    XX <- as.matrix( cbind(1,X) )
  }
  
  W <- as(W, "sparseMatrix")
  nmax <- function(x){ max(na.omit(x)) }   # new max function 
  
  ## Initial values
  M <- 20  # the number of initial values of k-means
  WSS <- c()
  CL <- list()
  for(k in 1:M){
    CL[[k]] <- kmeans(Sp, G)
    WSS[k] <- CL[[k]]$tot.withinss
  }
  Ind <- CL[[which.min(WSS)]]$cluster
  Pen <- rep(0, G)
  Beta <- matrix(0, p, G)
  dimnames(Beta)[[2]] <- paste0("G=",1:G)
  Sig <- rep(1, G)    
  
  ## iterative algorithm 
  val <- 0
  mval <- 0
  fit <- list()
  
  for(k in 1:maxitr){
    
    cat("*")
    cval <- val
    
    ## penalty term
    Ind.mat <- matrix(0, n, G)
    for(g in 1:G){
      Ind.mat[Ind==g, g] <- 1
    }
    Ind.mat <- as(Ind.mat, "sparseMatrix")
    Pen <- W%*%Ind.mat     # penalty term
    
    ## model parameters (clustered case)
    for(g in 1:G){
      if(length(Ind[Ind==g])>p+1){
        if(type=="lm"){
          fit[[g]] <- lm(Y[Ind==g]~X[Ind==g,])
          Beta[,g] <- as.vector( coef(fit[[g]]) )
          resid <- Y-as.vector(XX%*%Beta[,g])
          Sig[g] <- sqrt(mean(resid[Ind==g]^2))
          Sig[g] <- max(Sig[g], 0.1)
        }
        if (type=="lagsarlm"){
          XXX <- cbind(as.vector(W%*%Y),XX)
          listWg <- mat2listw(W[Ind==g,Ind==g])
          fit[[g]] <- lagsarlm(Y[Ind==g]~X[Ind==g,], listw=listWg, zero.policy=TRUE)
          Beta[,g] <- as.vector( coef(fit[[g]]) )
          Sig[g] <- fit[[g]]$s2
          Sig[g] <- max(Sig[g], 0.1)
        }
        if (type=="errorsarlm"){
          listWg <- mat2listw(W[Ind==g,Ind==g])
          fit[[g]] <- errorsarlm(Y[Ind==g]~X[Ind==g,], listw=listWg, zero.policy=TRUE)
          Beta[,g] <- as.vector( coef(fit[[g]]) )
          Sig[g] <- fit[[g]]$s2
          Sig[g] <- max(Sig[g], 0.1)
        }
        if (type=="lmSLX"){
          if (Wstand==FALSE){
            WI <- lag.listw(mat2listw(W), XX[,1], zero.policy = TRUE)
            WX <- spatialreg:::create_WX(XX[,-1],mat2listw(W), zero.policy=TRUE)
            XXX <- cbind(XX,WI,WX)
          } else {
            WX <- spatialreg:::create_WX(XX[,-1],mat2listw(W), zero.policy=TRUE)
            XXX <- cbind(XX,WX)
          }
          listWg <- mat2listw(W[Ind==g,Ind==g])
          fit[[g]] <- lmSLX(Y[Ind==g]~as.matrix(X[Ind==g,]), listw=listWg, zero.policy=TRUE)
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
    val <- sum(unlist(lapply(fit,function(x){as.numeric(x$LL)})))
    dd <- abs(cval-val)/abs(val)
    mval <- max(mval, cval)
    if( dd<ep | abs(mval-val)<ep ){ break }
  }
  
  ## varying parameters
  sBeta <- t(Beta[,Ind])
  sSig <- Sig[Ind]     # location-wise error variance
  
  ## maximum likelihood 
  if(type=="lm"){
    hmu <- apply(XX*sBeta, 1, sum)
    ML <- sum( dnorm(Y, hmu, sSig, log=T) ) 
  }
  if(type=="lagsarlm"){
    ML <- val
  }
  if(type=="errorsarlm"){
    ML <- val 
  }
  if(type=="lmSLX"){
    ML <- val
  }
  
  ## Final regression results
  
  if (type=="lm"){
    for (g in 1:G) {
      fit[[g]] <- lm(Y[Ind==g]~X[Ind==g,])
    }
  }
  if (type=="lagsarlm"){
    for (g in 1:G) {
      listWg <- mat2listw(W[Ind==g,Ind==g])
      fit[[g]] <- lagsarlm(Y[Ind==g]~X[Ind==g,], listw=listWg, zero.policy=TRUE)
    }
  }
  if (type=="errorsarlm"){
    for (g in 1:G) {
      listWg <- mat2listw(W[Ind==g,Ind==g])
      fit[[g]] <- errorsarlm(Y[Ind==g]~X[Ind==g,], listw=listWg, zero.policy=TRUE)
    }
  }
  if (type=="lmSLX"){
    for (g in 1:G) {
      listWg <- mat2listw(W[Ind==g,Ind==g])
      fit[[g]] <- lmSLX(Y[Ind==g]~X[Ind==g,], listw=listWg, zero.policy=TRUE)
    }
  }
  
  ## Results
  result <- list(Results=fit, Beta=Beta, Sig=Sig, group=Ind, sBeta=sBeta, sSig=sSig, ML=ML, itr=k)
  return(result)
  
}



########## G_select
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
G_select <- function(Y, X, W, Sp, Wstand=FALSE, Phi=1, maxitr=50, G.set=NULL, type=c("lm","lagsarlm","errorsarlm","lmSLX")){
  
  ## Preparations
  if(is.null(G.set)){
    G.set <- seq(2, 4, by=1)
    }
  L <- length(G.set)
  
  ## computing information criteria
  BIC <- c()
  for(l in 1:L) {
    try({
      cat(paste0("g=",l+1))
      fit <- Cl_spatialReg(Y, X, W, Sp, G=G.set[l], Phi=Phi, maxitr=maxitr, type=type)
      ## Number of regression parameters (covariates + intercept + spatial parameter rho + residual variance sigma squared)
      pp <- summary(fit$Results[[1]])$parameters
      # pp <- dim(fit$Beta)[1]
      ## Total number of estimated parameters (covariates + group-wise variances)*groups
      # k <- (pp + G.set[l])*G.set[l]
      k <- G.set[l]*pp
      ## Total number of observations/units
      n <- dim(fit$sBeta)[1]
      ## Clustering BIC
      # "Knee Point Detection in BIC for Detecting the Number of Clusters" (Qinpei Zhao, Ville Hautamaki & Pasi Fränti)
      BIC[l] <- -2*fit$ML + k*log(n)
    },
    silent=TRUE)
  }
  
  ## selection
  hG <- G.set[which.min(BIC)]
  
  ## result
  Result <- list(BIC=BIC, G=hG)
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
elbow_finder <- function(x_values, y_values) {
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
