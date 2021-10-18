#--------- Supplementary R code for the paper "Longitudinal Functional Principal Component Analysis ------------------#

##########################################################################################
### R function LFPCA implementing LFPCA for model (2)                                  ###
##########################################################################################

LFPCA <- function(Y,               # an n x D matrix with rows=subject-visits and columns=locations (along curve)
                  subject,         # a vector of length n containing subject identifiers for rows of Y
                  T,               # a vector of length n containing the covariate for the random slope
                  time = NULL,
                  L = 0.90,        # the pre-specified level of variance explained, determines number of components   
                  N.X = NA,        # the number of components to keep for X; N.X and N.U override L if not NA
                  N.U = NA,        # the number of components to keep for U; N.X and N.U override L if not NA
                  smoothing = FALSE,  # smooth=TRUE adds smoothing of the covariance matrices not done for smooth=FALSE  
                  smooth.method = "sf",
                  bf = 10          # the number of basis functions used per dimension for all smooths
                  ){
    require(mgcv)
    require(Matrix)
    
    ### checks for consistency of input ###
    if (nrow(Y) != length(subject)){
       stop("The number of rows in Y needs to agree with the length of the subject vector")
    }
    if (nrow(Y) != length(T)){
       stop("The number of rows in Y needs to agree with the length of the time vector")
    }
    if (!is.na(L)){
       if (L > 1 | L < 0){
          stop("The level of explained variance needs to be between 0 and 1")
       }
    }
    if (is.na(L) & (is.na(N.X) | is.na(N.U) | N.X < 0 | N.U < 0)){
       stop("If L is NA, both N.X and N.U need to be integers")
    }
    if (is.na(N.X) & !is.na(N.U)){
       warning("As N.X is missing, N.U will not be used. Default to variance explained.")
    }
    if (is.na(N.U) & !is.na(N.X)){
       warning("As N.U is missing, N.X will not be used. Default to variance explained.")
    }
    if (!is.na(N.X)){
       if (N.X > floor(N.X)){
          warning("N.X has to be an integer. Will use rounded value.")
          N.X <- round(N.X)
       }
    }
    if (!is.na(N.U)){
       if (N.U > floor(N.U)){
          warning("N.U has to be an integer. Will use rounded value.")
          N.U <- round(N.U)
       }
    }
    if (!is.na(L) & (!is.na(N.X) | !is.na(N.U))){
       warning("N.X and N.U will override choice of variance explained L")
    }
  
    Y_raw = Y
    if(smoothing == TRUE & smooth.method == "sf"){
      Y_sm = matrix(NA,ncol=ncol(Y),nrow=nrow(Y))
      tY = 1:ncol(Y_sm)
      
      for(i in 1:nrow(Y_sm)){
        Y_sm[i,] = gam(Y[i,] ~ s(tY,k = bf))$fitted.values
      }
      Y = Y_sm
    }
    mu = apply(Y, 2,mean)
       
    ### set up ###
    D        <- ncol(Y)                        # number of points per curve
    n        <- nrow(Y)                        # overall number of visits
    time     <- (T-mean(T)) / sqrt(var(T))     # standardize the time variable
    #time     <- T
    if (is.null(time)){
      time     <- (T-mean(T)) / sqrt(var(T))
      #time = T
    }
    d.vec    <- rep(1:D, each = n)                
    time.vec <- rep(time, D)
    J.vec    <- sapply(unique(subject[order(subject)]), function(subj){sum(subject == subj)})
                                               # number of visits for each subject
    I        <- length(J.vec)                  # number of subjects
    m        <- sum(J.vec^2)                   # overall number of visit pairs within subjects
                                                             
    ### estimate mean function eta(d, T_ij) and overall mean eta0(d) ###
    #gam1       <- gamm(as.vector(Y) ~ te(d.vec, time.vec, k = bf))  # overall fixed effects surface eta(d, T)
    #gam0       <- gamm(as.vector(Y) ~ s(d.vec, k = bf))             # time invariant mean function eta(d) for plotting
    #eta.matrix <- matrix(predict(gam0$gam), n, D)                 # mean function eta(d, T) at grid points for plotting
    eta.matrix <- matrix(rep(mu,n), n, D,byrow = TRUE) 
    Y.tilde    <- Y - eta.matrix                                  # centered Y (residuals)
                                                             
    ### estimate covariance functions using least squares ###
    G.0 <- G.01 <- H.01 <- G.1 <- G.U <- matrix(0, D, D)          # set up empty covariance matrices
    diago <- rep(NA, D)                                           # set up empty diagonal for K.U
    i1 <- function(i){                                            # find all j for indicator vectors for all i-(j,k) pairs
       before <- sum(J.vec[1:(i-1)]) * (i > 1)
       rep(before + (1:J.vec[i]), each = J.vec[i])
    } 
    i2 <- function(i){                                            # find all k for indicator vectors for all i-(j,k) pairs
       before <- sum(J.vec[1:(i-1)]) * (i > 1)
       rep(before + (1:J.vec[i]), J.vec[i])
    }
    ind1 <- as.vector(unlist(sapply(1:I, i1)))                    # indicator vectors for j in all i-(j,k) pairs
    ind2 <- as.vector(unlist(sapply(1:I, i2)))                    # indicator vectors for k in all i-(j,k) pairs
    
    X <- cbind(rep(1, m), time[ind2], time[ind1], time[ind1] * time[ind2], (ind1 == ind2) * 1)
                                            # set up design matrix for linear regression according to (5)
    c <- matrix(unlist(lapply(1:(D-1), FUN = function(s){Y.tilde[ind1, s] * Y.tilde[ind2, (s+1):D]})), length(ind1), D * (D-1) / 2)
                                            # set up outcome vector (empirical covariances) for linear regression      
    beta <- solve(crossprod(X), t(X) %*% c)   # estimate covariances K.U(s,t) and K.X(s,t) (containing K.0(s,t), K.1(s,t), K.01(s,t) and K.01(t,s)) for all s != t using linear regression based on (5)
    
    Xss <- cbind(X[, 1], X[, 2] + X[, 3], X[, 4], X[, 5])[(ind1 >= ind2), ] # set up design matrix according to (5)
    css <- sapply(1:D, FUN = function(s){Y.tilde[ind1[(ind1 >= ind2)], s] * Y.tilde[ind2[(ind1 >= ind2)], s]}) # set up outcome vector
    betass <- solve(crossprod(Xss), t(Xss) %*% css)   # estimate covariances K.U(s,s) and K.X(s,s) (containing K.0(s,s), K.1(s,s) and K.01(s,s)) for all s using linear regression based on (5)
    
    G.0[outer(1:D, 1:D, FUN = function(s, t){(s > t)})]  <- beta[1, ]   # estimates for K.0(s,t) are in the 1st column
        G.0 <- G.0 + t(G.0)                                         # symmetry constraints yield K.0(s,t) = K.0(t,s)
    G.1[outer(1:D, 1:D, FUN = function(s, t){(s > t)})]  <- beta[4, ]   # estimates for K.1(s,t) are in the 4th column
        G.1 <- G.1 + t(G.1)                                         # symmetry constraints yield K.1(s,t) = K.1(t,s)
    G.U[outer(1:D, 1:D, FUN = function(s, t){(s > t)})]  <- beta[5, ]   # estimates for K.U(s,t) are in the 5th column
        G.U <- G.U + t(G.U)                                         # symmetry constraints yield K.U(s,t) = K.U(t,s)
    H.01[outer(1:D, 1:D, FUN = function(s, t){(s > t)})] <- beta[2, ]   # estimates for K.10(s,t) are in the 2nd column
    G.01[outer(1:D, 1:D, FUN = function(s, t){(s > t)})] <- beta[3, ]   # estimates for K.01(s,t) are in the 3rd column
        G.01 <- G.01 + t(H.01)                                      # symmetry constraints yield K.01(s,t) = K.10(t,s)
    diag(G.0)  <- betass[1, ]                                       # estimates for K.0(s,s) are in the 1st column
    diag(G.1)  <- betass[3, ]                                       # estimates for K.1(s,s) are in the 3rd column
    diag(G.01) <- betass[2, ]                                       # estimates for K.01(s,s) are in the 2nd column
    diago      <- betass[4, ]                                       # estimates for K.U(s,s)+sigma^2 are in the 4th column
    diag(G.U)  <- rep(NA, D)                                        # do not use diagonal K.U(s,s)+sigma^2 in smoothing K.U
    
    ### smoothing of covariance functions, estimation of sigma^2 ###                                        
    row.vec <- rep(1:D, each = D)     # set up row variable for bivariate smoothing                             
    col.vec <- rep(1:D, D)          # set up column variable for bivariate smoothing
    if (smoothing == TRUE & smooth.method == "sc"){      # if smoothing is selected (corresponds to method described in the paper):                       
       K.0  <- matrix(predict(gamm(as.vector(G.0)  ~ te(row.vec, col.vec, k = bf))$gam), D, D) # smooth K.0
          K.0 <- (K.0 + t(K.0)) / 2                                                            # after smoothing, symmetrize K.0 
       K.1  <- matrix(predict(gamm(as.vector(G.1)  ~ te(row.vec, col.vec, k = bf))$gam), D, D) # smooth K.1
          K.1 <- (K.1 + t(K.1)) / 2                                                            # after smoothing, symmetrize K.1
       K.01 <- matrix(predict(gamm(as.vector(G.01) ~ te(row.vec, col.vec, k = bf))$gam), D, D) # smooth K.01
       K.U  <- matrix(predict(gamm(as.vector(G.U)  ~ te(row.vec, col.vec, k = bf))$gam,
                           newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)  # smooth K.U
          K.U <- (K.U + t(K.U)) / 2                                                            # after smoothing, symmetrize K.U
    }
    else {            # if no smoothing is selected (faster):
          K.U <- G.U    # do not smooth K.U (off-diagonal) 
          K.0 <- G.0    # do not smooth K.0
          K.01 <- G.01  # do not smooth K.01
          K.1 <- G.1    # do not smooth K.1
          diag(K.U) <- predict(gamm(as.vector(K.U) ~ te(row.vec, col.vec, k = bf))$gam, newdata = data.frame(row.vec = 1:D, col.vec = 1:D))                   # only separate diagonal K.U(s,s) + sigma^2 into K.U(s,s) and sigma^2 using bivariate smoothing 
    }
    
    ## comment1: change to K.X/2
    K.X <- rbind(cbind(K.0, K.01), cbind(t(K.01), K.1))   # put together K.X, containing K.0, K.01, K.10 and K.1          
    #sigma2.hat <- max(mean(diago - diag(K.U)), 0)         # estimate sigma^2 as mean difference between diagonal
    sigma2.hat <- 0                                                   # K.U(s,s) + sigma^2 and smoothed K.U(s,s) (0, if smaller 0 otherwise)
    
    ### estimate eigenfunctions, compute variance explained, N.X and N.U ###
    lambda.hat <- eigen(K.X, symmetric = TRUE, only.values = TRUE)$values # eigenvalues of K.X yield estimates for the lambda_k
    nu.hat     <- eigen(K.U, symmetric = TRUE, only.values = TRUE)$values # eigenvalues of K.U yield estimates for the nu_k
    
    total.variance <- sum(lambda.hat * (lambda.hat > 0)) + sum(nu.hat * (nu.hat > 0)) + sigma2.hat
               # the total average variance is the sum of all variance terms lambda_k, nu_k and sigma^2 according to Lemma 1 
    if (is.na(N.X) | is.na(N.U)){  # if no values for the numbers of principal components N.X and N.U are specified:
       prop <- N.X <- N.U <- 0
       while(prop < L){              # add components for X or U with decreasing variance,
                                   # until level L of explained average variance is reached
          if (lambda.hat[N.X + 1] >= nu.hat[N.U + 1]){
             N.X <- N.X + 1
          }
          else {
                N.U <- N.U + 1
          }
          prop <- (sum(lambda.hat[1:N.X]) + sum(nu.hat[1:N.U])) / total.variance  # update explained average variance
       }
    }
    explained.variance <- (sum(lambda.hat[1:N.X]) + sum(nu.hat[1:N.U])) / total.variance  # explained average variance 
                                                                    
    lambda.hat <- lambda.hat[1:N.X]                          # keep first N.X values in lambda.hat
    nu.hat     <- nu.hat[1:N.U]                              # keep first N.U values in nu.hat
    phi.X      <- eigen(K.X, symmetric = TRUE)$vectors[, 1:N.X]  # eigenvectors of K.X yield estimated for eigenfunctions phi.X_k
    phi.U      <- eigen(K.U, symmetric = TRUE)$vectors[, 1:N.U]  # eigenvectors of K.U yield estimated for eigenfunctions phi.U_k
    ###     The eigen vectosrs are unique only up to a change of signs.
    ###     Select the signs of eigenfunctions so that the integration over the domain
    ###     is non-negative
    for(i in 1:N.X) {
      v2 = phi.X[,i]
      tempsign = sum(v2)
      phi.X[,i] = ifelse(tempsign<0, -1,1) * v2
    }
    for(i in 1:N.U) {
      v2 = phi.U[,i]
      tempsign = sum(v2)
      phi.U[,i] = ifelse(tempsign<0, -1,1) * v2
    }
    phi.0      <- phi.X[1:D, ]                                # phi.X_k = (phi.0_k, phi.1_k)
    phi.1      <- phi.X[(D+1):(2*D), ]
                                                                                  
    
    ### estimate scores ###
    subject.ind <- unlist(sapply(subject, function(su){which(unique(subject[order(subject)]) == su)}))
    Z.X         <- kronecker(Diagonal(I)[subject.ind, ]                   , phi.0) +
                   kronecker(Diagonal(n, time) %*% Diagonal(I)[subject.ind, ], phi.1)        # set up the matrices for BLUP   
    Z.U         <- kronecker(Diagonal(n)                                 , phi.U)        #   estimation using the Woodbury 
    Z           <- cbind(Z.X, Z.U)                                                       #   formula
    D.inv       <- Diagonal(N.X*I + N.U*n, c(rep(1 / lambda.hat, I), rep(1 / nu.hat, n)))
    
    ## comment2: change to BLUP estimation (sigma2.hat = 0)
    b.hat       <- solve(crossprod(Z), t(Z) %*% as.vector(t(Y.tilde)))
    #b.hat        <- solve(crossprod(Z%*%sqrt(D.inv)), t(Z %*% D.inv) %*% as.vector(t(Y.tilde)))
    #b.hat       <- solve(crossprod(Z) + sigma2.hat * D.inv, t(Z) %*% as.vector(t(Y.tilde))) # estimate scores according to Thm. 2
    xi.hat      <- t(matrix(b.hat[1:(N.X*I)], N.X, I))                                    #   b.hat contains first xi.hat, 
    zeta.hat    <- t(matrix(b.hat[(N.X*I) + (1:(N.U*n))], N.U, n))                          #   then zeta.hat
    
    ### return results ###
    results <- list(Y = Y_raw, Y_sm = Y, subject = subject, time = time, #eta = gam1$gam, 
                    #eta0 = gam0$gam, eta.matrix = eta.matrix,
                    eta0 = mu, eta.matrix = eta.matrix,
                    phi.0 = phi.0, phi.1 = phi.1, phi.U = phi.U, K.X = K.X, K.U = K.U,
                    sigma2 = sigma2.hat, lambda = lambda.hat, nu = nu.hat,
                    xi = xi.hat, zeta = zeta.hat, N.X = N.X, N.U = N.U,
                    L = L, totvar = total.variance, exvar = explained.variance
                   )
    return(results)
}




##########################################################################################
### R function plot.LFPCA plotting results from fitting model (2) using function LFPCA ###
##########################################################################################

plot.LFPCA <- function(res,                      # a results file from LFPCA
                       outpdf,                   # a pdf file to which to plot the results
                       group = rep(1, nrow(res$Y)) # a grouping variable for the boxplots of the eigenscores (optional)
                       ){                                        
                                        
    D      <- ncol(res$Y)                                      # the number of points per curve
    group2 <- unique(cbind(res$subject, group))[, 2]            # the grouping variable per subject (not visit)
    eta    <- predict(res$eta0, newdata = data.frame(d.vec = 1:D))  # the time invariant mean function for plotting
    lu     <- range(cbind(matrix(eta, D, res$N.X) - 2*res$phi.0 %*% diag(sqrt(res$lambda)),  # limits for plotting
                          matrix(eta, D, res$N.X) + 2*res$phi.0 %*% diag(sqrt(res$lambda)),
                          matrix(eta, D, res$N.X) - 2*res$phi.1 %*% diag(sqrt(res$lambda)),
                          matrix(eta, D, res$N.X) + 2*res$phi.1 %*% diag(sqrt(res$lambda)),
                          matrix(eta, D, res$N.U) - 2*res$phi.U %*% diag(sqrt(res$nu    )),
                          matrix(eta, D, res$N.U) + 2*res$phi.U %*% diag(sqrt(res$nu    ))))
    
    plot.ef <- function(ev, ef, yl, var){  # for given eigenfunctions and corresponding eigenvalues:
                  plot(1:D, eta, type = "l", xlab = "d", ylab = yl, ylim = lu,
                       main = paste(round(var), "% variance", sep = ""))
                                            # plot mean function
                  points(eta + 2*sqrt(ev)*ef, pch = "+")   # plus / minus 2 times the standard deviation times the 
                  points(eta - 2*sqrt(ev)*ef, pch = "-")   #    corresponding eigenfunction
                                     }
    pdf(file = outpdf)
      plot(res$eta, main = "Estimated mean function", xlab = "d", ylab = "T")   # plot estimated mean profile eta(d,T)
      plot(res$eta0, main = "Estimated time-constant mean function", xlab = "d", ylab = expression(eta[0](d)))
                                                                       # plot estimated time-constant mean function eta(d)
      par(mfrow = c(ceiling(sqrt(res$N.X)), ceiling(res$N.X / ceiling(sqrt(res$N.X)))), mar = c(5.1, 5.1, 1.1, 1.1))
      for (k in 1:res$N.X){  # for each eigenfunction, plot boxplots of subject-specific scores by the grouping variable
        boxplot(res$xi[, k] ~ group2, ylab = eval(substitute(expression(hat(xi)[i][j]), list(j = k))))
        abline(h = 0, col = 8)
      }
      par(mfrow = c(ceiling(sqrt(res$N.X)), ceiling(res$N.X / ceiling(sqrt(res$N.X)))), mar = c(5.1, 5.1, 1.1, 1.1))
      for (k in 1:res$N.X){ # for each phi.0_k, plot mean function +/- 2 the standard deviation sqrt(lambda_k) times phi.0_k
        plot.ef(ev = res$lambda[k], ef = res$phi.0[, k], yl = substitute(hat(phi)[k]^0, list(k = k)),
                res$phi.0[, k] %*% res$phi.0[, k] * res$lambda[k] / res$totvar * 100)
      }
      par(mfrow = c(ceiling(sqrt(res$N.X)), ceiling(res$N.X / ceiling(sqrt(res$N.X)))), mar = c(5.1, 5.1, 1.1, 1.1))
      for (k in 1:res$N.X){ # for each phi.1_k, plot mean function +/- 2 the standard deviation sqrt(lambda_k) times phi.1_k
        plot.ef(ev = res$lambda[k], ef = res$phi.1[, k], yl = substitute(hat(phi)[k]^1, list(k = k)),
                res$phi.1[, k] %*% res$phi.1[, k] * res$lambda[k] / res$totvar * 100)
      }
      par(mfrow = c(ceiling(sqrt(res$N.X)), ceiling(res$N.X / ceiling(sqrt(res$N.X)))), mar = c(5.1, 5.1, 1.1, 1.1))
      for (k in 1:res$N.X){ # plot each phi.0_k
        plot(res$phi.0[, k], type = "l", xlab = "d", ylab = substitute(hat(phi)[k]^0, list(k = k)),
             main = paste(round(res$phi.0[, k] %*% res$phi.0[, k] * res$lambda[k] / res$totvar * 100), "% variance", sep = ""))
      }
      par(mfrow = c(ceiling(sqrt(res$N.X)), ceiling(res$N.X / ceiling(sqrt(res$N.X)))), mar = c(5.1, 5.1, 1.1, 1.1))
      for (k in 1:res$N.X){ # plot each phi.1_k
        plot(res$phi.1[, k], type = "l", xlab = "d", ylab = substitute(hat(phi)[k]^1, list(k = k)),
             main = paste(round(res$phi.1[, k] %*% res$phi.1[, k] * res$lambda[k] / res$totvar * 100), "% variance", sep = ""))
      }
      par(mfrow = c(ceiling(sqrt(res$N.U)), ceiling(res$N.U / ceiling(sqrt(res$N.U)))), mar = c(5.1, 5.1, 1.1, 1.1))
      for (k in 1:res$N.U){ # for each eigenfunction, plot boxplots of subject-visit-specific scores by the grouping variable
        boxplot(res$zeta[, k] ~ group, ylab = eval(substitute(expression(hat(zeta)[ij][l]), list(l = k))))
        abline(h = 0, col = 8)
      }
      par(mfrow = c(ceiling(sqrt(res$N.U)), ceiling(res$N.U / ceiling(sqrt(res$N.U)))), mar = c(5.1, 5.1, 1.1, 1.1))
      for (k in 1:res$N.U){ # for each phi.U_k, plot mean function +/- 2 the standard deviation sqrt(nu_k) times phi.U_k
        plot.ef(ev = res$nu[k], ef = res$phi.U[, k], yl = substitute(hat(phi)[k]^U, list(k = k)), res$nu[k] / res$totvar * 100)
      }
      par(mfrow = c(ceiling(sqrt(res$N.U)), ceiling(res$N.U / ceiling(sqrt(res$N.U)))), mar = c(5.1, 5.1, 1.1, 1.1))
      for (k in 1:res$N.U){ # plot each phi.U_k
        plot(res$phi.U[, k], type = "l", xlab = "d", ylab = substitute(hat(phi)[k]^U, list(k = k)),
             main = paste(round(res$nu[k] / res$totvar * 100), "% variance", sep = ""))
      }
      par(mfrow = c(1, 1))
      # plot all variances in one plot
      plot(res$lambda, ylim = c(0, max(res$lambda)), ylab = "Estimated variance components", xaxt = "n", xlab = "") 
      points(res$nu, pch = "+")
      points(res$sigma2, pch = 15)
      legend("topright", pch = c(1, 3, 15), legend = expression(hat(lambda), hat(nu), hat(sigma)^2))
    dev.off()
}






