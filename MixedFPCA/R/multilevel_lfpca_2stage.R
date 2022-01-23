LFPCA_3_level = function(Y = NULL, id = NULL, visit=NULL, day = NULL, 
                         time = NULL,
                         ##twoway = FALSE,  ## only nested model considered here
                         cov.method = c("m1","m2"), pve = 0.9,
                         K1 = NA, K2 = NA, K3 = NA,
                         tlength = ncol(Y),
                         smoothing = FALSE, smooth.method = c("sf","sc") ,nk = 15){
  # 1. Preliminary check on all possible mistakes ---------------------------
  stopifnot((!is.null(Y) && !is.null(id)))
  if(length(id) != nrow(Y)){
    stop("number of id visits should be the same as number of rows of Y")
  }
  
  cov.method = match.arg(cov.method)
  if(smoothing == TRUE){
    tlength = 1
  }
  
  if(smoothing  == FALSE){
    smooth.method = "nothing"
    nk = NULL
  }
  
  if(smoothing == TRUE & is.null(smooth.method)){
    stop("choosing the smoothing method")
  }
  if(smoothing == TRUE){
    cov.method = match.arg(cov.method)
  }
  
  Y_raw = Y
  if(smoothing == TRUE & smooth.method == "sf"){
    Y_sm = matrix(NA,ncol=ncol(Y),nrow=nrow(Y))
    tY = 1:ncol(Y_sm)
    
    # for(i in 1:nrow(Y_sm)){
    #   Y_sm[i,] = gam(Y[i,] ~ s(tY,k = nk),na.action = na.omit)$fitted.values
    # }
    Y_sm = t(apply(Y, 1, function(x) gam(x~s(tY,k = nk),na.action = na.omit)$fitted.values))
    Y = Y_sm
  }
  
  # 2. Check the dimension and preliminary ----------------------------------
  uni.id = unique(id) ## unique id
  uni.visit = unique(visit) ## unique visit
  id.visit <- paste0(id,"_",visit) ## measure subject-visit inputs
  uni.id.visit <- unique(id.visit)
  if (!is.null(day)){
    day = as.integer(factor(day))
  }else{ day = as.integer(ave(id.visit, id.visit, FUN=seq_along))}
  
  Y.df = data.frame(id = id, visit = visit, day = day, id.visit = id.visit)
  Y.df.visit = Y.df %>% group_by(id, visit) %>% slice(1) %>% as.data.frame()
  #J.vec = sapply(unique(Y.df$id[order(Y.df$id)]), function(subj){sum(Y.df$id == subj)})
  J.vec = sapply(unique(Y.df.visit$id[order(Y.df.visit$id)]), function(subj){sum(Y.df.visit$id == subj)})
  K.vec = sapply(unique(id.visit[order(id.visit)]), function(subj){sum(id.visit == subj)})
  id.visit.vec = sapply(unique(id.visit), function(subj){sum(id.visit == subj)})
  mm        <- sum(J.vec^2)
  Y.df$Y = Y
  
  if (is.null(time)){
    time     <- (visit-mean(visit)) / sqrt(var(visit))
    #time=visit
  }
  #time <- rep(0, length(visit))
  K = max(unique(Y.df$day))   # largest number of daily inputs
  J = length(unique(Y.df$visit)) # largest number of visit
  M = length(unique(Y.df$id))    # number of subjects
  N = NCOL(Y.df$Y)               # number of grid points/variables
  I = length(uni.id.visit)               # number of subject-visits
  S = nrow(Y.df)                # number of subject-visits-days
  
  num_visit = c(t(table(Y.df$id))) # number of daily visit for each subjects
  num_day = c(id.visit.vec) # number of day for each subjects_visit
  cum_visit_number = cumsum(num_visit)
  cum_day_number = cumsum(num_day)
  
  kx = tlength/N
  
  #########################################################################################
  
  # estimate grand mean mu
  mu = apply(Y, 2,mean)
  #mu = rep(0,N)
  
  # estimate eta_j
  Y.df.new = Y.df[order(Y.df$id, Y.df$visit,Y.df$day),  ]
  Y.df.new$Y.tilde = Y.df.new$Y - matrix(mu, S, N, byrow = TRUE)
  Y.tilde = Y.df.new$Y.tilde
  
  #########################################################################################
  
  
  # 4. Estimate the covariates matrix/surface Gbv and Gbd ---------------------
  ###     Estimate the three covariance functions: overall covariance G,
  ###     between visit covariance Gbv and between day covariance Gbd
  
  # covariance estimation method 1: all J_i terms inside sum i term, see Wang 2018
  if(cov.method == "m1"){
    print("not supported now")
  }
  
  #  covariance estimation method: hierarchical modeling.
  if(cov.method == "m2"){
    # 1. estimate between day covariance matrix K.W and ij-level covaraince operator
    num_day_gt1 = num_day#[which(num_day > 1)]
    i1_d <- function(i){                                            # find all j for indicator vectors for all i-(j,k) pairs
      before <- sum(id.visit.vec[1:(i-1)]) * (i > 1)
      rep(before + (1:id.visit.vec[i]), each = id.visit.vec[i])
    }
    i2_d <- function(i){                                            # find all k for indicator vectors for all i-(j,k) pairs
      before <- sum(id.visit.vec[1:(i-1)]) * (i > 1)
      rep(before + (1:id.visit.vec[i]), id.visit.vec[i])
    }
    
    ind1_d <- as.vector(unlist(sapply(1:I, i1_d)))                    # indicator vectors for j in all i-(j,k) pairs
    ind2_d <- as.vector(unlist(sapply(1:I, i2_d)))                    # indicator vectors for k in all i-(j,k) pairs
    
    Mij2 = sum(num_day_gt1^2)
    K.W = matrix(0, nrow = N, ncol = N)
    K.XU.ij = K.XU.ij1 = array(0, dim = c(length(num_day_gt1),N,N))
    num_day_sq = c(0, num_day_gt1^2)
    Gbd = matrix(0, N, N)
    for(ij in 1:length(num_day_gt1))
    {
      ij.ind = (sum(num_day_sq[1:ij])+1) : sum(num_day_sq[1:(ij+1)])
     # Y.tilde.ij = Y.tilde[unique(ind1_d[ij.ind]),]
      Y.tilde.ij = subset(Y.df.new, id.visit == uni.id.visit[ij])$Y.tilde
      Y.tilde.ij.sum = t(Y.tilde[ind1_d[ij.ind],]) %*% Y.tilde[ind2_d[ij.ind],]
      K.XU.ij[ij,,] = Y.tilde.ij.sum
      if(num_day_gt1[ij] > 1)
      {
        #K.W = K.W + (-1/num_day_gt1[ij] * Y.tilde.ij.sum+
        #               t(Y.tilde.ij) %*% Y.tilde.ij)
        K.W = K.W + (-1 * Y.tilde.ij.sum +
                       num_day_gt1[ij] * t(Y.tilde.ij) %*% Y.tilde.ij)
        ij.denom2 =  num_day_gt1[ij]^2 - num_day_gt1[ij]
        Y.tilde.ij.t1 = Y.tilde[ind1_d[ij.ind],][which(ind1_d[ij.ind] != ind2_d[ij.ind]),]
        Y.tilde.ij.t2 = Y.tilde[ind2_d[ij.ind],][which(ind1_d[ij.ind] != ind2_d[ij.ind]),]
        #K.XU.ij[ij,,] = 1/ij.denom2^2 * colSums(Y.tilde.ij.t) %*% t(colSums(Y.tilde.ij.t))
        #K.XU.ij[ij,,] = 1/(ij.denom2) * t(Y.tilde.ij.t1) %*% Y.tilde.ij.t2
        K.XU.ij[ij,,] = 1/(ij.denom2) * t(Y.tilde.ij.t1) %*% Y.tilde.ij.t2
      }
      #resid.pair.diff = pair.diff(Y.tilde.ij)
      #K.XU.ij1[ij,,] = ij.denom2 * t(Y.tilde.ij) %*% Y.tilde.ij - colSums(Y.tilde.ij) %*% t(colSums(Y.tilde.ij))
      #K.XU.ij[ij,,] = 1/num_day_gt1[ij]^2 * colSums(Y.tilde.ij) %*% t(colSums(Y.tilde.ij)) -  crossprod(resid.pair.diff)/num_day_gt1[ij]
      #K.XU.ij1[ij,,] = K.XU.ij[ij,,]
      # Y.tilde.i = subset(Y.df.new, id.visit == uni.id.visit[ij])$Y.tilde
      # Ji = nrow(Y.tilde.i)
      # if(Ji > 1){
      #   resid.pair.temp = Y.tilde.i[c(1:Ji), ]
      #   resid.pair.diff = pair.diff(resid.pair.temp)
      #   Gbd = Gbd + crossprod(resid.pair.diff) * 2
      # }
    }
    ij.denom1 = sum(num_day_gt1*(num_day_gt1-1))
    K.W = 1 / ij.denom1 * K.W
    
    # denom_1 = sum(num_day_gt1^2) - sum(num_day_gt1)
    # Gbd = 0.5/denom_1 * Gbd
    #K.W = Gbd
    
    # for(ij in 1:length(num_day_gt1))
    # {
    #   ij.denom2 =  num_day_gt1[ij]^2 - num_day_gt1[ij]
    #   #K.XU.ij[ij,,] = 1/ij.denom2 * (K.XU.ij[ij,,] - num_day_gt1[ij]/ij.denom1 * K.W)
    #   K.XU.ij[ij,,] = 1/ij.denom2 * (K.XU.ij[ij,,] - num_day_gt1[ij] * K.W)
    # }
    
    i1 <- function(i){                                            # find all j for indicator vectors for all i-(j,k) pairs
      before <- sum(J.vec[1:(i-1)]) * (i > 1)
      rep(before + (1:J.vec[i]), each = J.vec[i])
    } 
    i2 <- function(i){                                            # find all k for indicator vectors for all i-(j,k) pairs
      before <- sum(J.vec[1:(i-1)]) * (i > 1)
      rep(before + (1:J.vec[i]), J.vec[i])
    }
    
    ind1 <- as.vector(unlist(sapply(1:M, i1)))                    # indicator vectors for j in all i-(j,k) pairs
    ind2 <- as.vector(unlist(sapply(1:M, i2)))                    # indicator vectors for k in all i-(j,k) pairs
    
    K.XU = array(0, dim = c(length(ind1),N,N))
    ij.ind = which(ind1 == ind2)
    K.XU[ij.ind,,] = K.XU.ij
    for(i in 1:length(ind1))
    {
      #print(i)
      if(ind1[i] != ind2[i])
      {
        ij1j2.denom = id.visit.vec[ind1[i]] * id.visit.vec[ind2[i]]
        Y.tilde.ij1 = Y.tilde[(id.visit==names(id.visit.vec)[ind1[i]]),]
        Y.tilde.ij2 = Y.tilde[(id.visit==names(id.visit.vec)[ind2[i]]),]
        #K.XU[i,,] = 1/ij1j2.denom * matrix(colSums(kronecker(Y.tilde.ij1,Y.tilde.ij2)),N,N)
        if(ij1j2.denom==1)
          K.XU[i,,] = 1/ij1j2.denom * matrix(kronecker(Y.tilde.ij1,Y.tilde.ij2),N,N)
        else
          K.XU[i,,] = 1/ij1j2.denom * matrix(colSums(kronecker(Y.tilde.ij1,Y.tilde.ij2)),N,N)
        #K.XU[i,,] = 1/ij1j2.denom * colSums(Y.tilde.ij1) %*% t(colSums(Y.tilde.ij2))
      }
    }
    
    Y.tilde.t = Y.tilde[!duplicated(paste0(id,visit)),]
    
    time.t = time[!duplicated(paste0(id,visit))]
    X <- cbind(rep(1, mm), time.t[ind2], time.t[ind1], time.t[ind1] * time.t[ind2], (ind1 == ind2) * 1)
    beta <- solve(crossprod(X), t(X) %*% array(K.XU,c(length(ind1),N*N)))
    #beta <- ginv(crossprod(X)) %*% t(X) %*% array(K.XU,c(length(ind1),N*N))
    
    G.0 <- G.01 <- H.01 <- G.1 <- G.U <- matrix(0, N, N)
    G.0  <- array(beta[1,],dim = c(N,N))   # estimates for K.0(s,t) are in the 1st column
    #G.0 <- (G.0 + t(G.0))/2                                        # symmetry constraints yield K.0(s,t) = K.0(t,s)
    G.1 <- array(beta[4,],dim = c(N,N))  # estimates for K.1(s,t) are in the 4th column
    #G.1 <- (G.1 + t(G.1))/2                                       # symmetry constraints yield K.1(s,t) = K.1(t,s)                                   
    G.U <- array(beta[5,],dim = c(N,N))   # estimates for K.U(s,t) are in the 5th column
    #G.U <- (G.U + t(G.U))/2                                       # symmetry constraints yield K.U(s,t) = K.U(t,s)
    G.10 <- array(beta[2,],dim = c(N,N))   # estimates for K.10(s,t) are in the 2nd column
    G.01 <- array(beta[3,],dim = c(N,N))   # estimates for K.01(s,t) are in the 3rd column
    #G.01 <- (G.01 + t(G.01))/2                              # symmetry constraints yield K.01(s,t) = K.10(t,s)
    
    row.vec <- rep(1:N, each = N)     # set up row variable for bivariate smoothing                             
    col.vec <- rep(1:N, N)          # set up column variable for bivariate smoothing
    if (smoothing == TRUE & smooth.method == "sc"){      # if smoothing is selected (corresponds to method described in the paper):                       
      K.0  <- matrix(predict(gamm(as.vector(G.0)  ~ te(row.vec, col.vec, k = bf))$gam), N, N) # smooth K.0
      K.0 <- (K.0 + t(K.0)) / 2                                                            # after smoothing, symmetrize K.0 
      K.1  <- matrix(predict(gamm(as.vector(G.1)  ~ te(row.vec, col.vec, k = bf))$gam), N, N) # smooth K.1
      K.1 <- (K.1 + t(K.1)) / 2                                                            # after smoothing, symmetrize K.1
      K.01 <- matrix(predict(gamm(as.vector(G.01) ~ te(row.vec, col.vec, k = bf))$gam), N, N) # smooth K.01                                                          # after smoothing, symmetrize K.U
    }
    else {            # if no smoothing is selected (faster):
      K.U <- G.U    # do not smooth K.U (off-diagonal) 
      
      K.0 <- G.0   # do not smooth K.0
      K.01 <- (G.01 + t(G.10))/2 # do not smooth K.01
      K.1 <- G.1  # do not smooth K.1
    }
    K.X <- rbind(cbind(K.0, K.01), cbind(t(K.01), K.1))   # put together K.X, containing K.0, K.01, K.10 and K.1         
    sigma2.hat <- 0 
  }
  
  e1 = eigen(K.X, symmetric = TRUE)
  e2 = eigen(K.U, symmetric = TRUE)
  e3 = eigen(K.W, symmetric = TRUE)
  
  fpca1.value <- e1$values # eigenvalues of K.X yield estimates for the lambda_k
  fpca2.value   <- e2$values # eigenvalues of K.U yield estimates for the nu_k
  fpca3.value    <- e3$values # eigenvalues of K.W yield estimates for the xi_k
  
  ###     Keep only non-negative eigenvalues
  fpca1.value = ifelse(fpca1.value>=0, fpca1.value, 0)
  fpca2.value = ifelse(fpca2.value>=0, fpca2.value, 0)
  fpca3.value = ifelse(fpca3.value>=0, fpca3.value, 0)
  total.variance <- sum(fpca1.value) + sum(fpca2.value) + sum(fpca3.value)
  
  ###     Calculate the percentage of variance that are explained by the components
  percent1 = (fpca1.value)/sum(fpca1.value)
  percent2 = (fpca2.value)/sum(fpca2.value)
  percent3 = (fpca3.value)/sum(fpca3.value)  
  
  #########################################################################################
  
  # 6. Create level 1/2 eigenvectors/functions and eigen values -------------
  
  
  ###     Decide the number of components that are kept at level 1 and 2. The general
  ###     rule is to stop at the component where the cumulative percentage of variance
  ###     explained is greater than 90% and the variance explained by any single component
  ###     after is less than 1/N.
  if(is.na(K1) | is.na(K2) | is.na(K3))
  {
    #K1 = max( which(cumsum(percent1) < pve | percent1 > 1/N ) + 1)
    #K2 = max( which(cumsum(percent2) < pve | percent2 > 1/N ) + 1)
    #K3 = max( which(cumsum(percent3) < pve | percent3 > 1/N ) + 1)
    prop <- K1 <- K2 <- K3 <- 0
    while(prop < pve){              # add components for X or U with decreasing variance,
      # until level L of explained average variance is reached
      if (max(fpca1.value[K1 + 1],fpca2.value[K2 + 1],fpca3.value[K3 + 1]) == fpca1.value[K1 + 1]){
        K1 <- K1 + 1
      }
      else if (max(fpca1.value[K1 + 1],fpca2.value[K2 + 1],fpca3.value[K3 + 1]) == fpca2.value[K2 + 1])
      {
        K2 <- K2 + 1
      }
      else {
        K3 <- K3 + 1
      }
      prop <- (sum(fpca1.value[1:K1]) + sum(fpca2.value[1:K2]) + sum(fpca3.value[1:K3])) / total.variance  # update explained average variance
    }
  }
  ###     estimate eigen vectors
  fpca1.vectors = as.matrix(e1$vectors[, 1:K1])# * sqrt(1/kx)
  fpca2.vectors = as.matrix(e2$vectors[, 1:K2])# * sqrt(1/kx)
  fpca3.vectors = as.matrix(e3$vectors[, 1:K3])# * sqrt(1/kx)
  
  ###     Estimate amount of variability explained by between-day, 
  ###     between-visit and between-subject
  pbs = sum(fpca1.value)/total.variance
  pbv = sum(fpca2.value)/total.variance
  pbd = 1 - pbs - pbv
  #if(all(fpca3.value==0))
  #  fpca3.vectors = matrix(0, nrow = nrow(fpca3.vectors),ncol = K3)
  
  
  ###     The eigen vectosrs are unique only up to a change of signs.
  ###     Select the signs of eigenfunctions so that the integration over the domain
  ###     is non-negative
  for(i in 1:K1) {
    v1 = fpca1.vectors[,i]
    tempsign = sum(v1)
    fpca1.vectors[,i] = ifelse(tempsign<0, -1,1) * v1
  }
  for(i in 1:K2) {
    v2 = fpca2.vectors[,i]
    tempsign = sum(v2)
    fpca2.vectors[,i] = ifelse(tempsign<0, -1,1) * v2
  }
  
  for(i in 1:K3) {
    v3 = fpca3.vectors[,i]
    tempsign = sum(v3)
    fpca3.vectors[,i] = ifelse(tempsign<0, -1,1) * v3
  }  
  
  phi.0      <- fpca1.vectors[1:N, ]                                # phi.X_k = (phi.0_k, phi.1_k)
  phi.1      <- fpca1.vectors[(N+1):(2*N), ]
  
  #########################################################################################
  
  # 7. Get level 1/2 PC scores based on the projection methods.
  s1 = matrix(NA, M, K1)
  s2 = matrix(NA, I, K2)
  s3 = matrix(NA, S, K3)
  
  subject.ind <- unlist(sapply(id, function(su){which(unique(id[order(id)]) == su)}))
  subject.visit.ind <- unlist(sapply(id.visit, function(su){which(unique(id.visit[order(id,visit)]) == su)}))
  id.visit.mat <- data.frame(Y.df.new[,c("id","visit")] %>% group_by(id,visit) %>% slice(1))
  #Z.X = kronecker(Diagonal(M)[subject.ind, ], phi.0) + 
  #  kronecker(Diagonal(S, time) %*% Diagonal(M)[subject.ind, ], phi.1)  
  #Z.U = kronecker(diag(rep(1,J)),kronecker(t(rep(1,K)),fpca2.vectors))
  for(m in 1:M)
  {
    #print(m)
    index.m = which(Y.df.new$id == m)
    index.j = which(id.visit.mat$id == m)
    n_i = length(index.m)
    n_ij = length(unique(Y.df.new$visit[which(Y.df.new$id == m)]))
    B.X = kronecker(rep(1,n_i), phi.0) +
      #kronecker(Diagonal(length(index.m), time[index.m]) %*% rep(1,length(index.m)) , phi.1)
      kronecker(time[index.m], phi.1)
    U.ind = rep(1:n_ij,table(subject.visit.ind[which(subject.ind==m)]))
    B.U = kronecker(Diagonal(n_ij)[U.ind, ],fpca2.vectors)
    #B.U1 = kronecker(Diagonal(2),kronecker(rep(1,3),fpca2.vectors))
    B.W = kronecker(Diagonal(n_i), fpca3.vectors)
    B = cbind(B.X,B.U,B.W)
    xi_temp = solve(crossprod(B),t(B) %*% as.vector(c(t(Y.tilde[index.m,]))))
    s1[m,] = t(matrix(xi_temp[1:K1], K1, 1))
    s2[index.j,] = t(matrix(xi_temp[(K1+1):(K1+K2*n_ij)],  K2, n_ij))
    s3[index.m,] = t(matrix(xi_temp[(K1+K2*n_ij+1):(K1+K2*n_ij+K3*n_i)],K3,n_i))
  }
  
  
  #########################################################################################
  # 8. Create the result list as a class ------------------------------------
  
  npc = list(K1,K2,K3)
  evalues = list(fpca1.value[1:K1], fpca2.value[1:K2], fpca3.value[1:K3])
  efunctions = list(fpca1.vectors, fpca2.vectors, fpca3.vectors)
  scores = list(s1, s2, s3)
  pctvar = list(percent1[1:K1], percent2[1:K2], percent3[1:K3])
  varofTot = list(pbs, pbv, pbd)
  Y.df = Y.df.new
  
  names(efunctions) = names(evalues) = names(npc) = names(scores) = names(pctvar) = names(varofTot) = c("level1", "level2","level3")
  ret.objects = c("Y_raw","Y.df","time", "mu","npc", "pctvar","varofTot","efunctions", "evalues", "scores")
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  class(ret) = "mpca"
  return(ret)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}