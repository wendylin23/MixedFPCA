multilevel_pca_new = function(Y = NULL, id = NULL, visit=NULL, day = NULL, 
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
    
    for(i in 1:nrow(Y_sm)){
      Y_sm[i,] = gam(Y[i,] ~ s(tY,k = nk),na.action = na.omit)$fitted.values
    }
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
  Y.df$Y = Y
  
  K = max(unique(Y.df$day))   # largest number of daily inputs
  J = length(unique(Y.df$visit)) # largest number of visit
  M = length(unique(Y.df$id))    # number of subjects
  N = NCOL(Y.df$Y)               # number of grid points/variables
  I = length(uni.id.visit)               # number of subject-visits
  S = nrow(Y.df)                # number of subject-visits-days
  
  num_visit = c(t(table(Y.df$id))) # number of daily visit for each subjects
  num_day = c(t(table(id.visit))) # number of day for each subjects_visit
  cum_visit_number = cumsum(num_visit)
  cum_day_number = cumsum(num_day)
  
  kx = tlength/N
  
  #########################################################################################
  
  # estimate grand mean mu
  mu = apply(Y, 2,mean)
  
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
    Gt = matrix(0, N, N)
    for(m in 1:M){
      Y.tilde.m = subset(Y.df.new, id == uni.id[m])
      Jm = length(unique(Y.tilde.m$visit))
      Gt.temp = matrix(0, N, N)
      for(j in 1:length(Jm))
      {
        Y.tilde.mj = subset(Y.tilde.m, visit == unique(Y.tilde.m$visit)[j])$Y.tilde
        Gt.temp = Gt.temp + t(Y.tilde.mj) %*% Y.tilde.mj/nrow(Y.tilde.mj)
      }
      Gt = Gt + Gt.temp/Jm
    }
    Gt = Gt/M
    
    Gbv =  matrix(0, N, N)
    for (m in 1: M){
      Y.tilde.m = subset(Y.df.new, id == uni.id[m])
      Jm = length(unique(Y.tilde.m$visit))
      Gbv.temp1 =  matrix(0, N, N)
      if(Jm > 1){
        for(j in 1 : length(Jm))
        {
          Y.tilde.mj = subset(Y.tilde.m, visit == unique(Y.tilde.m$visit)[j])$Y.tilde
          #Gw.temp2 =  matrix(0, N, N)
          Jmk = nrow(Y.tilde.mj)
          if(Jmk >1)
          {
            for(k in 1:Jmk)
            {
              Y.tilde.nmj = subset(Y.tilde.m, visit != unique(Y.tilde.m$visit)[j] & day !=k)$Y.tilde
              Y.tilde.mj = rbind(Y.tilde.mj,Y.tilde.nmj)
              resid.pair.temp = Y.tilde.mj
              resid.pair.diff = pair.diff(resid.pair.temp)
              Gbv.temp1 = Gbv.temp1 + crossprod(resid.pair.diff)/nrow(resid.pair.diff)
              #if(nrow(Y.tilde.nmj)!=0)
              #  Gw.temp1 = Gw.temp1 + t(matrix(Y.tilde.mj[k,],nrow=nrow(Y.tilde.nmj),ncol=length(Y.tilde.mj[k,]),byrow=TRUE) )%*% Y.tilde.nmj/(Jmk*(Jmk-1))
            }
            #Gw.temp1 = Gw.temp1 + Gw.temp2
          }
        }
        Gbv = Gbv + Gbv.temp1/(Jm*(Jm-1))
      }
    }
    Gbv = Gbv/(2*M)
    
    Gbd =  matrix(0, N, N)
    for (m in 1: M){
      Y.tilde.m = subset(Y.df.new, id == uni.id[m])
      Jm = length(unique(Y.tilde.m$visit))
      Gbd.temp1 =  matrix(0, N, N)
      for(j in 1:length(Jm))
      {
        Y.tilde.mj = subset(Y.tilde.m, visit == unique(Y.tilde.m$visit)[j])$Y.tilde
        Jmk = nrow(Y.tilde.mj)
        if(Jmk >1)
        {
          resid.pair.temp = Y.tilde.mj[c(1:Jmk), ]
          resid.pair.diff = pair.diff(resid.pair.temp)
          Gbd.temp1 = Gbd.temp1 + crossprod(resid.pair.diff)/nrow(resid.pair.diff)
        }
      }
      Gbd = Gbd + Gbd.temp1/Jm
    }
    Gbd = Gbd/(2*M)
    
    Hu = Gt - Gbd
    Gbv = Hu - Gbs
  }
  
  #  covariance estimation method 2: based on MoM proposed in Koch 1968, Shou 2013.
  # 1. estimate between day covariance matrix Gbd
  if(cov.method == "m2"){
    num_day_gt1 = num_day[which(num_day > 1)]
    Gbd = matrix(0, N, N)
    for(i in 1:I){
      Y.tilde.i = subset(Y.df.new, id.visit == uni.id.visit[i])$Y.tilde
      Ji = nrow(Y.tilde.i)
      if(Ji > 1){
        # for(j1 in 1:Ji){
        #   for(j2 in setdiff(c(1:Ji),j1)){
        #     diff.j1j2 = Y.tilde.m[j1,] - Y.tilde.m[j2,]
        #     Gw = Gw + diff.j1j2 %*% t(diff.j1j2)
        #   }
        # }
        resid.pair.temp = Y.tilde.i[c(1:Ji), ]
        resid.pair.diff = pair.diff(resid.pair.temp)
        Gbd = Gbd + crossprod(resid.pair.diff) * 2
      }
    }
    
    denom_1 = sum(num_day_gt1^2) - sum(num_day_gt1)
    Gbd = 0.5/denom_1 * Gbd
    
    # 2. estimate between visit covariance matrix Gbv
    h_m <- matrix(0,N,N)
    for(m in 1:M){
      h_m1 <- matrix(0,N,N)
      Y.tilde.m <- subset(Y.df.new, id == uni.id[m])
      Jik <- nrow(Y.tilde.m)
      Ji <- length(unique(Y.tilde.m$visit))
      for(j in 1:Ji)
      { 
        temp_j <- matrix(0,N,N)
        Y.tilde.mj <- subset(Y.tilde.m, visit == unique(Y.tilde.m$visit)[j])$Y.tilde
        Kmj <- nrow(Y.tilde.mj)
        mult <- Jik-Kmj
        for (k in 1:Kmj) {
          temp_j <- temp_j + Y.tilde.mj[k,] %*% t(Y.tilde.mj[k,])
        }
        h_m1 <- h_m1 + mult * temp_j
      }
      
      g1 <- colSums(Y.tilde.m$Y.tilde) %*% t(colSums(Y.tilde.m$Y.tilde))
      g2 <- matrix(0,N,N)
      for(j in 1: Ji)
      {
        Y.tilde.mj <- subset(Y.tilde.m, visit == unique(Y.tilde.m$visit)[j])$Y.tilde
        y.ijk.sum <- colSums(Y.tilde.mj)
        g2 <- g2 + y.ijk.sum %*% t(y.ijk.sum)
      }
      gm <- g1 - g2
      h_m <- h_m + h_m1 - gm
    }
    denom_2 <- sum(num_visit^2) - sum(num_day^2)
    Hu = 2 * h_m/denom_2
    Gbv = Hu/2 - Gbd
    
    # 3. estimate between subject covariance matrix Gbs
    h1 = matrix(0,N,N)
    for(m in 1:M){
      Y.tilde.m = subset(Y.df.new, id == uni.id[m])$Y.tilde
      Ji = nrow(Y.tilde.m)
      temp = matrix(0,N,N)
      #mult = I - Ji
      mult = S - Ji
      for(j in 1:Ji){
        temp = temp + Y.tilde.m[j,] %*% t(Y.tilde.m[j,])
      }
      h1 = h1 + mult * temp
    }
    
    g1 = colSums(Y.df.new$Y.tilde) %*% t(colSums(Y.df.new$Y.tilde))
    g2 = matrix(0,N,N)
    for(m in 1:M){
      Y.tilde.m = subset(Y.df.new, id == uni.id[m])$Y.tilde
      y.ij.sum = colSums(Y.tilde.m)
      g2 = g2 + y.ij.sum %*% t(y.ij.sum)
    }
    gm = g1 - g2
    #denom_3 = I^2 - sum(num_visit^2)
    denom_3 = S^2 - sum(num_visit^2)
    Hz = 2 * h1/denom_3 - 2 * gm/denom_3
    
    Gbs = Hz/2 - Hu/2
    Gt = Hz/2
  }
  
  #########################################################################################
  
  
  # 5. Smoothing the covariance surface Gbd, Gbv and Gbs for functional dat --------
  
  ###     If "smoothing=TRUE" (which is recommended when functions are measured with noise),
  ###     this option will smooth Gbd(s,t), Gbv(s,t) and Gbs(s,t) before carrying out principal component
  ###     analysis. This step will reduce the bias in estimating eigenvalues and generate
  ###     smooth eigenfunctions.
  
  if(smoothing == TRUE & smooth.method == "sc") {
    stop("Not supported yet")
  }
  
  
  #########################################################################################
  
  # 6. Create level 1/2 eigenvectors/functions and eigen values -------------
  ###     Estimate eigen values and eigen functions at two levels
  e1 = eigen(Gbs)
  e2 = eigen(Gbv)
  e3 = eigen(Gbd)
  
  ###     Estimate amount of variability explained by between-day, 
  ###     between-visit and between-subject
  pbd = sum(diag(Gbd))/sum(diag(Gt))
  pbv = sum(diag(Gbv))/sum(diag(Gt))
  pbs = 1 - pbd - pbv
  
  
  
  ###    get eigen values
  fpca1.value = e1$values# * kx ### kx = tlength/N
  fpca2.value = e2$values# * kx
  fpca3.value = e3$values# * kx
  
  
  
  ###     Keep only non-negative eigenvalues
  fpca1.value = ifelse(fpca1.value>=0, fpca1.value, 0)
  fpca2.value = ifelse(fpca2.value>=0, fpca2.value, 0)
  fpca3.value = ifelse(fpca3.value>=0, fpca3.value, 0)
  
  ###     Calculate the percentage of variance that are explained by the components
  percent1 = (fpca1.value)/sum(fpca1.value)
  percent2 = (fpca2.value)/sum(fpca2.value)
  percent3 = (fpca3.value)/sum(fpca3.value)
  
  ###     Decide the number of components that are kept at level 1 and 2. The general
  ###     rule is to stop at the component where the cumulative percentage of variance
  ###     explained is greater than 90% and the variance explained by any single component
  ###     after is less than 1/N.
  if(is.na(K1) | is.na(K2) | is.na(K3))
  {
    K1 = max( which(cumsum(percent1) < pve | percent1 > 1/N ) + 1)
    K2 = max( which(cumsum(percent2) < pve | percent2 > 1/N ) + 1)
    K3 = max( which(cumsum(percent3) < pve | percent3 > 1/N ) + 1)
  }
  ###     estimate eigen vectors
  fpca1.vectors = e1$vectors[, 1:K1]# * sqrt(1/kx)
  fpca2.vectors = e2$vectors[, 1:K2]# * sqrt(1/kx)
  fpca3.vectors = e3$vectors[, 1:K3]# * sqrt(1/kx)
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
  
  #########################################################################################
  
  # 7. Get level 1/2 PC scores based on the projection methods.
  
  id.visit.day.dat = Y.df.new[,c(1:3)]
  id.visit.day.dat = as.data.frame(cbind(id.visit.day.dat,Y.df.new$Y.tilde))
  
  id.visit.day.all = data.frame(id = rep(uni.id, each = J*K), 
                                visit = rep(rep(c(1:J), each = K),M), 
                                day = rep(c(1:K),M*J) )
  
  dat = merge(x = id.visit.day.all, y = id.visit.day.dat, all.x = T)
  dat = dat[order(dat$id, dat$visit),]
  row.index.add = which(rowSums(is.na(dat)) == N)
  dat[row.index.add,-c(1:3)] = 0
  resid.temp2 = as.matrix(dat[,-c(1:3)]) ## create a complete residual matrix
  
  ### Create cross product of eigen vectors
  ## Matrix notation if not
  #Z.X = kronecker(Diagonal(M)[id.visit.day.all$id, ], fpca1.vectors)
  #Z.U = kronecker(Diagonal(M), kronecker(Diagonal(J)[id.visit.day.all$visit[id.visit.day.all$id==uni.id[1]],], fpca2.vectors))
  #Z.W = kronecker(Diagonal(M),kronecker(Diagonal(J*K), fpca3.vectors))
  #Z = cbind(Z.X, Z.U, Z.W)

  #Z.X = kronecker(Diagonal(1)[id.visit.day.all$id[1:(J*K)], ], fpca1.vectors)
  #Z.U = kronecker(Diagonal(1), kronecker(Diagonal(J)[id.visit.day.all$visit[id.visit.day.all$id==uni.id[1]],], fpca2.vectors))
  #Z.W = kronecker(Diagonal(1),kronecker(Diagonal(J*K), fpca3.vectors))
  #Z = cbind(Z.X, Z.U, Z.W)
  #ui = cbind(fpca)

  # int1 = matrix(0, M*J*K, K1)
  # int2 = matrix(0, M*J*K, K2)
  # int3 = matrix(0, M*J*K, K3)
  # for(i in 1:(M*J*K))   {
  #   for(j in 1:K1) int1[ i ,j] = sum( resid.temp2[i,] * fpca1.vectors[,j] ) * kx
  #   for(j in 1:K2) int2[ i ,j] = sum( resid.temp2[i,] * fpca2.vectors[,j] ) * kx
  #   for(j in 1:K3) int3[ i ,j] = sum( resid.temp2[i,] * fpca3.vectors[,j] ) * kx
  # }
  #
  # s1 = matrix(NA, M*J*K, K1)
  # s2 = matrix(NA, M*J*K, K2)
  # s3 = matrix(NA, M*J*K, K3)
  # 
  # resid.temp2.svd = svd(t(resid.temp2))
  # V_t = resid.temp2.svd$u
  # Sigma_u = diag(resid.temp2.svd$d)
  # U_t = t(resid.temp2.svd$v)
  
  # cross.integral12 = t(t(V_t)%*%fpca1.vectors)%*%(t(V_t)%*%fpca2.vectors)# * kx # Phi_1 * t(Phi_2)
  # cross.integral13 = t(t(V_t)%*%fpca1.vectors)%*%(t(V_t)%*%fpca3.vectors)# * kx # Phi_1 * t(Phi_3)
  # cross.integral23 = t(t(V_t)%*%fpca2.vectors)%*%(t(V_t)%*%fpca3.vectors)# * kx # Phi_2 * t(Phi_3)
  # cross.integral12 = t(fpca1.vectors)%*%(fpca2.vectors)# * kx # Phi_1 * t(Phi_2)
  # cross.integral13 = t(fpca1.vectors)%*%(fpca3.vectors)# * kx # Phi_1 * t(Phi_3)
  # cross.integral23 = t(fpca2.vectors)%*%(fpca3.vectors)# * kx # Phi_2 * t(Phi_3)
  # 
  # 
  # D11 = diag(rep(1,K1))*J*K
  # D12 = kronecker(t(rep(1,J)), cross.integral12 * K)
  # D13 = kronecker(t(rep(1,J*K)),cross.integral13)
  # D1 = cbind(D11, D12, D13)
  # #D22 = kronecker(diag(rep(1,J)), diag(rep(1,K2))*J*K) 
  # D22 = kronecker(diag(rep(1,J)), diag(rep(1,K2))*K) 
  # D23 = kronecker(diag(rep(1,J)),kronecker(t(rep(1,K)),cross.integral23))
  # D2 = cbind(t(D12), D22, D23)
  # D33 = kronecker(diag(rep(1,J*K)),diag(rep(1,K3)))
  # D3 = cbind(t(D13), t(D23), D33)
  # D = rbind(D1, D2, D3)
  
  # for(m in 1:M)
  # {
  #   index.m = ( (m-1) * J * K  + 1 ) : (m*J*K)
  #   # Y_hat1 = t(t(V_t)%*%fpca1.vectors) %*% Sigma_u %*% t(U_t[index.m,]) %*% rep(1,J*K)
  #   # Y_hat2 = matrix(t(t(V_t)%*%fpca2.vectors) %*% Sigma_u %*% t(U_t[index.m,]) %*% kronecker(rep(1,K),diag(rep(1,J))),K2*J, 1, byrow = TRUE)
  #   # Y_hat3 = matrix(t(t(V_t)%*%fpca3.vectors) %*% Sigma_u %*% t(U_t[index.m,]), K3 * J *K, 1, byrow = TRUE)
  #   Y_hat1 = t(fpca1.vectors) %*% t(resid.temp2[index.m,]) %*% rep(1,J*K)
  #   Y_hat2 = matrix(t(fpca2.vectors) %*% t(resid.temp2[index.m,]) %*% kronecker(rep(1,K),diag(rep(1,J))),K2*J, 1, byrow = TRUE)
  #   Y_hat3 = matrix(t(fpca3.vectors) %*% t(resid.temp2[index.m,]) , K3 * J *K, 1, byrow = TRUE)
  #   xi.temp = ginv(D) %*% rbind(Y_hat1, Y_hat2, Y_hat3)
  #   s1[index.m,] = matrix(rep(xi.temp[1:K1], each=J*K), nrow=J*K)
  #   #s2[index.m,] = matrix(rep(xi.temp[(K1+1):(K1+K2*J)], each=K), nrow=J*K)
  #   s2[index.m,] = apply(matrix(xi.temp[(K1+1):(K1+K2*J)], nrow=J,byrow = TRUE),2,function(x) rep.row(x,K))
  #   s3[index.m,] = matrix(xi.temp[(K1+K2*J+1):(K1+K2*J+K3*J*K)],nrow = J*K)
  # }
  
  s1 = matrix(NA, M, K1)
  s2 = matrix(NA, I, K2)
  s3 = matrix(NA, S, K3)
  subject.ind <- unlist(sapply(id, function(su){which(unique(id[order(id)]) == su)}))
  #subject.visit.ind <- unlist(sapply(id.visit, function(su){which(unique(id.visit[order(id.visit)]) == su)}))
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
    B.X = kronecker(rep(1,n_i), fpca1.vectors)# +
      #kronecker(Diagonal(length(index.m), time[index.m]) %*% rep(1,length(index.m)) , phi.1)
      #kronecker(time[index.m], phi.1)
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
  
  # s1 = as.data.frame(cbind(dat[,1:3],s1)) %>% group_by(id) %>% slice(1)
  # s1$visit = NULL
  # s1$day = NULL
  # names(s1)[2:ncol(s1)] = paste0("lv1_",names(s1)[2:ncol(s1)])
  # 
  # s2 = as.data.frame(cbind(dat[,1:3],s2)) %>% group_by(id, visit) %>% slice(1)
  # s2$day = NULL
  # names(s2)[3:ncol(s2)] = paste0("lv2_",names(s2)[3:ncol(s2)])
  # 
  # s3 = as.data.frame(cbind(dat[,1:3],s3))
  # s3 = merge(x = Y.df.new[,c(1:3)], y = s3, all.x = T)
  # names(s3)[4:ncol(s3)] = paste0("lv3_",names(s3)[4:ncol(s3)])
  
  
  #########################################################################################
  # 8. Create the result list as a class ------------------------------------
  
  npc = list(K1,K2,K3)
  evalues = list(fpca1.value[1:K1], fpca2.value[1:K2], fpca3.value[1:K3])
  efunctions = list(fpca1.vectors[,1:K1], fpca2.vectors[,1:K2], fpca3.vectors[,1:K3])
  scores = list(s1, s2, s3)
  pctvar = list(percent1[1:K1], percent2[1:K2], percent3[1:K3])
  varofTot = list(pbs, pbv, pbd)
  Y.df = Y.df.new
  
  names(efunctions) = names(evalues) = names(npc) = names(scores) = names(pctvar) = names(varofTot) = c("level1", "level2","level3")
  ret.objects = c("Y_raw","Y.df", "mu","npc", "pctvar","varofTot","efunctions", "evalues", "scores")
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  class(ret) = "mpca"
  return(ret)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}