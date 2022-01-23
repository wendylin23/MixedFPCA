rm(list=ls())
library(dplyr)
#library(refund)
library(ggplot2)
library(data.table)
library(hms)
library(reshape2)
library(mgcv)
library(ICSNP)
#library(SemiPar)
library(MASS)
#library(fpca)
library(fda)
library(Matrix)
library(orthopolynom)

setwd("~/Dropbox/UCSD/2019Spring/Rotation")
source("R/basic/multilevel_lfpca_2stage.R")
source("R/basic/multilevel_lfpca_new.R")
#source("R/basic/multilevel_pca_2level.R")
source("R/basic/multilevel_pca_3level.R")
#source("R/basic/fpca.sc.R")
#source("R/basic/quadWeights.R")
source("R/MixedFPCA_simulation/utils.R")

###########################################################
########## Simulations of LFPCA (2-level) #################
###########################################################

###### Simulation parameter setting ######
M = 100 ## Number of subjects
J = 3 ## Number of visits
D = 600 ## Number of grid points
Nsim = 100
Ngrid = seq(0,D,length.out = D)
N.U = N.V = N.W = 4 ## number of PC numbers
balanced = "ba"
S_vec = c(3,5,7) ## Number of days
S = S_vec[1]

input.data = expand.grid(1:M,1:J, 1:S)
names(input.data) = c("id","visit","day")
input.data = input.data[order(input.data$id,input.data$visit),]
id.visit <- paste0(input.data$id,"_",input.data$visit)
id.visit.vec = sapply(unique(id.visit), function(subj){sum(id.visit == subj)})
n.obs = dim(input.data)[1]
ids = input.data$id                                           
visit = input.data$visit
time = (visit-mean(visit)) / sqrt(var(visit))

## two level indicator
ids.t = ids[!duplicated(paste0(ids,visit))]
visit.t = visit[!duplicated(paste0(ids,visit))]
time.t = time[!duplicated(paste0(ids,visit))]
id.visit.t = rep(1:length(unique(paste0(ids,visit))),each=S)

## Set up basis for each level
phiW.matrix <-                matrix(NA, N.W, D)
phiV.matrix <-                matrix(NA, N.V, D)
phi1.matrix <- phi0.matrix <- matrix(NA, N.U, D)
for(k in 1:N.U)
{
  #phi0.matrix[k, ] = p.basis(Ngrid, 2*N.U, k)
  #phi1.matrix[k, ] = p.basis(Ngrid, 2*N.U, k+4)
  phi0.matrix[k, ] = f1.basis(Ngrid, k)
  phi1.matrix[k, ] = f1.basis(Ngrid, k+4)
}

for(k in 1:N.V)
{
  #phiV.matrix[k, ] = p.basis(Ngrid,2*k)
  phiV.matrix[k, ] = f2.basis(Ngrid,k)
}

# for(k in 1:N.U)
# {
#   #phi0.matrix[k, ] = f1.basis(Ngrid, k)
#   if(k==1)
#     phi1.matrix[k, ] = phiV.matrix[k, ] * sqrt(4)
#   else
#     phi1.matrix[k, ] = phi0.matrix[k-1, ]* sqrt(4/3)
# }

for(k in 1:N.W)
{
  #phiW.matrix[k, ] = p.basis(Ngrid,2*k-1)
  phiW.matrix[k, ] = f2.basis(Ngrid,k+4)
}

phiU.matrix <- cbind(phi0.matrix, phi1.matrix)

basis_com = "f0ff"

## Normalize basis
normU = sqrt(diag(phiU.matrix%*%t(phiU.matrix)))
phiU.matrix = phiU.matrix/normU
phi0.matrix = phiU.matrix[,1:D]
phi1.matrix = phiU.matrix[,(D+1):(2*D)]
normV = sqrt(diag(phiV.matrix%*%t(phiV.matrix)))
phiV.matrix = phiV.matrix/normV
normW = sqrt(diag(phiW.matrix%*%t(phiW.matrix)))
phiW.matrix = phiW.matrix/normW

## Set up egien functions
eigenv = c("equal","lv2","lv3")
eigenv.opt1 = 0.5^(0:3)
eigenv.opt2 = c(4, 0.5^(1:3))
eigenv.balanced = eigenv[1]

if(eigenv.balanced=="equal")
{
  lambda1 = eigenv.opt1
  lambda2 = eigenv.opt1
  lambda3 = eigenv.opt1
}else if(eigenv.balanced=="lv2")
{
  lambda1 = eigenv.opt1
  lambda2 = eigenv.opt2
  lambda3 = eigenv.opt1
}else if(eigenv.balanced=="lv3")
{
  lambda1 = eigenv.opt1
  lambda2 = eigenv.opt1
  lambda3 = eigenv.opt2
}

N1 = length(lambda1)## Number of basis for each level
N2 = length(lambda2)
N3 = length(lambda3)
sigma_x = 0.05  

## adding error terms
nbasis = 13
err.sd = c(0.3,0.1,0)
err.add = FALSE
err.name = "no"

dat_name = paste0("lv3_",M,"_",J,"_",S,"_",balanced,"_",basis_com,"_",err.name,"_",eigenv.balanced,".rdata")

######### 2 level model simulation #########
## recovered signal saving
Ui     = Ui.hat  = Ui.hat.lv1 = array(NA, dim = c(M, N.U, Nsim))                  # set up empty vectors and matrices for all
Vij    = Vij.hat = array(NA, dim = c(M*J, N.V, Nsim))                             # for all quantities to be estimated in 
Wijk   = Wijk.hat = array(NA, dim = c(n.obs, N.W, Nsim))                             # simulations
phiW.hat   <- phiW.hat.lv2   <- array(NA, dim = c(D, N.W, Nsim))
phiV.hat   <- phiV.hat.lv2   <- array(NA, dim = c(D, N.V, Nsim))                  
phi0.hat   <- phi1.hat <- array(NA, dim = c(D, N.U, Nsim))
phi0.hat.lv2   <- phi1.hat.lv2 <- array(NA, dim = c(D, N.U, Nsim))
X.matrix <- array(NA, dim = c(n.obs, D, Nsim))
#sigma2.hat <- rep(NA, rep)
lambda1.hat <- lambda1.hat.lv2 <- matrix(NA, N.U, Nsim)
lambda2.hat <- lambda2.hat.lv2 <- matrix(NA, N.V, Nsim)
lambda3.hat <- lambda3.hat.lv2 <- matrix(NA, N.W, Nsim)

set.seed(2021)
for(nsim in 1:Nsim)
{
  if (nsim%%(Nsim / 10) == 0){
    print(paste("Simulation", nsim * 100 / Nsim, "% done"))                   # give out progress report
  }
  
  for(n1 in 1:N.U)
  {
    Ui[,n1,nsim] = rnorm(M,0,sqrt(lambda1[n1]))
  }
  
  for(n2 in 1:N.V)
  {
    Vij[,n2,nsim] = rnorm(M*J,0,sqrt(lambda2[n2]))
  }
  
  for(n3 in 1:N.W)
  {
    Wijk[,n3,nsim] = rnorm(n.obs,0,sqrt(lambda3[n3]))
  }
  
  U.0 <- Ui[, , nsim] %*% phi0.matrix
  U.0 <- U.0[rep(seq(nrow(U.0)), each = J*S), ]
  U.1 <- Ui[, , nsim] %*% phi1.matrix
  U.1 <- U.1[rep(seq(nrow(U.1)), each = J*S), ]
  V   <- Vij[, , nsim]             %*% phiV.matrix
  V   <- V[rep(seq(nrow(V)), each = S),]
  
  if(err.add == TRUE)
  {
    err.name = err.sd[1]
    err.sm = err.generator(nbasis, Ngrid, n.obs,err.sd[1])
    X.matrix[,,nsim] <- U.0 + U.1 * matrix(rep(time, D), n.obs, D) + V + err.sm
  }else{
    W   <- Wijk[, , nsim]            %*% phiW.matrix 
    eps <- matrix(rnorm(D * n.obs, 0, sigma_x), n.obs, D)                        # generate measurement error
    #Y.matrix[,,r] <- X.0 + X.1 * matrix(rep(time, D), n, D) + U + W + eps       # compute observations according to model (2)
    X.matrix[,,nsim] <- U.0 + U.1 * matrix(rep(time, D), n.obs, D) + V + W + eps
  }
}
sim_dat = list(X.matrix = X.matrix,ids = ids, visit = visit,
               phi0.matrix = phi0.matrix, phi1.matrix = phi1.matrix,
               phiV.matrix = phiV.matrix, phiW.matrix = phiW.matrix,
               Ui = Ui, Vij = Vij,Wijk = Wijk)
#save(sim_dat,file = paste0("R/MixedFPCA_simulation/sim_dat/sim_dat_",dat_name))

###########################################################
####### Analysis with FPCA, MFPCA and LFPCA   #############
###########################################################

## Load data (if needed)
#load(paste0("R/MixedFPCA_simulation/sim_dat/",sim_dat_name))
# X.matrix = sim_dat$X.matrix
# ids = sim_dat$ids
# visit = sim_dat$visit
# phi0.matrix = sim_dat$phi0.matrix
# phi1.matrix = sim_dat$phi1.matrix
# phiV.matrix = sim_dat$phiV.matrix
# Ui = sim_dat$Ui
# Vij = sim_dat$Vij

## saving estimation results
mu.hat = mu.hat.m  = mu.hat.lv2 = array(NA, dim = c(D, Nsim))
Ui.hat = Ui.hat.m  = Ui.hat.lv2 = array(NA, dim = c(M, N.U, Nsim))                  # set up empty vectors and matrices for all
Vij.hat = Vij.hat.m = Vij.hat.lv2 = array(NA, dim = c(M*J, N.V, Nsim))  #    for all quantities to be estimated in 
Wijk.hat = Wijk.hat.m = array(NA, dim = c(n.obs, N.V, Nsim))
phiW.hat   = phiW.hat.m = array(NA, dim = c(D, N.W, Nsim)) 
phiV.hat   = phiV.hat.m = phiV.hat.lv2 = array(NA, dim = c(D, N.V, Nsim))                  #    simulations
phi0.hat   = phi0.hat.m = phi0.hat.lv2 = array(NA, dim = c(D, N.U, Nsim))
phi1.hat = phi1.hat.lv2 = array(NA, dim = c(D, N.U, Nsim))
X.matrix.est = array(NA, dim = c(n.obs, D, Nsim))
#sigma2.hat <- rep(NA, rep)
lambda1.hat = lambda1.hat.m = lambda1.hat.lv2 = matrix(NA, N.U, Nsim)
lambda2.hat = lambda2.hat.m = lambda2.hat.lv2 = matrix(NA, N.V, Nsim)
lambda3.hat = lambda3.hat.m = matrix(NA, N.W, Nsim)
exp.val = exp.val.m = exp.val.lv2 = rep(0, Nsim)

missing.prop = 0.8
set.seed(2021)
for(nsim in 1:Nsim)
{
  if (nsim%%(Nsim / 10) == 0){
    print(paste("Simulation", nsim * 100 / Nsim, "% done"))                   # give out progress report
  }
  X = X.matrix[,,nsim]
  id.m3 = sample(1:M, size = 0.2 * M )
  id.m2 = sample(1:M, size = 0.2 * M)
  ind.m3 = c((9*id.m3-2),(9*id.m3-1),9*id.m3)
  ind.m2 = c((9*id.m2-5),(9*id.m2-4),(9*id.m2-3))
  X.m = X[-c(ind.m3,ind.m2),]
  ids.m = ids[-c(ind.m3,ind.m2)]
  visit.m = visit[-c(ind.m3,ind.m2)]
  ret.lfpca = LFPCA_3_level(X.m, ids.m, visit.m,
                            K1=N.U, K2=N.V, K3=N.W,
                            cov.method = "m2",
                            smoothing = TRUE,# bf=2
                            smooth.method = "sf") #LFPCA
  ret.mfpca = multilevel_pca_new(X.m,ids.m,visit.m,
                                 K1=N.U, K2=N.V, K3=N.W,
                                 cov.method = "m2",
                                 smoothing = TRUE, smooth.method = "sf"
                                 )
  ## two level indicator
  ids.t = ids.m[!duplicated(paste0(ids.m,visit.m))]
  visit.t = visit.m[!duplicated(paste0(ids.m,visit.m))]
  #time.t = time[!duplicated(paste0(ids,visit))]
  id.visit.t = rep(1:length(unique(paste0(ids.m,visit.m))),each=S)
  X_lv2 = as.matrix(aggregate(X.m, list(id.visit.t), mean)[,-1])
  ret.lv2 = LFPCA(Y=X_lv2, subject = ids.t, T = visit.t,
                  N.X=N.U,N.U=N.V,
                  smoothing = TRUE,# bf=2
                  smooth.method = "sf")
  
  ## Save results
  ### 1. LFPCA
  mu.hat[,nsim] = ret.lfpca$mu
  Ui.hat[,,nsim] = as.matrix(ret.lfpca$scores$level1)
  Vij.hat[,,nsim] = as.matrix(ret.lfpca$scores$level2)
  Wijk.hat[,,nsim] = as.matrix(ret.lfpca$scores$level3)
  
  phi0.hat[,,nsim] = ret.lfpca$efunctions$level1[1:D,]
  phi1.hat[,,nsim] = ret.lfpca$efunctions$level1[(D+1):(2*D),]
  phiV.hat[,,nsim] = ret.lfpca$efunctions$level2
  phiW.hat[,,nsim] = ret.lfpca$efunctions$level3
  
  lambda1.hat[,nsim] = ret.lfpca$evalues$level1
  lambda2.hat[,nsim] = ret.lfpca$evalues$level2
  lambda3.hat[,nsim] = ret.lfpca$evalues$level3
  
  #exp.val[nsim] = ret.lfpca
  
  ### 2. MFPCA
  mu.hat.m[,nsim] = ret.mfpca$mu
  Ui.hat.m[,,nsim] = as.matrix(ret.mfpca$scores$level1)
  Vij.hat.m[,,nsim] = as.matrix(ret.mfpca$scores$level2)
  Wijk.hat.m[,,nsim] = as.matrix(ret.mfpca$scores$level3)
  
  phi0.hat.m[,,nsim] = ret.mfpca$efunctions$level1
  phiV.hat.m[,,nsim] = ret.mfpca$efunctions$level2
  phiW.hat.m[,,nsim] = ret.mfpca$efunctions$level3
  
  lambda1.hat.m[,nsim] = ret.mfpca$evalues$level1
  lambda2.hat.m[,nsim] = ret.mfpca$evalues$level2
  lambda3.hat.m[,nsim] = ret.mfpca$evalues$level3
  #exp.val.m[nsim] = sum(Reduce("+",ret.mfpca$varofTot))
  
  ### 3. LFPCA-lv2
  mu.hat.lv2[,nsim] = ret.lv2$eta0
  Ui.hat.lv2[,,nsim] = ret.lv2$xi
  Vij.hat.lv2[,,nsim] = ret.lv2$zeta
  phi0.hat.lv2[,,nsim] = ret.lv2$phi.0
  phi1.hat.lv2[,,nsim] = ret.lv2$phi.1
  phiV.hat.lv2[,,nsim] = ret.lv2$phi.U
  
  lambda1.hat.lv2[,nsim] = ret.lv2$lambda
  lambda2.hat.lv2[,nsim] = ret.lv2$nu
  
  #exp.val.lv2[nsim] = ret.lv2$exvar
}

sim_res = list(mu.hat = mu.hat, Ui.hat =  Ui.hat, Vij.hat = Vij.hat, Wijk.hat = Wijk.hat,
               phi0.hat = phi0.hat, phi1.hat = phi1.hat,
               phiV.hat = phiV.hat, phiW.hat = phiW.hat,
               lambda1.hat = lambda1.hat, lambda2.hat = lambda2.hat,lambda3.hat = lambda3.hat)

sim_res_m = list(mu.hat = mu.hat.m, Ui.hat =  Ui.hat.m, Vij.hat = Vij.hat.m, Wijk.hat.m = Wijk.hat.m,
               phi0.hat = phi0.hat.m,
               phiV.hat = phiV.hat.m, phiW.hat = phiW.hat.m,
               lambda1.hat = lambda1.hat.m, lambda2.hat = lambda2.hat.m,lambda3.hat = lambda3.hat.m)
sim_res_lv2 = list(mu.hat = mu.hat.lv2, Ui.hat =  Ui.hat.lv2,Vij.hat = Vij.hat.lv2,
               phi0.hat = phi0.hat.lv2, phi1.hat = phi1.hat.lv2,phiV.hat = phiV.hat.lv2,
               lambda1.hat = lambda1.hat.lv2, lambda2.hat = lambda2.hat.lv2)

sim_res_all = list(sim_res = sim_res,
                   sim_res_m = sim_res_m,
                   sim_res_lv2 = sim_res_lv2)

#save(sim_res_all,file = paste0("R/MixedFPCA_simulation/sim_dat/pca_res_",dat_name))

## Score errors
### LFPCA
score.err = apply(sim_res_all$sim_res$Ui.hat - sim_dat$Ui,c(1,3), function(x) x/sqrt(lambda1))
dim(score.err) = c(4, M*Nsim)
boxplot(t(score.err),pch=19,cex=0.3,col="white",# outline = FALSE, ylim=c(-4,4),
        xlab = "Subject-level scores", ylab = "Normalized biases")
abline(h=0, col="red")

score.err = apply(sim_res_all$sim_res$Vij.hat - sim_dat$Vij,c(1,3), function(x) x/sqrt(lambda2))
dim(score.err) = c(4, M*J*Nsim)
boxplot(t(score.err),pch=19,cex=0.3,col="white",# outline = FALSE, ylim=c(-4,4),
        xlab = "Visit-level scores", ylab = "Normalized biases")
abline(h=0, col="red")

### MFPCA
score.m.err = apply(sim_res_all$sim_res_m$Ui.hat - sim_dat$Ui,c(1,3), function(x) x/sqrt(lambda1))
dim(score.m.err) = c(4, M*Nsim)
boxplot(t(score.m.err),pch=19,cex=0.3,col="white",# outline = FALSE, ylim=c(-4,4),
        xlab = "Subject-level scores", ylab = "Normalized biases")
abline(h=0, col="red")

score.err = apply(sim_res_all$sim_res_lv2$Vij.hat - sim_dat$Vij,c(1,3), function(x) x/sqrt(lambda2))
dim(score.err) = c(4, M*J*Nsim)
boxplot(t(score.err),pch=19,cex=0.3,col="white",# outline = FALSE, ylim=c(-4,4),
        xlab = "Visit-level scores", ylab = "Normalized biases")
abline(h=0, col="red")

### FPCA
score.lv2.err = apply(sim_res_all$sim_res_lv2$Ui.hat - sim_dat$Ui,c(1,3), function(x) x/sqrt(lambda1))
dim(score.lv2.err) = c(4, M*Nsim)
boxplot(t(score.lv2.err),pch=19,cex=0.3,col="white",# outline = FALSE, ylim=c(-4,4),
        xlab = "Subject-level scores", ylab = "Normalized biases")
abline(h=0, col="red")

score.err = apply(sim_res_all$sim_res_lv2$Vij.hat - sim_dat$Vij,c(1,3), function(x) x/sqrt(lambda2))
dim(score.err) = c(4, M*J*Nsim)
boxplot(t(score.err),pch=19,cex=0.3,col="white",# outline = FALSE, ylim=c(-4,4),
        xlab = "Visit-level scores", ylab = "Normalized biases")
abline(h=0, col="red")

###########################################################
###############   Regression Analysis   ###################
###########################################################

## Load data (if needed)
basis_com = "fff"
dat_name = paste0("lv3_",M,"_",J,"_",S,"_",balanced,"_",basis_com,"_",err.name,"_",eigenv.balanced,".rdata")
load(paste0("R/MixedFPCA_simulation/sim_dat_lv3/sim_dat_",dat_name))
load(paste0("R/MixedFPCA_simulation/sim_dat_lv3/pca_res_",dat_name))

X.matrix = sim_dat$X.matrix
ids = sim_dat$ids
visit = sim_dat$visit
phi0.matrix = sim_dat$phi0.matrix
phi1.matrix = sim_dat$phi1.matrix
phiV.matrix = sim_dat$phiV.matrix
phiW.matrix = sim_dat$phiW.matrix
Ui = sim_dat$Ui
Vij = sim_dat$Vij
Wijk = sim_dat$Wijk
tps = 1:D
Tps = matrix(tps,M*J,D,byrow=T)
subjf =  as.factor(ids.t)

b.option = c("gamma","linear","constant")
b.type = b.option[1]
y.outcome = matrix(0, nrow = M*J, ncol = Nsim) # y_i = \int b(s)U_i(s)ds + \epsilon
var.epsilon = 0.05
b0.err = sqrt(0.1)
b0 = matrix(0, nrow = M*J,ncol = Nsim)
mse.mat = mse.m.mat = mse.lv1.mat = rep(0, Nsim)

set.seed(2021)
for(nsim in 1:Nsim)
{
  if (nsim%%(Nsim / 10) == 0){
    print(paste("Simulation", nsim * 100 / Nsim, "% done"))                   # give out progress report
  }
  
  if(b.type == "gamma")
    btrue = 10*dgamma(Ngrid, 5, 1/5)
  else if(b.type == "linear")
    btrue = 5e-4 * Ngrid
  else if(b.type == "constant")
    btrue = rep(0.5, D)
  
  ## Simulated scalar outcomes
  epsilon = rnorm(M,0, var.epsilon)
  
  U.0 <- Ui[, , nsim] %*% phi0.matrix
  U.0 <- U.0[rep(seq(nrow(U.0)), each = J), ]
  U.1 <- Ui[, , nsim] %*% phi1.matrix
  U.1 <- U.1[rep(seq(nrow(U.1)), each = J), ]
  V   <- Vij[, , nsim]             %*% phiV.matrix
  b0_t <- rnorm(M,mean=0,sd=b0.err)
  b0[,nsim] <- rep(b0_t, each=J)
  #y.outcome[,nsim] = b0[,nsim] + (U.0 + U.1 * matrix(rep(time.t, D), M*J , D) + V) %*% btrue + epsilon
  y.outcome[,nsim] = b0[,nsim] + (U.0 + V) %*% btrue + epsilon
  
  id.m2 = sim_dat$id.m2[,nsim]
  id.m3 = sim_dat$id.m3[,nsim]
  ind.t = (!(ids.t%in%id.m2 & visit.t==2))&(!(ids.t%in%id.m3 & visit.t==3))
  Tps.t = Tps[ind.t,]
  subjf.t = subjf[ind.t]
  y.t = y.outcome[ind.t,nsim]
  
  ## LFPCA prediction
  U0.est = sim_res_all$sim_res$Ui.hat[,,nsim] %*% t(sim_res_all$sim_res$phi0.hat[,,nsim])
  U0.est = U0.est[rep(seq(nrow(U0.est)), each = J), ]
  U1.est = sim_res_all$sim_res$Ui.hat[,,nsim] %*% t(sim_res_all$sim_res$phi1.hat[,,nsim])
  U1.est = U1.est[rep(seq(nrow(U1.est)), each = J), ]
  U.est = U0.est + U1.est * matrix(rep(time.t, D), M*J , D)
  U.est = U.est[ind.t,]
  V.est = sim_res_all$sim_res$Vij.hat[,,nsim] %*% t(sim_res_all$sim_res$phiV.hat[,,nsim])
  gam.fit = gamm(y.t ~ s(Tps.t,by = U.est, bs="ps", m=1) + s(Tps.t,by = V.est, bs="ps", m=1), random=list(subjf.t=~1))
  mse.mat[nsim] = mean((y.t - predict(gam.fit$lme))^2)
  
  ## MFPCA prediction
  U.est = sim_res_all$sim_res_m$Ui.hat[,,nsim] %*% t(sim_res_all$sim_res_m$phi0.hat[,,nsim])
  U.est = U.est[rep(seq(nrow(U.est)), each = J), ]
  U.est = U.est[ind.t,]
  V.est = sim_res_all$sim_res_m$Vij.hat[,,nsim] %*% t(sim_res_all$sim_res_m$phiV.hat[,,nsim])
  gam.fit = gamm(y.t ~ s(Tps.t,by = U.est, bs="ps", m=1) + s(Tps.t,by = V.est, bs="ps", m=1), random=list(subjf.t=~1))
  mse.m.mat[nsim] = mean((y.t - predict(gam.fit$lme))^2)
  
  ## LFPCA-lv2 prediction
  U0.est = sim_res_all$sim_res_lv2$Ui.hat[,,nsim] %*% t(sim_res_all$sim_res_lv2$phi0.hat[,,nsim])
  U0.est = U0.est[rep(seq(nrow(U0.est)), each = J), ]
  U1.est = sim_res_all$sim_res_lv2$Ui.hat[,,nsim] %*% t(sim_res_all$sim_res_lv2$phi1.hat[,,nsim])
  U1.est = U1.est[rep(seq(nrow(U1.est)), each = J), ]
  U.est = U0.est + U1.est * matrix(rep(time.t, D), M*J , D)
  U.est = U.est[ind.t,]
  V.est = sim_res_all$sim_res_lv2$Vij.hat[,,nsim] %*% t(sim_res_all$sim_res_lv2$phiV.hat[,,nsim])
  gam.fit = gamm(y.t ~ s(Tps.t,by = U.est, bs="ps", m=1) + s(Tps.t,by = V.est, bs="ps", m=1), random=list(subjf.t=~1))
  mse.lv1.mat[nsim] = mean((y.t - predict(gam.fit$lme))^2)
}

regress_err = data.frame(LFPCA = mse.mat, MFPCA = mse.m.mat, FPCA = mse.lv1.mat)
boxplot(regress_err)

regress.res = list(y.outcome = y.outcome, regress_err = regress_err)
save(regress.res, file = paste0("R/MixedFPCA_simulation/sim_dat_lv3/regress_res_",b.type,"_",dat_name))

###########################################################
###############   Code checking plots   ###################
###########################################################

rets = ret.lfpca
id.visit <- paste0(ids.m,"_",visit.m)
id.visit.vec = sapply(unique(id.visit), function(subj){sum(id.visit == subj)})
time <- ret3$time
mu.hat <- ret3$mu
Xi.hat <- ret3$scores$level1
Uij.hat <- ret3$scores$level2
Wijk.hat <- ret3$scores$level3

phi0.hat <- ret3$efunctions$level1[1:D,]
phi1.hat <- ret3$efunctions$level1[(D+1):(2*D),]
phiU.hat <- ret3$efunctions$level2
phiW.hat <- ret3$efunctions$level3
lambda1.hat <- ret3$evalues$level1
lambda2.hat <- ret3$evalues$level2
lambda3.hat <- ret3$evalues$level3

id_counts = data.frame(table(ids.m))$Freq
rep_create1<-function(i1){
  rep(Xi.hat[i1,],id_counts[i1])
}
score1 <- matrix(unlist(sapply(1:length(unique(ids.m)),rep_create1)),ncol=dim(Xi.hat)[2],byrow=TRUE)
#visit_counts = data.frame(table(paste0(id,visit)))$Freq
visit_counts = id.visit.vec
rep_create2<-function(i1){
  rep(Uij.hat[i1,],visit_counts[i1])
}
score2 <- matrix(unlist(sapply(1:length(unique(paste0(ids.m,visit.m))),rep_create2)),ncol=dim(Uij.hat)[2],byrow=TRUE)
score3 <- Wijk.hat

n = dim(X.m)[1]
mu.mat <- matrix(mu.hat, n, D, byrow = TRUE)
W.hat   <- score3 %*% t(phiW.hat)                 # generate day-specific deviations
U.hat   <- score2 %*% t(phiU.hat)  
X.0.hat <- score1 %*% t(phi0.hat)                  # generate random functional intercept
X.1.hat <- score1 %*% t(phi1.hat)                  # generate random functional slope
Y.hat   <- mu.mat + X.0.hat + X.1.hat * matrix(rep(time, D), n, D) + U.hat + W.hat     # compute observations according to model (2)

subject1 = ids.m[which.max(score1[,1]-score2[,1])]
subject2 = ids.m[which.max(score2[,1]-score1[,1])]

Ji = length(unique(visit.m[ids.m==subject1]))
par(oma=c(2,2,0,4),mar=c(3,3,2,0),mfrow = c(Ji,1))
for(j in 1:Ji)
{
  Y.plot <- X.m[(ids.m == subject1) & (visit.m==j), ]
  Y.plot.sm <- ret3$Y.df$Y[(ids.m == subject1) & (visit.m==j),]
  Y.plot.hat <- Y.hat[(ids.m == subject1) & (visit.m==j), ]
  plot(c(1, D), c(min(Y.plot), max(Y.plot)), col = 0, xlab = "", ylab = "", main = paste("True_Subject ", subject1," Visit",j, sep = ""))
  for (s in 1: dim(Y.plot)[1]){
    lines(Y.plot[s, ], col = s,lwd = 0.5)
    lines(Y.plot.sm[s,], col = s, lwd = 1)
    lines(Y.plot.hat[s, ], col = s,lwd = 2,lty=2)
  }
  #legend(x = "bottomleft", col = c(1:S), lty = 1, legend = 1:S, title = "day")
  
  # plot(c(1, D), c(min(Y.plot.hat), max(Y.plot.hat)), col = 0, xlab = "", ylab = "", main = paste("Estimate_Subject ", i," Visit",j, sep = ""))
  # for (s in 1: dim(Y.plot.hat)[1]){
  #   lines(Y.plot.hat[s, ], type = "l", col = s)
  # }
}
legend(x = "bottomleft", col = c(1:dim(Y.plot.hat)[1]), cex=0.4, bty="n",
       lty = 1, legend = 1:dim(Y.plot.hat)[1], title = "day")
