rm(list=ls())
library(dplyr)
library(refund)
library(ggplot2)
library(data.table)
library(hms)
library(reshape2)
library(mgcv)
library(ICSNP)
library(SemiPar)
library(MASS)
library(fpca)
library(fda)
library(Matrix)
library(orthopolynom)

setwd("~/Dropbox/UCSD/2019Spring/Rotation")
source("R/basic/multilevel_lfpca_2stage.R")
source("R/basic/multilevel_lfpca_new.R")
source("R/basic/multilevel_pca_2level.R")
source("R/basic/multilevel_pca_3level.R")
source("R/basic/fpca.sc.R")
source("R/basic/quadWeights.R")
source("R/MixedFPCA_simulation/utils.R")

###########################################################
########## Simulations of LFPCA (2-level) #################
###########################################################

###### Simulation parameter setting ######
M = 100 ## Number of subjects
D = 600 ## Number of grid points
Nsim = 100
Ngrid = seq(0,D,length.out = D)
N.U = N.V = 4 ## number of PC numbers
balanced = "ba"
J_vec = c(3,5,7) ## Number of visits
J = J_vec[1]

input.data = expand.grid(1:M,1:J)
names(input.data) = c("id","visit")
input.data = input.data[order(input.data$id,input.data$visit),]
n.obs = dim(input.data)[1]
ids = input.data$id                                           
visit = input.data$visit
time = (visit-mean(visit)) / sqrt(var(visit))
#time = rep(0, n.obs)

## Set up basis for each level
phiV.matrix <-                matrix(NA, N.V, D)
phi1.matrix <- phi0.matrix <- matrix(NA, N.U, D)
for(k in 1:N.U)
{
  phi0.matrix[k, ] = f1.basis(Ngrid, k)
  phi1.matrix[k, ] = f1.basis(Ngrid, k+4)
}
phiU.matrix <- cbind(phi0.matrix, phi1.matrix)

for(k in 1:N.V)
{
  #phiV.matrix[k, ] = p.basis(Ngrid,k)
  phiV.matrix[k, ] = f2.basis(Ngrid,k)
}

basis_com = "fff"

## Normalize basis
normU = sqrt(diag(phiU.matrix%*%t(phiU.matrix)))
phiU.matrix = phiU.matrix/normU
phi0.matrix = phiU.matrix[,1:D]
phi1.matrix = phiU.matrix[,(D+1):(2*D)]
normV = sqrt(diag(phiV.matrix%*%t(phiV.matrix)))
phiV.matrix = phiV.matrix/normV

## Set up egien functions
eigenv = c("equal","lv1","lv2")
eigenv.opt1 = 0.5^(0:3)
eigenv.opt2 = c(4, 0.5^(1:3))
eigenv.balanced = eigenv[1]

if(eigenv.balanced=="equal")
{
  lambda1 = eigenv.opt1
  lambda2 = eigenv.opt1
}else if(eigenv.balanced=="lv1")
{
  lambda1 = eigenv.opt1
  lambda2 = eigenv.opt2
}else if(eigenv.balanced=="lv2")
{
  lambda1 = eigenv.opt2
  lambda2 = eigenv.opt1
}

N1 = length(lambda1)## Number of basis for each level
N2 = length(lambda2)
sigma_x = 0.05  

## adding error terms
nbasis = 13
err.sd = c(0.5,0.3,0)
err.add = FALSE
err.name = "no"

######### 2 level model simulation #########
## recovered signal saving
Ui     = Ui.hat  = Ui.hat.lv1 = array(NA, dim = c(M, N.U, Nsim))                  # set up empty vectors and matrices for all
Vij    = Vij.hat = array(NA, dim = c(n.obs, N.V, Nsim))  #    for all quantities to be estimated in 
phiV.hat   <- phiV.hat.lv2   <- array(NA, dim = c(D, N.V, Nsim))                  #    simulations
phi0.hat   <- phi1.hat <- array(NA, dim = c(D, N.U, Nsim))
phi0.hat.lv2   <- phi1.hat.lv2 <- array(NA, dim = c(D, N.U, Nsim))
X.matrix <- array(NA, dim = c(n.obs, D, Nsim))
#sigma2.hat <- rep(NA, rep)
lambda1.hat <- lambda1.hat.lv2 <- matrix(NA, N.U, Nsim)
lambda2.hat <- lambda2.hat.lv2 <- matrix(NA, N.V, Nsim)

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
  
  U.0 <- Ui[, , nsim] %*% phi0.matrix
  U.0 <- U.0[rep(seq(nrow(U.0)), each = J), ]
  U.1 <- Ui[, , nsim] %*% phi1.matrix
  U.1 <- U.1[rep(seq(nrow(U.1)), each = J), ]
  
  if(err.add == TRUE)
  {
    err.name = err.sd[1]
    err.sm = err.generator(nbasis, Ngrid,n.obs,err.sd[1])
    X.matrix[,,nsim] <- U.0 + U.1 * matrix(rep(time, D), n.obs, D) + err.sm
  }else{
    V   <- Vij[, , nsim]             %*% phiV.matrix
    eps <- matrix(rnorm(D * n.obs, 0, sigma_x), n.obs, D)                        # generate measurement error
    #Y.matrix[,,r] <- X.0 + X.1 * matrix(rep(time, D), n, D) + U + W + eps       # compute observations according to model (2)
    X.matrix[,,nsim] <- U.0 + U.1 * matrix(rep(time, D), n.obs, D) + V + eps
  }
}

dat_name = paste0("lv2_",M,"_",J,"_",balanced,"_",basis_com,"_",err.name,"_",eigenv.balanced,".rdata")
sim_dat = list(X.matrix = X.matrix,ids = ids, visit = visit,
               phi0.matrix = phi0.matrix, phi1.matrix = phi1.matrix,
               phiV.matrix = phiV.matrix,
               Ui = Ui, Vij = Vij)
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
mu.hat = mu.hat.m  = mu.hat.lv1 = array(NA, dim = c(D, Nsim))
Ui.hat = Ui.hat.m  = Ui.hat.lv1 = array(NA, dim = c(M, N.U, Nsim))                  # set up empty vectors and matrices for all
Vij.hat = Vij.hat.m = array(NA, dim = c(n.obs, N.V, Nsim))  #    for all quantities to be estimated in 
phiV.hat   = phiV.hat.m = array(NA, dim = c(D, N.V, Nsim))                  #    simulations
phi0.hat   = phi0.hat.m = phi0.hat.lv1 = array(NA, dim = c(D, N.V, Nsim))
phi1.hat = array(NA, dim = c(D, N.U, Nsim))
X.matrix.est = array(NA, dim = c(n.obs, D, Nsim))
#sigma2.hat <- rep(NA, rep)
lambda1.hat = lambda1.hat.m = lambda1.hat.lv1 = matrix(NA, N.U, Nsim)
lambda2.hat = lambda2.hat.m = matrix(NA, N.V, Nsim)
exp.val = exp.val.m = exp.val.lv1 = rep(0, Nsim)

set.seed(2021)
for(nsim in 1:Nsim)
{
  if (nsim%%(Nsim / 10) == 0){
    print(paste("Simulation", nsim * 100 / Nsim, "% done"))                   # give out progress report
  }
  X = X.matrix[,,nsim]
  ret.lfpca = LFPCA(X, subject = ids,
                    #time = time,
                    T = visit,
                    N.X=N.U,N.U=N.V,
                    smoothing = TRUE,# bf=2
                    smooth.method = "sf") #LFPCA
  ret.mfpca = multilevel_pca(X, ids, visit,
                             cov.method = "m1",K1=N.U, K2=N.V) # MFPCA
  X_lv1 = data.frame(X)
  X_lv1$id = ids
  X_lv1 = as.data.frame(X_lv1 %>% group_by(id) %>% summarise_all("mean"))
  X.t = as.matrix(X_lv1[,-1])
  ret.lv1 = fpca.sc.new(X.t,npc = N.U, var = TRUE)
  
  ## Save results
  ### 1. LFPCA
  mu.hat[,nsim] = ret.lfpca$eta0
  Ui.hat[,,nsim] = ret.lfpca$xi
  Vij.hat[,,nsim] = ret.lfpca$zeta
  
  phi0.hat[,,nsim] = ret.lfpca$phi.0
  phi1.hat[,,nsim] = ret.lfpca$phi.1
  phiV.hat[,,nsim] = ret.lfpca$phi.U
  
  lambda1.hat[,nsim] = ret.lfpca$lambda
  lambda2.hat[,nsim] = ret.lfpca$nu
  
  exp.val[nsim] = ret.lfpca$exvar
  
  ### 2. MFPCA
  mu.hat.m[,nsim] = ret.mfpca$mu
  Ui.hat.m[,,nsim] = as.matrix(ret.mfpca$scores$level1[,-1])
  Vij.hat.m[,,nsim] = as.matrix(ret.mfpca$scores$level2[,-c(1,2)])
  
  phi0.hat.m[,,nsim] = ret.mfpca$efunctions$level1
  phiV.hat.m[,,nsim] = ret.mfpca$efunctions$level2
  
  lambda1.hat.m[,nsim] = ret.mfpca$evalues$level1
  lambda2.hat.m[,nsim] = ret.mfpca$evalues$level2
  exp.val.m[nsim] = sum(Reduce("+",ret.mfpca$pctvar))
  
  ### 3. FPCA
  mu.hat.lv1[,nsim] = ret.lv1$mu
  Ui.hat.lv1[,,nsim] = ret.lv1$scores
  phi0.hat.lv1[,,nsim] = ret.lv1$efunctions
  lambda1.hat.lv1[,nsim] = ret.lv1$evalues
  exp.val.lv1[nsim] = ret.lv1$expvar[N.U]
}

sim_res = list(mu.hat = mu.hat, Ui.hat =  Ui.hat, Vij.hat = Vij.hat,
               phi0.hat = phi0.hat, phi1.hat = phi1.hat,
               phiV.hat = phiV.hat, exp.val = exp.val,
               lambda1.hat = lambda1.hat, lambda2.hat = lambda2.hat)

sim_res_m = list(mu.hat = mu.hat.m, Ui.hat =  Ui.hat.m, Vij.hat = Vij.hat.m,
               phi0.hat = phi0.hat.m,
               phiV.hat = phiV.hat.m, exp.val = exp.val.m,
               lambda1.hat = lambda1.hat.m, lambda2.hat = lambda2.hat.m)
sim_res_lv1 = list(mu.hat = mu.hat.lv1, Ui.hat =  Ui.hat.lv1,
               phi0.hat = phi0.hat.lv1,
               exp.val = exp.val.lv1,
               lambda1.hat = lambda1.hat.lv1)

sim_res_all = list(sim_res = sim_res,
                   sim_res_m = sim_res_m,
                   sim_res_lv1 = sim_res_lv1)

#save(sim_res_all,file = paste0("R/MixedFPCA_simulation/sim_dat/pca_res_",dat_name))

## Score errors
### LFPCA
score.err = apply(sim_res_all$sim_res$Ui.hat - sim_dat$Ui,c(1,3), function(x) x/sqrt(lambda1))
dim(score.err) = c(4, M*Nsim)
boxplot(t(score.err),pch=19,cex=0.3,col="white", outline = FALSE, ylim=c(-2,2),
        xlab = "Subject-level scores", ylab = "Normalized biases")
abline(h=0, col="red")

### MFPCA
score.m.err = apply(sim_res_all$sim_res_m$Ui.hat - sim_dat$Ui,c(1,3), function(x) x/sqrt(lambda1))
dim(score.m.err) = c(4, M*Nsim)
boxplot(t(score.m.err),pch=19,cex=0.3,col="white", outline = FALSE, ylim=c(-2,2),
        xlab = "Subject-level scores", ylab = "Normalized biases")
abline(h=0, col="red")

### FPCA
score.lv1.err = apply(sim_res_all$sim_res_lv1$Ui.hat - sim_dat$Ui,c(1,3), function(x) x/sqrt(lambda1))
dim(score.lv1.err) = c(4, M*Nsim)
boxplot(t(score.lv1.err),pch=19,cex=0.3,col="white", outline = FALSE, ylim=c(-2,2),
        xlab = "Subject-level scores", ylab = "Normalized biases")
abline(h=0, col="red")

###########################################################
###############   Regression Analysis   ###################
###########################################################

## Load data (if needed)
basis_com = "f0f"
dat_name = paste0("lv2_",M,"_",J,"_",balanced,"_",basis_com,"_",err.name,"_",eigenv.balanced,".rdata")
load(paste0("R/MixedFPCA_simulation/sim_dat_lv2/sim_dat_",dat_name))
load(paste0("R/MixedFPCA_simulation/sim_dat_lv2/pca_res_",dat_name))

X.matrix = sim_dat$X.matrix
ids = sim_dat$ids
visit = sim_dat$visit
phi0.matrix = sim_dat$phi0.matrix
phi1.matrix = sim_dat$phi1.matrix
phiV.matrix = sim_dat$phiV.matrix
Ui = sim_dat$Ui
Vij = sim_dat$Vij
tps = 1:D
Tps = matrix(tps,M,D,byrow=T)

b.option = c("gamma","linear","constant")
b.type = b.option[1]
y.outcome = matrix(0, nrow = M, ncol = Nsim) # y_i = \int b(s)U_i(s)ds + \epsilon
var.epsilon = 0.05
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
  
  y.outcome[,nsim] = Ui[, , nsim] %*% phi0.matrix %*% btrue + epsilon
  y.t = y.outcome[,nsim]
  
  ## LFPCA prediction
  U.est = sim_res_all$sim_res$Ui.hat[,,nsim] %*% t(sim_res_all$sim_res$phi0.hat[,,nsim])
  gam.fit = gam(y.t ~ s(Tps,by = U.est, bs="ps", m=1))
  mse.mat[nsim] = mean((y.t - predict(gam.fit))^2)
  
  ## MFPCA prediction
  U.est = sim_res_all$sim_res_m$Ui.hat[,,nsim] %*% t(sim_res_all$sim_res_m$phi0.hat[,,nsim])
  gam.fit = gam(y.t ~ s(Tps,by = U.est, bs="ps", m=1))
  mse.m.mat[nsim] = mean((y.t - predict(gam.fit))^2)
  
  ## FPCA prediction
  U.est = sim_res_all$sim_res_lv1$Ui.hat[,,nsim] %*% t(sim_res_all$sim_res_lv1$phi0.hat[,,nsim])
  gam.fit = gam(y.t ~ s(Tps,by = U.est, bs="ps", m=1))
  mse.lv1.mat[nsim] = mean((y.t - predict(gam.fit))^2)
}

regress_err = data.frame(LFPCA = mse.mat, MFPCA = mse.m.mat, FPCA = mse.lv1.mat)
boxplot(regress_err)

regress.res = list(y.outcome = y.outcome, regress_err = regress_err)
save(regress.res, file = paste0("R/MixedFPCA_simulation/sim_dat_lv2/regress_res_",b.type,"_",dat_name))

