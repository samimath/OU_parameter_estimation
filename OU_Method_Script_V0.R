#Inverse of the complete-data covariance matrix, matrix square root and related derivatives

#Explicit form of B^-1

#This function computes the exact inverse of the covariance matrices used in the OU process

OU_COVINV <- function(b){

# INPUT:
#       b = vector of elements that form the sq matrx.
#       e.g. b(1:N)=exp(?mu.*nu(1:N));
#       Assign values to border elements:

  N <- length(b)
  B_inv <- matrix(NA, nrow=N, ncol=N)
  B_inv[1,1] <- 1/(1-(b[2])^2)
  B_inv[N,N] <- 1/(1-(b[N])^2)
  for (k in 2:(N-1)){
    B_inv[k,k]=1/(1-(b[k+1])^2)+1/(1-(b[k])^2)-1
  }
  for (k in 2:N){
    B_inv[k, (k-1)]=-b[k]/(1-(b[k])^2)
    B_inv[(k-1),k]=-b[k]/(1-(b[k])^2)
  }
  B_inv[is.na(B_inv)] <- 0
  return(B_inv)
}



#Bidiagonal matrix square root of B^-1

# This function finds the lower bidiagonal
# matrix square root for a symmatric tridiagonal matrix
# INPUT:
  # G = symmetric tridiagonal matric
# OUTPUT:
  # L = lower bidiagonal matrix such that G=L*L'


OU_SQRTM <- function(G){
  #in matlab: norm(X)return 2-norm by default, G' means transposition of a matrix; 
  #in R: norm(X)returns 1-norm by default, t(G) means transposition of a matrix
  try (if (norm(G-t(G),type="2") > 0) stop("G must be symmetric and tridiagonal"))
  K <- dim(G)[1]  # K=size(G,1), number of rows (dimension 1) in Matlab 
  L <- matrix(0,nrow=K,ncol=K) # L=zeros(K,K), generate a zero matrix of K rows and K columns
  L[1,1] <- sqrt(G[1,1])
  L[K,K] <- sqrt(G[K,K])
  for (j in 2:K){
    L[j,(j-1)] <- G[(j-1),j]/L[(j-1),(j-1)]
    L[j,j] <- sqrt(G[j,j]-L[j,(j-1)]^2)
  }
  return(L)
}


#Computing DμB^?1 and D2μB^?1
# This function computes the derivative of
# the inverse of the covariance matrix B??1(nmu):
  # function DB inv = OU DBINV(mu,v)

OU_DBINV <- function(mu,v){
# INPUT:
  # mu = parameter value for B(mu)
  # v = Vector for the vertical coordinate of the random field
# OUTPUT:
  # DB inv = A matrix of derivatives for B??1(nmu)
  # Assign values to border elements:
  N <- length(v)
  DB_inv<- matrix(0,nrow=N,ncol=N)
# standaridze vectors into column vectors
  v <- matrix(v)
  
  if (dim(v)[2] != 1){  #??I'm not sure what this step does
    v <- t(v)
  }
  nu <- rbind(0,abs(v[2:N, ,drop=FALSE]-v[1:(N-1), ,drop=FALSE]))
  b <- exp(-mu*nu)
  b_sq <- b^2   # b_sq=b.^2
  DB_inv <- matrix(NA, nrow=N, ncol=N)
  DB_inv[1,1] <- -2*nu[2]*b_sq[2]/((1-b_sq[2])^2)
  DB_inv[N,N] <- -2*nu[N]*b_sq[N]/((1-b_sq[N])^2)
  for (j in 2:(N-1)){
    DB_inv[j,j] <- -2*(nu[j]*b_sq[j]*(1+b_sq[j]))^2+nu[j+1]*b_sq[j+1]/((1-b_sq[j+1])^2)
  }
  for (j in 2:N){
    DB_inv[(j-1),j] <- nu[j]*b[j]*(1+b_sq[j])/((1-b_sq[j])^2)
    DB_inv[j,(j-1)] <- DB_inv[(j-1),j]
  }
  DB_inv[is.na(DB_inv)] <- 0
  return(DB_inv)
}


# This function computes the 2nd derivtive
# of the inverse of the covariance matrix B??1(nmu):
  # function DB2 inv = OU DBINV(mu,v)

OU_D2BINV <- function(mu,v){
# INPUT:
  # mu = parameter value for B(mu)
  # v = Vector for the vertical coordinate of the random field
# OUTPUT:
  # DB inv = A matrix of derivatives for B??1(nmu)
  # Assign values to border elements:
  
  N <- length(v)
  DB2_inv <- matrix(0,N,N)
# standaridze vectors into column vectors
  v <- matrix(v)
  
  if (dim(v)[2] != 1){  #??I'm not sure what this step does
    v <- t(v)
  }
  nu <- rbind(0,abs(v[2:N, ,drop=FALSE]-v[1:(N-1), ,drop=FALSE]))
  b <- exp(-mu*nu)
  b_sq <- b^2
  One <- matrix(1,nrow=N, ncol=1)
  b_star <- (nu*b_sq*(One+b_sq))/((One-b_sq)^3)
  DB2_inv[1,1] <- 4*b_star[2]
  DB2_inv[N,N] <- 4*b_star[N]
  for (j in 2:(N-1)){
    DB2_inv[j,j] <- 4*(b_star[j]+b_star[j+1])
  }
  for (j in 2:N){
    DB2_inv[(j-1),j]<-((nu[j])^2*b[j]*(1+6*(b[j])^2+(b[j])^4))/((1-b_sq[j])^2)
    DB2_inv[j,(j-1)]<-DB2_inv[(j-1),j]
  }
  DB2_inv[is.na(DB2_inv)] <- 0
  return (DB2_inv)
}


# This function finds the lower?bidiagonal
# matrix square?root for a symmatric tridiagonal matrix
# INPUT:
  # G = symmetric tridiagonal matric
# OUTPUT:
  # L = lower bidiagonal matrix such that G=L*L'

OU_SQRTM <- function(G){

  try (if (norm(G-t(G),type="2") > 0) stop("G must be symmetric and tridiagonal"))
  K <- dim(G)[1]
  L <- matrix(0,nrow=K, ncol=K)
  L[1,1] <- sqrt(G[1,1])
  L[N,N] <- sqrt(G[K,K])
  for (j in 2:K){
    L[j,(j-1)] <- G[(j-1),j]/L[(j-1),(j-1)]
    L[j,j] <- sqrt(G[j,j]-L[j,(j-1)]^2)
  }
  return(L)
}
  
  
# Estimation schemes, likelihood functions and related derivatives
  
  # Complete-data likelihood function
  
  # OU LIKE.m
  # Function to evaluate the ?2 log likelihood function
  # of the parameter values given an observation of
  # the OU process
  # Author : Sami Cheong
  # Date : 7/29/14
  # Version : 1
  
#  function likelihood = OU LIKE(lambda,mu,sigma2,u,v,X)  
  
  OU_LIKE <- function (lambda,mu,sigma2,u,v,X){
  
    N <- dim(X)[1]
    M <- dim(X)[2]
    u <- matrix(u)
    if (dim(u)[2] != 1){
      u <- t(u)
    }
    v <- matrix(v)
    if (dim(v)[2] != 1){
      v <- t(v)
    }
    eta <- rbind(0,abs(u[2:M, ,drop=FALSE]-u[1:(M-1), ,drop=FALSE]))
    nu <- rbind(0, abs(v[2:N, ,drop=FALSE]-v[1:(N-1), ,drop=FALSE]))
    b <- exp(-mu*nu)
    
    # define the different terms in the likelihood function
    term_pi <- M*N*log(2*pi)
    term_sigma <- M*N*log(sigma2)
    term_lambda <- N*sum(log(1-exp(-2*mu*nu[2:N])))
    term_mu <- M*sum(log(1-exp(-2*mu*nu[2:N])))
    
    # initialize the terms for the quadratic form n
    long_term <- matrix(0, nrow=M, ncol=1)
    Xa <- matrix(0, nrow=N, ncol=M)  # long_term=zeros(size(X)), the same size of X
    Binv <- OU_COVINV(b)
    Xa[,1] <- X[,1]
    for (i in 2:M){
      Xa[,i] <- X[,i] - (exp(-lambda*eta[i])*X[,(i-1)])
      long_term[i] <- (t(Xa[,i])%*%Binv%*%Xa[,i])/(1-exp(-2*lambda*eta[i]))
    }
    likelihood = term_pi+term_sigma+term_lambda+term_mu+
      ((t(X[,1])%*%Binv%*%X[,1]+sum(long_term[2:M])))/sigma2
    return(likelihood)
  }
  


# Parameter updates based on the Hessian matrix of the completedata likelihood function
  # function hess = OU LIKE HESS UPDATE(lambda,mu,sigma2,u,v,X)
  # function to evaluate the hessian matrix
  # of the complete data likelihood
  # generate structure :H=[h lmblmb, h lmu; h mul, h mumu]
  # Uses other functinos : OU DBINV(mu,v) , OU D2BINV(mu,v)

OU_LIKE_HESS_UPDATE <- function(lambda, mu, sigma2, u, v, X){
  # H=[H 11, H 12; H 21, H 22]
  # define the components of the Hessian matrix:
  M <- dim(X)[2]
  N <- dim(X)[1]
  
  # standardize things to be column vectors :
  u <- matrix(u)
  if (dim(u)[2] != 1){
    u <- t(u)
  }
  
  v <- matrix(v)
  if (dim(v)[2] != 1){
    v <- t(v)
  }
  
  eta <- rbind(0,abs(u[2:length(u), ,drop=FALSE]-u[1:(length(u)-1), ,drop=FALSE]))
  zeta <- rbind(0,abs(v[2:length(v), ,drop=FALSE]-v[1:(length(v)-1), ,drop=FALSE]))
  zeta_expmu <- exp(-mu*zeta)
  zeta_expmu2 <- 1 - zeta_expmu^2
  eta_explmb <- exp(-lambda*eta)
  eta_explmb2 <- 1 - eta_explmb^2
  
  # compute B?1:
  B_inv <- OU_COVINV(zeta_expmu)
  B_inv[is.na(B_inv)] <- 0
  # initialize stuff
  # X q=NaN(size(X));
  X_q <- matrix(0,N,M)
  h_lmblmb_1 <- numeric(M)
  h_lmblmb_2 <- numeric(M)
  h_lmblmb_3 <- numeric(M)
  h_lmblmb_4 <- numeric(M)
  h_lmblmb_5 <- numeric(M)
  h_lmbmu_1 <- numeric(M)
  h_lmbmu_2 <- numeric(M)
  h_mulmb_1 <- numeric(M)
  h_mulmb_2 <- numeric(M)
  h_mumu_star <- numeric(M)
  dl_mu_star <- numeric(M)
  dl_lmb_1 <- numeric(M)
  dl_lmb_2 <- numeric(M)
  dl_lmb_3 <- numeric(M)
  for (j in 2:M){
    X_q[,j]<-X[,j,drop=FALSE] - eta_explmb[j]*X[,(j-1),drop=FALSE]
    h_lmblmb_1[j] <- -4*(eta[j]*eta_explmb[j])^2/(eta_explmb2[j]^2)
    h_lmblmb_2[j] <- (eta[j]^2*eta_explmb[j])*(1+(eta_explmb[j])^2)*
      (t(X[,(j-1),drop=FALSE]) %*% B_inv %*% X_q[,j,drop=FALSE])/((eta_explmb2[j])^2)
    h_lmblmb_3[j] <- (eta[j]*eta_explmb[j])^2*(t(X[,(j-1)])) %*% B_inv %*% X[,(j-1)]/eta_explmb2[j]
    h_lmblmb_4[j] <- 2*(eta[j]*eta_explmb[j])^2*(1+(eta_explmb[j])^2)*
      t(X_q[,j]) %*% B_inv %*% X_q[,j]/(eta_explmb2[j]^3)
    h_lmblmb_5[j] <- 2*(eta[j]*eta_explmb[j])^2*eta_explmb[j]*
      t(X[,(j-1)]) %*% B_inv %*% X_q[,j]/(eta_explmb2[j]^2)
    
    # split the H 12 = h_lambdamu:
    h_lmbmu_1[j] <- eta[j]*eta_explmb[j]*t(X[,(j-1)])%*%
      OU_DBINV(mu,v) %*% X_q[,j]/eta_explmb2[j]
    h_lmbmu_2[j] <- eta[j]*eta_explmb[j]^2*t(X_q[,j])%*%
      OU_DBINV(mu,v) %*% X_q[,j]/(eta_explmb2[j]^2)
    h_mulmb_1[j] <- eta[j]*(eta_explmb[j]^2)*t(X[,(j-1)])%*%
      OU_DBINV(mu,v) %*% X_q[,j]/eta_explmb2[j]
    h_mulmb_2[j] <- eta[j]*(eta_explmb[j]^2)*t(X_q[,j])%*%
      OU_DBINV(mu,v) %*% X_q[,j]
    
    # define the lambda term in h mumu:
    h_mumu_star[j] <- (t(X_q[,j]) %*% OU_D2BINV(mu,v) %*% X_q[,j])/eta_explmb2[j]
    dl_mu_star[j] <- (t(X_q[,j]) %*% OU_DBINV(mu,v) %*% X_q[,j])/eta_explmb2[j]
    
    # define the terms in the first derivative of l(theta):
    dl_lmb_1[j] <- 2*eta[j]*(eta_explmb[j]^2)/eta_explmb2[j]
    dl_lmb_2[j] <- eta[j]*eta_explmb[j]*t(X[,(j-1)]) %*% B_inv %*% X_q[,j]/eta_explmb2[j]
    dl_lmb_3[j] <- eta[j]*(eta_explmb[j]/eta_explmb2[j])^2*t(X[,(j-1)]) %*% B_inv %*% X_q[,j]
  }
  h_mumu_1 <- numeric(N)
  dl_mu_1 <- numeric(N)
  for (k in 2:N){
    h_mumu_1[k] <- -4*(zeta[k]*zeta_expmu[k])^2/zeta_expmu2[k]^2
    dl_mu_1[k] <- 2*zeta[k]*(zeta_expmu[k])^2/zeta_expmu2[k]
  }
  
  # Put all the terms together:
  h_lmblmb <- sum(N*h_lmblmb_1[2:M])-(2/sigma2)*(sum(h_lmblmb_2[2:M]-
                h_lmblmb_3[2:M]-h_lmblmb_4[2:M]+h_lmblmb_5[2:M]))
                                                     
  h_lmbmu <- (2/sigma2)*sum(h_lmbmu_1[2:M]-h_lmbmu_2[2:M])
  h_mulmb <- (2/sigma2)*(sum(h_mulmb_1[2:M]-h_mulmb_2[2:M]))
  h_mumu <- sum(M*h_mumu_1[2:N])+(1/sigma2)*(t(X[,1]) %*% OU_D2BINV(mu,v) %*% X[,1]+
                                               sum(h_mumu_star[2:M]))
  # Define the Hessian matrix:
  H <- rbind(cbind(h_lmblmb,h_lmbmu),cbind(h_mulmb,h_mumu))
  
  # Partial derivatives:
  dl_lmb <- sum(N*dl_lmb_1[2:M]+(2/sigma2)*dl_lmb_2[2:M]-(2/sigma2)*dl_lmb_3[2:M])
  dl_mu <- sum(M*dl_mu_1[2:N])+(1/sigma2)*t(X[,1]) %*% OU_DBINV(mu,v) %*% X[,1]+sum(dl_mu_star[2:M])
  update <- solve(H,rbind(dl_lmb,dl_mu))
  return (update)
}

#Estimation schemes

#Obtaining MLE of sigma^2 using existing lambda hat and mu hat
# OU SIG LIKE.m
# function sigma2 hat = OU SIG LIKE(lambda hat,mu hat,u,v,X)
# This function returns the likelihood estimate of sigma2
# evaluated with mu hat and lambda hat:
  # function sigma2 hat = OU SIG LIKE(lambda hat,mu hat,u,v,X);
# INPUT :
  # lambda hat, mu hat : current estimate of lambda and mu.
# u , v : horizontal and vertical strip respectively.
# X : random field with complete data or missing data imputed.
# OUTPUT:
  # sigma2 hat = MLE of sigma based on input of lambda and mu
OU_SIG_LIKE <- function (lambda_hat,mu_hat,u,v,X){
  N <- dim(X)[1]
  M <- dim(X)[2]
  eta <- matrix(0, nrow=M, ncol=1)
  nu <- matrix(0, nrow=N, ncol=1)
  eta[1,1] <- u[1]
  nu[1,1] <- v[1]
  eta[2:M,1] <- abs(u[2:M]- u[1:(M-1)])
  nu[2:N,1] <- abs(v[2:N]-v[1:(N-1)])
  b <- c(0)
  b[2:N] <- exp(-mu_hat*nu[2:N,1])
  Xa <- matrix(0,nrow=N, ncol=M)
  long_term <- matrix(0, nrow=M, ncol=1)
  Xa[,1] <- X[,1]
  Binv <- OU_COVINV(b)
  for (i in 2:M){
    Xa[,i]<- X[,i]-(exp(-lambda_hat*eta[i,1])*X[,(i-1)])
    long_term[i,1] <- t(Xa[,i]) %*% Binv %*% Xa[,i]/(1-exp(-2*lambda_hat*eta[i,1]))
  }
  
  sigma2_hat <- (t(X[,1]) %*% Binv %*% X[,1] + sum(long_term[2:M,1]))/(M*N)
  
  return (sigma2_hat)
}
  


# Obtaining ALE of , ?? and 2 based on Markov property assumption
# With block-missing observations
# OU APPROXLIKE.m
# This code implements the approximate likelihood function
# based on th Markov property assumption in partitioned data
# Author: Sami Cheong
# Version : 0
# Date : 10/27/2014
# function likelihood=OU APPROXLIKE(X miss,lambda,mu,sigmas)
# The goal is to approximate l mn with l 1 + l 2,
# where l 1 is the likelihood for complete
# column observations, and l 2 is the likelihood for incomplete
# observations

OU_APPROXLIKE <- function(lambda,mu,sigmas,u,v,X){
  
  N <- dim(X)[1]
  M <- dim(X)[2]
  
  # Set of indices for the random field:
  J <- c(1:M)
  K <- c(1:N)
  
  #identify columns with missing data (NaN):
  Miss_mat <- which(is.na(X), arr.ind=T)
  Miss_col <- Miss_mat[,2]        #Miss_col=find(any(Miss mat))
  Miss_row <- Miss_mat[,1]             #Miss_row=find(isnan(X(:,Miss_col(1)))==1)
  Nprime <- length(Miss_row)
  Mprime <- length(Miss_col)
  
  # J is the horizontal axis of the field
  Jprime <- Miss_col
  J2prime <- Jprime
  
  # Kprime is the indices that for the usable information in the
  # missing block columns
  if (Miss_row[1] > ceiling(N/2)){
    Kprime <- c(1:(Miss_row[1]-1))
    K2prime <- Miss_row[length(Miss_row)]+c(1:N)
  }else{
    Kprime <- c((1+Miss_row[length(Miss_row)]):N)
    K2prime <- c(1:Miss_row[1])-1
  }
  
  eta = c()
  nu = c()
  
  if (u[1] == 0){
    eta[1] <- u[1]
  }else{
    eta[1] <- 0
  }
  
  if (v[1] == 0){
    nu[1] <- v[1]
  }else{
    nu[1] <- 0
  }
  
  eta[2:M] <- abs(u[2:M] - u[1:(M-1)])
  nu[2:N] <- abs(v[2:N] - v[1:(N-1)])
  
#covariance structure for the complete data:
  b <- exp(-mu*nu)
  Binv <- OU_COVINV(b)
  
#covariance structure for the incomplete data:
  b0 <- exp(-mu*nu[Kprime])
  B0inv <- OU_COVINV(b0)
  b1 <- exp(-mu*nu[K2prime])
  B1inv <- OU_COVINV(b1)
  
#Define indices and dimensions for l 1 and l 2
  J1 <- setdiff(J, Jprime)
  K1 <- setdiff(K, Kprime)
  J2 <- Jprime
  K2 <- Kprime
  J3 <- J2prime
  K3 <- K2prime
  M1 <- length(J1)
  N1 <- N
  M2 <- length(J2)
  N2 <- length(K2)
  M3 <- length(J3)
  N3 <- length(K3)
  
#Define quadratic term for l_1:
  J11 <- setdiff(J1,1)
  K11 <- setdiff(K1,1)
  Q1 <- 0
  for (j in setdiff(J11, Jprime[length(Jprime)]+1)){
    Xj_star <- X[,j] - exp(-lambda*eta[j])*X[,(j-1)]
    s1 <- (t(Xj_star) %*% Binv %*% Xj_star)/(1-exp(-2*lambda*eta[j]))
    Q1 <- s1 + Q1
  }
  l_1 <- M1*N1*log(2*pi*sigmas) + M1*sum(log(1-exp(-2*mu*nu[2:N]))) +
    N1*sum(log(1-exp(-2*lambda*eta[J11]))) + (t(X[,1]) %*% Binv %*% X[,1] + Q1)/sigmas
  
  #Define the usable observations:
  X0 <- X[K2,J2]
  
  #Define quadratic term for l_2:
  Q2 <- 0
  for (j in 2:M2){
    X0j_star <- X0[,j] - exp(-lambda*eta[J2[j]])*X0[,(j-1)]
    Q2 <- Q2 + ((t(X0j_star) %*% B0inv %*% X0j_star)/(1-exp(-2*lambda*eta[J2[j]])))
  }
  
  K22 <- setdiff(K2, 1)
  J22 <- setdiff(J2,1)
  l_2 <- M2*N2*log(2*pi*sigmas) + M2*sum(log(1-exp(-2*mu*nu[K22]))) +
    N2*sum(log(1-exp(-2*lambda*eta[J22]))) + (t(X0[,1]) %*% B0inv %*% X0[,1]+Q2)/sigmas
  
  #Define the usable observations:
  X1 <- X[K3,J3]
  Q3 <- 0
  for (j in 2:M3){
    X1j_star <- X1[,j] - exp(-lambda*eta[J3[j]])*X1[,(j-1)]
    Q3 <- Q3 + ((t(X1j_star) %*% B1inv %*% X1j_star)/(1-exp(-2*lambda*eta[J2[j]])))
  }
  
  K33 <- setdiff(K3,1)
  J33 <- setdiff(J3,1)
  l_3 <- M3*N3*log(2*pi*sigmas)+M3*sum(log(1-exp(-2*mu*nu[K33])))+
    N3*sum(log(1-exp(-2*lambda*eta[J33])))+(t(X1[,1]) %*% B1inv %*% X1[,1]+Q3)/sigmas
  
  #Sum up the likelihood functions from partitioned data
  likelihood <- l_1 + l_2 + l_3
  
  return(likelihood)
}

# OU_SIM.m
#
  # Function to simulate Gaussian random field
# with Ornstein-Uhlenbeck covariance function 
# 
  # Author  : Sami Cheong
# Date    : 7/29/14
# Version : 1
#
##############################################
# function [X,A,B,Gamma]=OU_SIM(u,v,lambda,mu,sigma2)
#
# INPUT: 
  #
# u                    = horizontal coordinate of input grid between [0,1]
# v                    = vertical coordinate of input grid between [0,1]
# cov(X(u,v),X(u',v')) = sigma2*exp(-lambda|u-u'|-mu|v-v'|)
#
# OUTPUT: 
  #
# X = an N-by-M array of a random field realization.
# A = Covariance matrix for the 'horizonal' part of the random field.
# B = Covariance matrix for the 'vertical' part of the random field.
# Gamma = Covariance matrix for the whole random field with complete
# observations.
#

#################################################
OU_SIM <- function(u,v,lambda,mu,sigma2){
  

# Define the length of each of the input vectors 
M <- length(u)
N <- length(v)

#Initialize and define the distance between consecutive input points:
eta <- matrix(0, nrow=M, ncol=1)
nu  <- matrix(0, nrow=N, ncol=1)

if(u[1]==0){
  eta[1] <- u[1]
}else{
  eta[1] <- 0
}

if (v[1]==0){
  nu[1] <- v[1]
}else{
  nu[1] <- 0
}

eta[2:M] <- abs(u[2:M - u[1:(M-1)]])
nu[2:N] <- abs(v[2:N] - v[1:(N-1)])

#Define the component of the OU covariance matrix 
#based on distance between grid points:

a <- exp(-lambda*eta)
b <- exp(-mu*nu)

#a(1:M) = exp(-lambda.*eta(1:M));
#b(1:N) = exp(-mu.*nu(1:N));

#Initialize and define the covariance function 
#for horizonal and vertical components:

A <- diag(M)
B <- diag(N)

#Assign elements to the covariance matrices 

for (i in 1:(M-1)){
  A[i,(i+1)] <- a[i+1]
  for (j in (i+1):M){
    A[i,j] <- prod(a[(i+1):j])
    A[j,i] <- A[i,j]
  }
}

for (i in 1:(N-1)){
  B[i,(i+1)] <- b[i+1]
  for (j in (i+1):N){
    B[i,j] <- prod(b[(i+1):j])
    B[j,i] <- B[i,j]
    }
}

#Covariance matrix for X as sigma2*Kronecker product of A and B,
#where X is represented as [X_1,...X_M]'. 
Gamma <- kronecker(A,B)
Gamma <-sigma2*Gamma

#Simulate Multivariate normal observations
#based on the covariance matrix Gamma
#The covariance matrix has dimension K-by-K, where K = M*N.
K <- M*N

#Generate K Normal(0,1) random vairables
Z <- matrix(rnorm(K,mean=0,sd=1), nrow=K, ncol=1)

#Cholesky decompision of Gamma, where L is a lower triangular
#matrix such that L'*L=Gamma
L <- chol(Gamma)

#X now is Normal with covariance matrix Gamma.
X <- t(L)%*%Z
X <- matrix(X, nrow=N, ncol=M)
out_list <- list("X"=X,"A"=A,"B"=B,"Gamma"=Gamma)
return(out_list)
}


#Obtaining MLE of lambda, ?? and sigma2 using Newton's method
#OU LIKE NEWTON.m
#function [theta] = OU LIKE NEWTON(theta0, delta)
# Function to implement Newton's method on l(X) to obtain MLE
# where X is an OU process with complete observation
# INPUT:
# lambda0,mu0,sigma20 = initial guess for the three parameter
# values
# u,v = input grid for the random field
# X = set of observations from which we wish to
# approximate the parameter
# delta = accuracy we set for the estimate,
# delta = j j theta???theta0 j j 2

OU_LIKE_NEWTON <- function(lambda0,mu0,sigma20,u,v,X, delta){
  #format long e
  lambda0 <- OU_SUB_reset_bound(lambda0)
  mu0 <- OU_SUB_reset_bound(mu0)
  sigma20 <- OU_SUB_reset_bound(sigma20)
  X[is.na(X)] <- 0
  #Evaluate the initial value wrt the complete data likelihood
  l_0 = OU_LIKE(lambda0,mu0,sigma20,u,v,X)
  if (abs(l_0) <= delta){
    #check to see if initial guess satisfies
    return(l_0)
  }
  
  #############################
  ##
  ## MAIN ROUTINE
  ##
  #############################
  
  max_iter <- 2000
  iter <- 0
  error <- 1000
  while (error > delta & iter < max_iter){
    
  
  l_0 <- OU_LIKE(lambda0, mu0, sigma20, u, v, X)
  #update parameters lambda and mu using Newton's method
  theta_update <- rbind(lambda0,mu0) - OU_LIKE_HESS_UPDATE(lambda0,mu0,sigma20,u,v,X)
  mu0 <- OU_SUB_reset_bound(theta_update[2])
  
  #update sigma2
  sigma20 <- OU_SIG_LIKE(lambda0,mu0,u,v,X)
  sigma2_update <- sigma20
  
  #measure error in the likelihood function
  error <- abs(l_0 - OU_LIKE(lambda0,mu0,sigma20,u,v,X))
  
  #update iteration
  iter = iter + 1
  # print stuff
  #fprintf('nn Newton iteration = %d, delta = %d,nn', iter,error)
  #fprintf('nn lambda = %d, mu = %d, sigma2 = %d, nn', ...
  # lambda0, mu0,sigma20)
  #theta = [theta update(1);theta update(2);sigma2 update];
  list_out <- list("lambda"=theta_update[1], "mu"=theta_update[2],"sigma2"=sigma2_update,"likelihood"=l_0)
  }
  
  return (list_out)
}
  
OU_SUB_reset_bound <-  function(x){
  if (x <= 1 | x >100){
    x <- 2 + sample(1:5,1)
    print ('/n Parameter values reset to default,value = %d, /n',x)
  }
  theta_reset <- x
  return(theta_reset)  
}


#Simulation of missing blocks in the random field

# OU BLOCKMISS.m
# Function to simulate a block of missing observations in a random
# field
# function [Y,K1,J1] = OU BLOCKMISS1(X,missing prop)
#
  # INPUT :
  # X = matrix of realization of a random field
# with complete observations
# missing prop = proportion of the observations that are missing

 OU_BLOCKMISS <- function(X,missing_prop){
   if (missing_prop <= 0 | missing_prop >=0.45)
     {print('missing level must be strictly between 0 and 0.45')}
   # get dimension:
   N <- dim(X)[1]
   M <- dim(X)[2]
   # Total number of observations
   # generate the block dimension according to level of missingness:
     # C is the area of the missing block,
   # we need to generate dimensions of the missing block
   Y <- X
   C <- round(1/missing_prop)
   factor <- c(9, 8, 7, 6, 5, 4, 3, 2)
   i <- 1
   
   while (i < length(factor)){
     if (C%%factor[i] == 0){
       C1 <- C/factor[i]
       C2 <- C/C1
       i <- i+1
     }else{
       C1 <- floor(C/3)
       C2 <- floor(C/C1)
       i <- i + 1
     }
   }
   
   N1 <- round(N/C1)
   M1 <- round(M/C2)
   
   #generate an index to start the missing block:
   k1 <- sample(3:floor(N/4),1)
   j1 <- sample(3:floor(M/4),1)
   K1 <- k1:min((k1+(N1-1)), N)
   J1 <- j1:min((j1+(M1-1)), M)
   
   #Assign NA to the resulting block
   Y[K1,J1] <- NA
   out_list <- list("Y"=Y,"K1"=K1,"J1"=J1)
   return(out_list)
 }
 
#Code validation Fei 9/20/2016
 library(raster)
 library(sp)
#Generate X
 M <- 39
 N <- 33
 MN <- M*N
 Z <- rnorm(MN, 0, 1)
 u <- seq(0,1,1/(M-1))
 v <- seq(0,1,1/(N-1))
 lambda <- 3.2
 mu <- 5.1
 sigma2 <- 4.18
 out_list<-OU_SIM(u,v,lambda,mu,sigma2)
 X <- out_list$X
 A_true <- out_list$A
 B_true <- out_list$B
 Gamma_true <- out_list$Gamma
 X_true <- X
 image(v,u,X)

 #Remove block
 missing_prop <- 0.04
 out_list <- OU_BLOCKMISS(X,missing_prop)
 Y <- out_list$Y
 K1 <- out_list$K1
 J1 <- out_list$J1
 image(v,u,Y)
 #write.csv(Y,"missingField.csv",row.names=FALSE)

 #Simulate missing data: estimate parameters
 
OU_LIKE_test<-function(x){
  OU_LIKE(lambda,x,sigma2,u,v,X)
}


