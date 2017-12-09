%% Evaluate the conditonal mean and variance given incomplete X
% OU_COND.m
% function [X_cond,S_oo,S_uu,S_ou,S_uo]=OU_COND(X,Gamma)
% INPUT:
% X        = N-by-M centered Gaussian random (OU) field 
% (with missing observations).
% Gamma    = Covariance structure of the OU field, 
% evaluated at the current parameter value.
%
% OUTPUT:
% X_cond   = Random field where missing observations are replaced with
% conditonal mean. 
% Cov_cond = conditional covariance matrix for the unobserved samples
function [X_cond,S_oo,S_uu,S_ou,S_uo]=OU_COND(X,Gamma)
% reshape X into a MN-by-1 vector:
N=size(X,1);
M=size(X,2);
X= reshape(X,[N*M,1]);
% get unobserved and observed indices:
[Unobs_ind] = find(isnan(X)==1);
[Obs_ind]=find(isnan(X)==0);
% Partition the covariance matrix S= [S_oo | S_ou; S_uo | S_uu]:
S_oo = Gamma(Obs_ind,Obs_ind);
S_uu = Gamma(Unobs_ind, Unobs_ind);
S_ou = Gamma(Obs_ind,Unobs_ind);
S_uo=S_ou';
X_obs=X(Obs_ind);

Xstar=S_oo\X_obs;
X_unobs = S_uo*Xstar ;
X_cond=X;
X_cond(Unobs_ind)=X_unobs;
X_cond=reshape(X_cond,[N,M]);
