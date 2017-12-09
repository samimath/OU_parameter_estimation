%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OU_MISS.m
%% Function to generate missing observations for OU_SIM
%% Author  : Sami Cheong
%% Date    : 7/29/14
%% Version : 1
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function X_miss= OU_MISS(X,miss_level)
% INPUT:
% X          = complete-observation Gaussian random field.
% miss_level = % of the observations that is missing 
% OUTPUT:
% X_miss     = X with randomly missing values according to miss_level
function [X_miss,miss_ind]= OU_MISS(X,miss_level)
N=size(X,1);
M=size(X,2);
X=reshape(X,[N*M,1]);
% Create missing observations using Bernoulli distribution with
% probability defined by 'miss_level':
miss_ind=binornd(1,miss_level,N*M,1);
% Or randomly permute the sampling sites:
% perm_ind = randperm(N*M);
% Choose the % of observations missing as represented by the index:
%miss=perm_ind(1:floor(miss_level*(N*M)));
X_miss=X;
X_miss(miss_ind==1)=NaN;
X_miss=reshape(X_miss,[N,M]);
miss_ind=reshape(miss_ind,[N,M]);