%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OU_SIM.m
%%
%% Function to simulate Gaussian random field
%% with Ornstein-Uhlenbeck covariance function 
%% 
%% Author  : Sami Cheong
%% Date    : 7/29/14
%% Version : 1
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [X,A,B,Gamma]=OU_SIM(u,v,lambda,mu,sigma2)
%
% INPUT: 
%
% u                    = horizontal coordinate of input grid between [0,1]
% v                    = vertical coordinate of input grid between [0,1]
% cov(X(u,v),X(u',v')) = sigma2*exp(-lambda|u-u'|-mu|v-v'|)
%
% OUTPUT: 
%
% X = an N-by-M array of a random field realization.
% A = Covariance matrix for the 'horizonal' part of the random field.
% B = Covariance matrix for the 'vertical' part of the random field.
% Gamma = Covariance matrix for the whole random field with complete
% observations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,A,B,Gamma]=OU_SIM(u,v,lambda,mu,sigma2)
% Define the length of each of the input vectors 
M = length(u);
N = length(v);

%% Initialize and define the distance between consecutive input points:
eta = zeros(M,1); 
nu  = zeros(N,1);
if u(1)==0
    eta(1) = u(1);
else
    eta(1) = 0;
end
if v(1)==0
    nu(1)  = v(1);
else
    nu(1)=0;
end
eta(2:M) = abs(u(2:M)-u(1:M-1));
nu(2:N)  = abs(v(2:N)-v(1:N-1));
%% Define the component of the OU covariance matrix 
%% based on distance between grid
% points:
a = exp(-lambda.*eta);
b = exp(-mu.*nu);
%a(1:M) = exp(-lambda.*eta(1:M));
%b(1:N) = exp(-mu.*nu(1:N));

%% Initialize and define the covariance function 
%% for horizonal and vertical components:
A = eye(M,M);
B = eye(N,N);
% Assign elements to the covariance matrices 
for i = 1:M-1
    A(i,i+1) = a(i+1);
    for j = i+1:M
        A(i,j) = prod(a(i+1:j));
        A(j,i) = A(i,j);
    end
end
for i = 1:N-1
    B(i,i+1) = b(i+1);
    for j = i+1:N
        B(i,j) = prod(b(i+1:j));
        B(j,i) = B(i,j);
    end
end
%% Covariance matrix for X as sigma2*Kronecker product of A and B,
% where X is represented as [X_1,...X_M]'. 
[Gamma] = kron(A,B);
Gamma=sigma2*Gamma;
%% Simulate Multivariate normal observations
% based on the covariance matrix Gamma
% The covariance matrix has dimension K-by-K, where K = M*N.
K = M*N;
% Generate K Normal(0,1) random vairables. 
Z = randn(K,1);
% Cholesky decompision of Gamma, where L is a lower triangular
% matrix such that L'*L=Gamma.
L = chol(Gamma);
% X now is Normal with covariance matrix Gamma.
X = L'*Z;
X=reshape(X,[N,M]);