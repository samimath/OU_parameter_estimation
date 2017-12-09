%% This function finds the lower-bidiagonal 
%% matrix square-root for a symmatric tridiagonal matrix 
% INPUT:
% G = symmetric tridiagonal matric
% OUTPUT:
% L = lower bidiagonal matrix such that G=L*L'
function [L]=OU_SQRTM(G)
if norm(G-G') > 0
    error('G must be symmetric and tridiagonal')
else
end
K=size(G,1);
L=zeros(K,K);
L(1,1) = sqrt(G(1,1));
L(K,K) = sqrt(G(K,K));
for j=2:K
   L(j,j-1)=G(j-1,j)/L(j-1,j-1);
   L(j,j)  =sqrt(G(j,j)-L(j,j-1)^2);
end