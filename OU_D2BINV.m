%% This function computes the 2nd derivtive 
%% of the inverse of the covariance matrix B^-1(\mu):
% function DB2_inv = OU_DBINV(mu,v)
function DB2_inv = OU_D2BINV(mu,v)
% INPUT:
% mu = parameter value for B(mu)
% v  = Vector for the vertical coordinate of the random field
% OUTPUT:
% DB_inv = A matrix of derivatives for B^-1(\mu) 
% Assign values to border elements:
N = length(v);
DB2_inv=zeros(N,N);
% standaridze vectors into column vectors 
if size(v,2)~=1
   v=v;
else
end
nu =[0; abs(v(2:N)-v(1:N-1))];
b = exp(-mu.*nu);
b_sq=b.^2;
b_star=(nu.*b_sq.*(ones(N,1)+b_sq))/((ones(N,1)-b_sq).^3);
DB2_inv(1,1) = 4*b_star(2);
DB2_inv(N,N) = 4*b_star(N);
for j=2:N-1
    DB2_inv(j,j)  = 4*(b_star(j)+b_star(j+1));
end
for j=2:N    
    DB2_inv(j-1,j) = ...
        ((nu(j))^2*b(j)*(1+6*(b(j))^2+(b(j))^4))/((1-b_sq(j))^2);
    DB2_inv(j,j-1) =DB2_inv(j-1,j);
end