%% function hess = OU_LIKE_HESSIAN(lambda,mu,sigma2,u,v,X)
%% function to evaluate the hessian matrix of the complete data likelihood 
%% generate structure :H=[h_lmblmb, h_lmu; h_mul, h_mumu]
%% Uses other functinos : OU_DBINV(mu,v) , OU_D2BINV(mu,v)
function [H] = OU_LIKE_HESSIAN(lambda,mu,sigma2,u,v,X)
% H=[H_11, H_12; H_21, H_22]
% define the components of the Hessian matrix: 
M=size(X,2);
N=size(X,1);
% standardize things to be column vectors : 
if size(u,2)~=1
u=u';
else
end
if size(v,2)~=1
v=v';
end

eta=[0;abs(u(2:end)-u(1:end-1))];
zeta=[0;abs(v(2:end)-v(1:end-1))];

zeta_expmu  = exp(-mu.*zeta);
zeta_expmu2 = ones(N,1)-(zeta_expmu.^2);
eta_explmb  = exp(-lambda.*eta);
eta_explmb2 = ones(M,1)-(eta_explmb.^2);

% compute B^-1:
B_inv=OU_COVINV(zeta_expmu);
% initialize stuff
%X_q=NaN(size(X));
for j=2:M
% define terms in H_11 = h_lambdalambda
X_q(:,j)=X(:,j)-eta_explmb(j)*X(:,j-1);
h_lmblmb_1(j) =  -4*(eta(j)*eta_explmb(j))^2/((eta_explmb2(j))^2);
h_lmblmb_2(j) = ((eta(j))^2*eta_explmb(j)*(1+(eta_explmb(j))^2)*(X(:,j-1)')*B_inv*X_q(:,j))/((eta_explmb2(j))^2);
h_lmblmb_3(j) = (eta(j)*eta_explmb(j))^2*(X(:,j-1)')*B_inv*X(:,j-1)/(eta_explmb2(j));
h_lmblmb_4(j) = 2*(eta(j)*eta_explmb(j))^2*(1+(eta_explmb(j))^2)*X_q(:,j)'*B_inv*X_q(:,j)/((eta_explmb2(j))^3);
h_lmblmb_5(j) = 2*(eta(j)*eta_explmb(j))^2*eta_explmb(j)*(X(:,j-1)')*B_inv*X_q(:,j)/((eta_explmb2(j))^2);

% split the H_12 = h_lambdamu:
h_lmbmu_1(j) = (eta(j)*eta_explmb(j)*X(:,j-1)'*OU_DBINV(mu,v)*X_q(:,j))/(eta_explmb2(j));

h_lmbmu_2(j) = ((eta(j)*(eta_explmb(j))^2)*X_q(:,j)'*OU_DBINV(mu,v)*X_q(:,j))/((eta_explmb2(j))^2);

h_mulmb_1(j) = (eta(j)*((eta_explmb(j))^2)*X(:,j-1)'*OU_DBINV(mu,v)*X_q(:,j))/(eta_explmb2(j));
h_mulmb_2(j) = (eta(j)*((eta_explmb(j))^2)*X_q(:,j)'*OU_DBINV(mu,v)*X_q(:,j));

% define the lambda term in h_mumu:
h_mumu_star(j) = (X_q(:,j)'*OU_D2BINV(mu,v)*X_q(:,j))/(eta_explmb2(j));
end 

for k=2:N
h_mumu_1(k)=-4*(zeta(k)*zeta_expmu(k))^2/((zeta_expmu2(k))^2);
end

h_lmblmb = N*sum(h_lmblmb_1(2:M))+...
          -(2/sigma2)*(sum(h_lmblmb_2(2:M)-h_lmblmb_3(2:M)-h_lmblmb_4(2:M)+h_lmblmb_5(2:M)));

h_lmbmu=(2/sigma2)*sum(h_lmbmu_1(2:M)-h_lmbmu_2(2:M));

h_mulmb = (2/sigma2)*(sum(h_mulmb_1(2:M)-h_mulmb_2(2:M)));

h_mumu = M*sum(h_mumu_1(2:N)) +...
       (1/sigma2)*(X(:,1)'*OU_D2BINV(mu,v)*X(:,1) + sum(h_mumu_star(2:M)));
H=[h_lmblmb,h_lmbmu; h_mulmb h_mumu];
