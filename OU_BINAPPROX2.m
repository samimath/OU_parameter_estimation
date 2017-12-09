%% function [approx_like,approx_err] = OU_BINAPPROX2(X_miss,u,v,A,B,mu,lambda,sigma2,num_term_approx)
%
% This codes implements the approximation likelihood function for a
% realization of an OU field with missing observations
% Definition: 
% f(Xo|\theta) = f(Xo_1)*\Pi_{j=2}^{m}f(Xo_j|Xo_{j-1})
% where each of Xo_j|Xo_{j-1} is multivariate normal with mean mo_j and
% covariance matrix Sigmao_j, where 
% mo_j         = e^{-\lambda|u_j-u_{j-1}|}B_{j,j-1}(B_{j-1,j-1})^{-1}Xo_{j-1}
% Sigmao_j     = \sigmas(B_{j,j}-e^{-2\lambda|u_j-u_{j-1}|}B_{j,j-1}(B_{j-1,j-1})^{-1}B_{j-1,j})
% We propose to approximate the inverse of Sigmao_j using power series
% expansion:
% Sinv_j       = (L_j)^{-1}(I + T_j + (T_j)^2 +(T_j)^3)(L_j)^{-1}
% where T_j    = e^{-2\lambda|u_j-u_{j-1}|}(L_jB_{j,j-1}L_{j-1})
% The negative log likelihood is defined as :
% l(Xo|\theta) = (sum_{j=1}^{m}n_j) log(2*pi*\sigmas)
%                 + sum_{j=1}^{m} log(det(Sigmao_j))
%                 + \frac{1}{\sigmas}(Xo_1^{T}(Sigmao_1)^{-1}Xo_1 
%                   + sum_{j=2}^{m} Xo_j^{T}(Sinv_j)Xo_j)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT: 
% X_miss            : set of available observations in the OU field arranged in a 2D
%                     matrix
% u, v              : horizontal and vertical input grids
% A, B              : Covariance matrices for the horizontal and vertical components
% mu,lambda,sigmas  : parameters for the model  
% num_term_approx   : number of terms to use in the power series
% row_col           : if row_col ='row', then partition is by row, if 'col' then partition is by col 
% OUTPUT:
% approx_like       : value of the approximated likelihood based on the
%                     description above 
% approx_error      : difference between the direct inverse and the
%                     approximation 
% Version : 1
% Date    : 7/23/15
% Author  : Sami Cheong 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [approx_log_like,approx_log_like_est,approx_err] = OU_BINAPPROX2(X_miss,u,v,A,B,mu,lambda,sigma2,num_term_approx)
% initialize components of the log likelihood 
approx_err_temp=0;

M = size(X_miss,2);
N = size(X_miss,1);
log_like_vec=NaN(M,1);
log_like_est=NaN(M,1);
time_cov_temp=NaN(M,1);
time_est_cov_temp=NaN(M,1);
quadratic_term = NaN(M,1);
num_of_obs = NaN(M,1);
% get set of indices with available observations:
v1_avail_ind = find(~isnan(X_miss(:,1)));
Xo_1    = X_miss(v1_avail_ind,1);
B11     = B(v1_avail_ind,v1_avail_ind);  
v1_avail     = v(v1_avail_ind);
% standardize vectors to be column vectors 
if size(v1_avail,2)~=1
    v1_avail = v1_avail';
    else
end
% The first column of X_miss is treated as a special case since it does not have
% conditional density
zeta_11 = [0;abs(v1_avail(2:end)-v1_avail(1:end-1))];
B11_inv = OU_COVINV(exp(-mu.*zeta_11));
num_of_obs(1)      = length(Xo_1);
quadratic_term(1)  = Xo_1'*B11_inv*Xo_1;
% assign values to the negative log-likelihood function 
log_like_vec(1)    = num_of_obs(1)*log(2*pi*sigma2) + ...
                     det(log(B11)) + quadratic_term(1);
log_like_vec_est(1)    = num_of_obs(1)*log(2*pi*sigma2) + ...
                     det(log(B11)) + quadratic_term(1);
%% For the rest of the columns:
for j=2:M
% find indices of available sites: 
vj_avail_ind  = find(~isnan(X_miss(:,j)));
vj1_avail_ind = find(~isnan(X_miss(:,j-1)));
% get the corresponding data: 
Xo_j          = X_miss(vj_avail_ind,j);
Xo_j1         = X_miss(vj1_avail_ind,j-1);

% covariance matrices:
B_jj    = B(vj_avail_ind,vj_avail_ind);
B_jj1   = B(vj_avail_ind,vj1_avail_ind);
B_j1j   = B_jj1';
% number of available observations for column j:
num_of_obs(j) = length(Xo_j);
vj_avail      = v(vj_avail_ind);
vj1_avail      = v(vj1_avail_ind);
if size(vj_avail,2)~=1
    vj_avail = vj_avail';
    else
end
if size(vj1_avail,2)~=1
    vj1_avail = vj1_avail';
    else
end
zeta_jj = [0;abs(vj_avail(2:end)-vj_avail(1:end-1))];
zeta_j1j1 = [0;abs(vj1_avail(2:end)-vj1_avail(1:end-1))];
% inverse of B_jj and B_j1j1:
B_jj_inv = OU_COVINV(exp(-mu.*zeta_jj));
B_j1j1_inv = OU_COVINV(exp(-mu.*zeta_j1j1));
% LU decomposition of B_j1j1_inv
L_j = OU_SQRTM(B_jj_inv);
L_j1= OU_SQRTM(B_j1j1_inv);

eta_j  =abs(u(j)-u(j-1));
% define conditional mean:
mo_j   = exp(-lambda*eta_j)*B_jj1*B_j1j1_inv*Xo_j1;
% define conditional covariance matrix:
So_j = B_jj - exp(-2*lambda*eta_j)*B_jj1*B_j1j1_inv*B_j1j;

% define the term used in the power series expansion:
T_j = L_j'*B_jj1*L_j1;
Tstar_j=exp(-2*lambda*eta_j)*T_j*(T_j');
Tstar_j_terms=NaN(size(Tstar_j,1),size(Tstar_j,2),num_term_approx);
for k=1:num_term_approx
Tstar_j_terms(:,:,k)=Tstar_j^k;
end
Tstar_j_sum=sum(Tstar_j_terms,3);


I_j=eye(num_of_obs(j),num_of_obs(j));

Ao_j = ((L_j')\(I_j-Tstar_j))/(L_j);

So_j_inv =inv(So_j);

% Approximate the inverse of the conditional variance:
Ao_j_inv = L_j*(I_j+Tstar_j_sum)*(L_j');
%{
err_SA=norm(So_j-Ao_j);
err_SA_inv=norm(So_j_inv-Ao_j_inv);
err_qterm=(Xo_j-mo_j)'*So_j_inv *(Xo_j-mo_j)-(Xo_j-mo_j)'*Ao_j_inv
*(Xo_j-mo_j);
%}
l1=num_of_obs(j)*log(2*pi*sigma2)+log(det(So_j))+(1/sigma2)*((Xo_j-mo_j)'*So_j_inv *(Xo_j-mo_j));
l2=num_of_obs(j)*log(2*pi*sigma2)+log(det(Ao_j))+(1/sigma2)*((Xo_j-mo_j)'*Ao_j_inv' *(Xo_j-mo_j));
%err_loglike(j)=l1-l2;

quadratic_term(j)=(Xo_j-mo_j)'*Ao_j_inv *(Xo_j-mo_j);
%log_like_vec_est(j)=l2;
log_like_vec(j)=l2;

approx_err_temp=approx_err_temp+norm(So_j_inv-Ao_j_inv);





end
approx_err=approx_err_temp
approx_log_like=sum(log_like_vec);
approx_log_like_est=sum(log_like_vec_est);
err_ll=approx_log_like-approx_log_like_est

