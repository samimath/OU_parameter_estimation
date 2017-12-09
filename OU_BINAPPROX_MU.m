%% OU_BINAPPROX.m
%% function [approx_like,approx_err] = 
%% OU_BINAPPROX_MU(X_miss,u,v,A,B,mu,lambda,sigma2,num_term_approx)
%
% This codes implements the approximation likelihood function for a
% realization of an OU field with missing observations
% uses horizontal partition instead of vertical partition 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT: 
% X_miss            : set of available observations in 
%                   the OU field arranged in a 2D matrix
% u, v              : horizontal and vertical input grids
% A, B              : Covariance matrices for the horizontal
%                   and vertical components
% mu,lambda,sigmas  : parameters for the model  
% num_term_approx   : number of terms to use in the power series
% calc_error        : indicates whether to return the error of
% approximating the inverse of the conditional covariance matrix
% OUTPUT:
% approx_like       : value of the approximated likelihood based on the
%                     description above 
% approx_error      : difference between the direct inverse and the
%                     approximation 
% Version : 1
% Date    : 8/17/15
% Author  : Sami Cheong 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [approx_log_like,approx_err] = OU_BINAPPROX_MU(X_miss,u,v,A,B,...
                                        mu,lambda,sigma2,...
                                        num_term_approx,calc_error)
% initialize components of the log likelihood 
approx_err_temp=0   ;
M = size(X_miss,2);
N = size(X_miss,1);
log_like_vec=NaN(N,1);
quadratic_term = NaN(N,1);
num_of_obs = NaN(N,1);
% get set of indices with available observations:
u1_avail_ind = find(~isnan(X_miss(1,:)));
Xo_1    = X_miss(1,u1_avail_ind);
% work with column vector:
if size(Xo_1,2)~=1
    Xo_1=Xo_1';
end
A11     = A(u1_avail_ind,u1_avail_ind);  

u1_avail     = u(u1_avail_ind);
% standardize vectors to be column vectors 
if size(u1_avail,2)~=1
    u1_avail = u1_avail';
    else
end
% The first column of X_miss is treated as a ...
% special case since it does not have conditional density on previous col
eta_11 = [0;abs(u1_avail(2:end)-u1_avail(1:end-1))];
A11_inv = OU_COVINV(exp(-lambda.*eta_11));
num_of_obs(1)      = length(Xo_1);
quadratic_term(1)  = Xo_1'*A11_inv*Xo_1;
% assign values to the negative log-likelihood function 
log_like_vec(1)    = num_of_obs(1)*log(2*pi*sigma2) + ...
                     det(log(A11)) + quadratic_term(1);
%% For the rest of the columns:
for j=2:N
% find indices of available sites: 
uj_avail_ind  = find(~isnan(X_miss(j,:)));
uj1_avail_ind = find(~isnan(X_miss(j-1,:)));
% get the corresponding data: 
Xo_j          = Sub_colvec(X_miss(j,uj_avail_ind));
Xo_j1         = Sub_colvec(X_miss(j-1,uj1_avail_ind));

% covariance matrices:
A_jj    = A(uj_avail_ind,uj_avail_ind);
A_jj1   = A(uj_avail_ind,uj1_avail_ind);
A_j1j   = A_jj1';
% number of available observations for column j:
num_of_obs(j) = length(Xo_j);

% get distance between each sampling sites:
uj_avail      = u(uj_avail_ind);
uj1_avail      = u(uj1_avail_ind);
% standardize the vectors to be column vectors:
if size(uj_avail,2)~=1
    uj_avail = uj_avail';
    else
end
if size(uj1_avail,2)~=1
    uj1_avail = uj1_avail';
    else
end
% components for the covariance function:
eta_jj = [0;abs(uj_avail(2:end)-uj_avail(1:end-1))];
eta_j1j1 = [0;abs(uj1_avail(2:end)-uj1_avail(1:end-1))];

% inverse of B_jj and B_j1j1:
A_jj_inv = OU_COVINV(exp(-lambda.*eta_jj));
A_j1j1_inv = OU_COVINV(exp(-lambda.*eta_j1j1));

% LU decomposition of B_j1j1_inv
L_j = OU_SQRTM(A_jj_inv);
L_j1= OU_SQRTM(A_j1j1_inv);
% Define the distance between each column
zeta_j  =abs(v(j)-v(j-1));
% Define conditional mean:
mo_j   = exp(-mu*zeta_j)*A_jj1*A_j1j1_inv*Xo_j1;
% Define conditional covariance matrix So_j:
So_j = A_jj - exp(-2*mu*zeta_j)*A_jj1*A_j1j1_inv*A_j1j;
So_j_inv=inv(So_j);
% Define the terms used in the power series expansion:
T_j = L_j'*A_jj1*L_j1;
Tstar_j=exp(-2*mu*zeta_j)*T_j*(T_j');
% initialize the power sum
Tstar_j_terms=NaN(size(Tstar_j,1),size(Tstar_j,2),num_term_approx);
% The power sum depends on user input num_term_approx:
for k=1:num_term_approx
Tstar_j_terms(:,:,k)=Tstar_j^k;
end
Tstar_j_sum=sum(Tstar_j_terms,3);
% Identify matrix used for the power series expansion
I_j=eye(num_of_obs(j),num_of_obs(j));
% So_j expressed in a different form:
Ao_j = ((L_j')\(I_j-Tstar_j))/(L_j);
% Approximate the inverse of the conditional variance:
Ao_j_inv = L_j*(I_j+Tstar_j_sum)*(L_j');
% Approximated quadratic form:
quadratic_term(j)=(Xo_j-mo_j)'*So_j_inv*(Xo_j-mo_j);
% log likelihood for the jthe column:
log_like_vec(j)=num_of_obs(j)*log(2*pi*sigma2)+log(det(So_j))+...
                (1/sigma2)*quadratic_term(j);
% keep track of error between true inverse and power series expansion
%
if calc_error == 1
approx_err_temp=approx_err_temp+mean(mean(So_j_inv-Ao_j_inv));
else 
end
end
approx_err=(approx_err_temp)/M;
approx_log_like=sum(log_like_vec);
end
    function v = Sub_colvec(v)
if size(v,2)~=1
    v=v';
end
return
end

