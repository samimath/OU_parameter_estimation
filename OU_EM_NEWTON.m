%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to implement the EM algorithm via Newton's method.
%% Author  : Sami Cheong
%% Date    : 7/29/14
%% Version : 1
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [theta_new,lmb_vec,sig_vec,mu_vec,iter_vec,lmn_vec] = OU_EM(X_miss,u,v,lambdap,mup,sigma2p)
% INPUT : 
% X_miss                 = OU field with randomly missing observations.
% u                      = horizontal coordinate of input grid between [0,1]
% v                      = vertical coordinate of input grid between [0,1]
% lambdap, mup, sigma2p  = initial parameter values for the algorithm.
% OUTPUT:
% theta_new              = EM estimation of the parameter lambda, mu, sigma2
% lmb_vec,sig_vec,mu_vec = Vector of the EM estimates at each iteration.
% iter_vec               = vector of indices that keeps track of steps.
%
function [theta_new,lmb_vec,sig_vec,mu_vec,iter_vec,lmn_vec] = OU_EM_NEWTON(X_miss,u,v,lambdap,mup,sigma2p)
% set tolerance and max number of iterations 
tol=0.0001;
max_step=100;
iter=1;
%sigma_int=sigma2p;
% Initial parameter values:
%theta_p=[lambdap;mup;sigma2p];
theta_new=zeros(3,1);
err=1000;
% Initialize:
lmb_vec = NaN(max_step,1);
sig_vec = NaN(max_step,1);
mu_vec  = NaN(max_step,1);
lmn_vec=NaN(max_step,1);
err_new=2000;
while abs(err_new-err) > tol && iter < max_step   
    err_new=err;
    % Evaluate covariance matrix wrt current parameter value:
    [~,~,~,Gamma]=OU_SIM(u,v,OU_SUB_reset_bound(lambdap),OU_SUB_reset_bound(mup),OU_SUB_reset_bound(sigma2p));
    % Generate conditional random field:
    [X_cond,~,~,~,~]=OU_COND(X_miss,Gamma);
    % Minimize the likelihood function with current parameter value:
   [theta_new] = OU_LIKE_NEWTON(lambdap,mup,sigma2p,u,v,X_cond,tol);

    lmb_vec(iter)=lambdap;
    mu_vec(iter)=mup;
    sig_vec(iter)=sigma2p;
    
    lambdap = theta_new(1);
    mup = theta_new(2);
    sigma2p = theta_new(3);

    
    % Keep track of the value of the likelihood function 
    lmn_vec(iter)=OU_LIKE(theta_new(1),theta_new(2),theta_new(3),u,v,X_cond);
    theta_p=[lambdap;sigma2p;mup]; 
    err = sum(abs(theta_new-theta_p));
     iter=iter+1;
    
    
end


[~,~,~,Gamma]=OU_SIM(u,v,theta_new(2),theta_new(3),theta_new(3));
% Generate conditional random field:
[X_cond,~,~,~,~]=OU_COND(X_miss,Gamma);
theta_new(3)= OU_SIG_LIKE(theta_new(1),theta_new(2),u,v,X_cond);
    
iter_vec=1:iter-1;
lmb_vec=lmb_vec(iter_vec);
mu_vec=mu_vec(iter_vec);
sig_vec=sig_vec(iter_vec);
lmn_vec=lmn_vec(iter_vec,:);

%% SUB function to control initial values:
function theta_reset = OU_SUB_reset_bound(x)
	if (x <= 1 || x > 100)
    x = 2+randi(5);
fprintf('\n Parameter values reset to default,value = %d, \n',x)
else
end
theta_reset = x;
                                      
return;

