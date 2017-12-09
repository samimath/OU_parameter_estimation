%% OU_LIKE_NEWTON.m
% function [theta] = OU_LIKE_NEWTON(theta0, delta)
%
% Function to implement Newton's method on l(X) to obtain MLE
% where X is an OU process with complete observation 
% INPUT:
% lambda0,mu0,sigma20 = initial guess for the three parameter
%                       values
% u,v                 = input grid for the random field 
% X                   = set of observations from which we wish to 
%                       approximate the parameter 
% delta               = accuracy we set for the estimate, 
%                       delta = ||theta-theta0||_2 
function [theta] = OU_LIKE_NEWTON(lambda0,mu0,sigma20,u,v,X, delta)
%format long e
lambda0=OU_SUB_reset_bound(lambda0);
mu0=OU_SUB_reset_bound(mu0);
sigma20=OU_SUB_reset_bound(sigma20);
% Evaluate the initial value wrt the complete data likelihood     
l_0 = OU_LIKE(lambda0,mu0,sigma20,u,v,X);                   
if abs(l_0) <= delta  
           %% check to see if initial guess satisfies
  return;                       %% convergence criterion.
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                    %%
%% MAIN ROUTINE                                                                                                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_iter=2000;
iter=0;
error = 1000;
while (error > delta && iter < max_iter)
  l_0 = OU_LIKE(lambda0,mu0,sigma20,u,v,X);
% update parameters lambda and mu using Newton's method 
  theta_update = [lambda0; mu0]-OU_LIKE_HESS_UPDATE(lambda0,mu0,sigma20,u,v,X);
  lambda0 = OU_SUB_reset_bound(theta_update(1));
  mu0= OU_SUB_reset_bound(theta_update(2));
% update sigma2
  sigma20 =  OU_SIG_LIKE(lambda0,mu0,u,v,X);
  sigma2_update = sigma20;
% measure error in the likelihood function 
  error= abs(l_0 - OU_LIKE(lambda0,mu0,sigma20,u,v,X));
% update iteration
  iter = iter +1;
% print stuff
fprintf('\n Newton iteration  = %d, delta = %d,\n', iter,error)
fprintf('\n lambda = %d, mu = %d, sigma2 = %d, \n', lambda0, mu0,sigma20)
theta = [theta_update(1);theta_update(2);sigma2_update];
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sub function to make sure the estimates don't go nuts.
function theta_reset = OU_SUB_reset_bound(x)
	if (x <= 1 || x > 100)
    x = 2+randi(5);
fprintf('\n Parameter values reset to default,value = %d, \n',x)
else
end
theta_reset = x;                                 
return;
