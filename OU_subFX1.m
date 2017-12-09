%% this sub-function calculates the exact form of -2logf(Xo_1) 
% where Xo_1 is a set of available observations inf the first column of an OU-field

function [log_like1]=OU_subFX1(X_miss,u,v,A,B,lambda,mu,sigma2)
% identify available sites:
v1_ind              = find(~isnan(X_miss(:,1)));
% get data from available sites:
Xo_1            = X_miss(v1_ind,1);
% get length of samples: 
n1              = length(v1_ind);
v1              = v(v1_ind);
v1=v1';
% define covariance matrix for the samples on the fixed column:
B11             = B(v1_ind,v1_ind);
% create vector that defines the inverse of B11:
b_v1            = [0;abs(v1(2:end)-v1(1:end-1))];
% define inverse of B11:
B11_inv         = OU_COVINV(exp(-mu.*b_v1));
% define terms of the determinant of B11:
%e1=exp(-2*mu.*b_v1(2:end))
det_terms_beta  = ones(n1-1,1)-exp(-2*mu.*b_v1(2:end));
% define the terms in the quadratic form:
% initialize the denominator terms
q_term_denom = NaN(n1,1);
% define the border terms (first and last):
q_term_denom(1)=det_terms_beta(1);
q_term_denom(end)=det_terms_beta(end);
% define the rest:
q_term_denom(2:end-1)=(1./det_terms_beta(2:end)) + (1./det_terms_beta(1:end-1))-ones(n1-2,1);
% put the terms together : 
q_term1=(Xo_1.^2)./q_term_denom;
% define data interaction term: 
x_term_int = Xo_1(2:end).*Xo_1(1:end-1);
cov_term_int = -exp(-mu.*b_v1(2:end))./det_terms_beta;
q_term2 = x_term_int.*cov_term_int;

% PUT ALL THE SCHITZ TOGETHER:
log_like1 = n1*log(2*pi*sigma2)+sum(log(det_terms_beta))+(1/sigma2)*(sum(q_term1)+2*sum(q_term2));
log_like2 = n1*log(2*pi*sigma2)+log(det(B11))+(1/sigma2)*(Xo_1'*B11_inv*Xo_1);
% err1=prod(det_terms_beta)-det(B11)
% 
% err2=(sum(q_term1)-2*sum(q_term2))-(Xo_1'*B11_inv*Xo_1)
% d1=prod(det_terms_beta)
% d2=det(B11)
% f1=log(prod(det_terms_beta))
% f2=log(det(B11))