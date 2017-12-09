%% This code implements the binomial approximation of the quadratic term in the likelihood function for OU process 
function log_like=OU_BINAPP2(X_miss,A,B,u,v,lambda,mu,sigma2)

% Initialize the terms:


N=size(X_miss,1);
M=size(X_miss,2);

Q      = NaN(M,1);
n      = NaN(M,1);
det_j  = NaN(M,1);
v1     = find(~isnan(X_miss(:,1)));
Xo_1   = X_miss(v1,1);
B11    = B(v1,v1);
b11    = [0,v(v1(2:end))-v(v1(1:end-1))];
b11=b11';
b11    = exp(-mu.*b11);
B11_inv= OU_COVINV(b11);

% Define the terms that go into the likelihood function for the first column:
% length of observed data
n(1)=length(Xo_1);
% determinant
det_j(1)=det(B11);
% quadratic form:
Q(1) = (Xo_1')*B11_inv*(Xo_1);
for j = 2:M
    % Define the observed data vector
    vj         = find(~isnan(X_miss(:,j)));
    vj1        = find(~isnan(X_miss(:,j-1)));
    Xo_j       = X_miss(vj,j);
    Xo_j1      = X_miss(vj1,j-1);
    n(j)       = length(vj);
    % define the sub matrices
    B_jj       = B(vj,vj);
    B_jj1      = B(vj,vj1);
    B_j1j      = B_jj1';
    
    % inverse of the matrices 
    % vector of site distances:
    b_jj       = [0,abs(v(vj(2:end))-v(vj(1:end-1)))];
    b_jj       = b_jj';
    B_jj_inv   = OU_COVINV(exp(-mu.*b_jj));
    b_j1j1     = [0,abs(v(vj1(2:end))-v(vj1(1:end-1)))];
    b_j1j1=b_j1j1';
    B_j1j1_inv = OU_COVINV(exp(-mu.*b_j1j1));
    
    % LU decomposition of each sub matrix B_jj and B_{j-1,j-1} which will
    % be used for approximation of (S_j)^{-1}
    Lj         = OU_SQRTM(B_jj_inv);
    Lj1        = OU_SQRTM(B_j1j1_inv);
    eta_j      = abs(u(j)-u(j-1));
    % conditional mean
    mo_j       = exp(-lambda*eta_j)*B_jj1*B_j1j1_inv*Xo_j1;
   % err = Xo_j-mo_j
   % conditional covariance
   % S_jj1      = exp(-2*lambda*eta_j)*B_jj1*(B_j1j1_inv)*B_j1j;
    S_j        = B_jj - exp(-2*lambda*eta_j)*B_jj1*(B_j1j1_inv)*B_j1j; 
    C_j        = Lj'*B_jj1*Lj1;
   % size_Cj=size(C_j)
    
    I_j        = eye(n(j),n(j));
   
    C_j_square=C_j*C_j';
    Approx_term = (I_j+exp(-2*lambda*eta_j)*C_j_square +...
        exp(-4*lambda*eta_j)*(C_j_square)^2+...
        exp(-6*lambda*eta_j)*(C_j_square)^3);
    %S_j_approx = ((Lj')\Approx_term)/(Lj);
    % Compute the determinant of the covariance matrix
    %det_j(j) = det(Lj)*det(Lj')*det(Approx_term);
    det_j(j)=det(S_j);
   % Q(j) = (Xo_j-mo_j)'*(S_j_approx)*(Xo_j-mo_j);
    Q(j) = ((Xo_j-mo_j)'/(S_j))*(Xo_j-mo_j);
    
    
    
end 
log_like=sum(n)*log(2*pi*sigma2)+sum(log(det_j))+(1/sigma2)*sum(Q);