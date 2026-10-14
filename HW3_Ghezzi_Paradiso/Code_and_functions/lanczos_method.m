function [T, Q, lambda, v, eigenvalues, eigenvectors] = lanczos_method(A, q1, tol)
% function implementation of the Lanczos method for eigenvalues, returning
% a tridiagonal matrix T similar to A
% input: 
% A : nxn symmetric positive definite matrix of which we want to compute eigenvalues 
% q1 : n dimensional starting vector of the sequence qk
% tol : tolerance for stopping criterion
%
% output:
% T: kxk tridiagonal matrix similar to A through Q
% Q : nxk matrix having as columns the vectors qk st T = Q'AQ
% eigenvalues : n dimensional vector containing the eigenvalues of A
% eigenvectors : nxn matrix of eigenvectors of A

n = size(A, 1);
alpha_n = zeros(n,1); %vector of all alphak (if we iterate until k=n)
beta_n = zeros(n-1, 1); % vector of betak (if we iterate until k=n)
Q_n = zeros(n,n); % initialize Q
q1 = q1 / norm(q1,2); %normalize q1
Q_n(:,1) = q1; % 1st column of Q = q1
qk = q1; % initialize sequence qk
alphak = qk' * A * qk; % initialize alpha1
alpha_n(1) = alphak; % 1st element of alpha_n = alpha1
rk = A* qk - alphak*qk; % initialize rk
betak = norm(rk, 2); % initialize beta1
beta_n(1) = betak; % 1st element of beeta_n = betak

k=2;
while k<=n && betak>tol
    qkm1= qk;
    qk = rk/betak;
    % save vector qk as kth column of Q
    Q_n(:, k) = qk;
    % compute rk
    rk = A*qk - betak*qkm1;
    % compute and save alphak
    alphak = qk' * A * qk;
    alpha_n(k) = alphak;
    %rk = rk - alphak*qk; 
    % with full reorthogonalization
    rk = rk - Q_n(:,1:k)*(Q_n(:, 1:k)'*rk);
    % compute and save betak
    betak = norm(rk, 2);
    beta_n(k) = betak;
    k = k+1;
end
k=k-1;
alpha = alpha_n(1:k);
beta = beta_n(1:k-1);
Q = Q_n(:, 1:k);
T = diag(alpha) + diag(beta,1) + diag(beta,-1);

% since the method returns a tridiagonal matrix similar to A, we use it to
% compute its eigenvalues and eigenvectors is a faster way
% compute greater and smallest eigenvalue with power and inverse power
maxIter = 10000;
relTol = 1e-12;


lambda=0;
v= zeros(n,1); 
eigenvalues=zeros(n,1);
eigenvectors=zeros(n,n);

% if strcmp(which, 'first')
%     %lambda1 greater eigenvalue of A and T
%     [lambda,v1_T, ~] = power_method(T, v0, maxIter, relTol);
%     % since v1_T corresponding eigenvector for T, let us retreive the one for A
%     % using Q
%     v = Q*v1_T;
% elseif strcmp(which, 'last')
%     % lambdan smallest eigenvalue of A and T
%     [lambda,vn_T, ~] = inverse_power_method(T, v0, maxIter, relTol);
%     % since v1_T corresponding eigenvector for T, let us retreive the one for A
%     % using Q
%     v = Q*vn_T;

% compute eigenvalues of T using deflated PM
eigenvalues= eig(T);
% check if they are different 
if length(eigenvalues)==length(unique(eigenvalues))
    % since Ak symmetric use shifting + PM for eigenvectors
    [eigenvalues, eigenvectors_T] = compute_eigenpairs(T, maxIter, relTol);
else % otherwise I use the matlab eig to find eigenvectors
    [eigenvectors_T, ~] = eig(T);
end
eigenvectors = Q*eigenvectors_T;

end