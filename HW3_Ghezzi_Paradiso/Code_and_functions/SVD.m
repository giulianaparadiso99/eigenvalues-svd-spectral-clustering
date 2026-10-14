function [U, Sigma, V, Q, d] = SVD(A, my_eigen)
% function to compute the SVD decomposition of a rectangular mxn (m>=n) matrix A
% A = USigmaV'
% the function can also be used to compute the eigenvalues and eigenvectors
% of A'A
%
% input : 
% A : mxn (m>=n) matrix of which we compute SVD
% my_eigen : string in {'QR','PM', else} telling whether to use: 'QR' the QR 
% method for eigenvalues, 'PM' the self implemented compute_eigenpairs 
% (for matrix with different eigenvalues), or  the matlab eig() to
% compute eigenvalues and eigenvectors of B'B (default is else)
%
% output:
% U : mxm orthogonal matrix of left singular values of A
% Sigma : mxn diagonal matrix having as diagonal elements the singular
% values of A in decreasing order
% V : nxn orthogonal matrix of right singular vectors of A
% Q : nxn matrix of eigenvectors of A'A
% d : n dimensional vector of eigenvalues of A'A 

%%% setting default value for both
if nargin < 2
    my_eigen = 'no';
end

[m, n] = size(A);
% check to see if A is a tall rectngular matrix (m>=n), otherwise we
% transpose it
if m<n
 error('not a tall matrix (m>=n), transpose it');
end

% turning A into bidiagonal form
[B, ~, H] = bidiagonal(A);

% computing eigenvectors and eigenvalues of B'B (also eigenvalues of A'A)
if strcmp(my_eigen, 'QR')
    [d, Q_tilde, ~ , ~] = QR_method(B'*B,'yes', 1e-10, 10000, 'up_hess');
elseif strcmp(my_eigen, 'PM')
    [d, Q_tilde] = compute_eigenpairs(B'*B, 10000, 1e-10);
else
    [Q_tilde, D] = eig(B'*B);
    d=diag(D);
end


% retrieving eigenvectors of A'A (being A'A = H(B'B)H' similar through H)
Q = H*Q_tilde;

C = A*Q;
% permuted QR decomposition of C st CP = UR
[U,R,P] = qr(C);
V = Q*P;
Sigma = abs(R);
end