function [eigenvalues, eigenvectors] = compute_eigenpairs(A, maxIter, relTol)
% function to compute eigenvalues of a nxn symmetric matrix A st 
% |lambda1| > |lambda2| > ... > |lambdan| (all different, in increasing
% order) iteratively using power method and shifting
%
% input:
% A : nxn symmetric matrix of which we want to compute eigenvalues
% maxIter : maximum numer of iterations for power method
% relTol : tolerance for power method
%
% output:
% eigenvalues : n dimensional vector of eigenvalues of A in decreasing
% order
% eigenvectors : nxn matrix having as columns the eigenvectors of A

n = size(A,1);
eigenvalues = zeros(n,1);
eigenvectors = zeros(n,n);
% compute the first eigenpair with power method
v0 = ones(n,1);
[lambda1,v1] = power_method(A, v0, maxIter, relTol);
v1 = v1 / norm(v1, 2);
eigenvalues(1) = lambda1;
eigenvectors(:, 1) = v1; 
% compute successive eigenpairs using shifting and then power method again
% on the resulting matrix
A_shift = A;
lambdalast = lambda1;
vlast = v1;
for i=2:n
    % apply shifting
    A_shift = A_shift - lambdalast*(vlast*vlast');
    % apply power method to A_shift
    [lambdalast, vlast] = power_method(A_shift, v0, maxIter, relTol);
    vlast  = vlast / norm(vlast, 2);
    eigenvalues(i) = lambdalast;
    eigenvectors(:,i) = vlast;
end

end