function [eigenvalues, eigenvectors, Ak , k] = QR_method(A,both, relTol, maxIter, improve)
% function implementation of the QR method for eigenvalues
%
% input:
% A : nxn matrix of which we want to compute eigenvalues
% both : string specifying whether to compute also the eigenvectors; if
% 'yes' the function computes both eigenvalues and eigenvectors, if 'no'
% it only computes eigenvalues
% relTol : tolerance for stopping criterion
% maxIter : maximum number of iterations
% improve : string specifying whether to apply the algorthm to an upper
% hessenberg matrix similar to A (if so, improve = 'up_hess')
%
% output:
% eigenvalues : n dimensional vector containing the eigenvaues of A in
% decreasing order
% eigenvectors : nxn matrix containing the eigenvectors of A
% Ak : matrix created by the QR iterations
% k : number of iterations performed 

%%% setting default value for both
if nargin < 5
    both = 'yes';
end

n = size(A,1);
Ak=A; % initializing sequence
if strcmp(improve, 'up_hess')
    % compute an upper hessemberg matrix similar to A
    Ak = similar_hessemberg(A);
end

Sk=diag(Ak, -1); % subdiagonal of Ak
Skm1=Sk+1; % subdiagonal of Ak-1
prodQk = eye(n,n); % product of matrices Qk to compute eigenvectors
k=1;
% stopping criterion : small difference of two subsequent subdiagonals
while k<maxIter && norm(Sk-Skm1, 2)>=relTol*norm(Sk, 2)
    % updating Sk-1
    Skm1 = Sk;
    % computing QR decomposition of Ak
    [Qk, Rk] = qr(Ak);
    % updating Ak+1
    Ak = Rk*Qk;
    % updating the the product of Qk
    prodQk = prodQk*Qk;
    % updating Sk
    Sk = diag(Ak,-1);
    k = k+1;
end
k=k-1;
% check if Ak is upper triangular (otherwise it's upper Hessemberg)
% since it is upper triangular with some tolerance we check if diags are
% almost zero
triu =  1;
for i=1:n-1
    % check if any element of the lower triangular part of A is greater than relTol
    if any(abs(diag(Ak, -i)) > relTol)
        triu = 0;
    end
end
if triu == 1
    eigenvalues = diag(Ak); % returning eigenvalues of A 
else 
    % if Ak not upper triangular, it's upper Hessemberg having non zero
    % subdiagonal elements correponding to the first repeated eigenvalue (in abs)
    % we have to compute the eigenvalues of the 2x2 blocks of the upper
    % hessemberg part

    % Find the indices of non-zero subdiagonal elements (indices of columns of Ak)
    nonzero_indices = find( abs(diag(Ak, -1)) > relTol);
    eigenvalues = diag(Ak);

    for i=1:length(nonzero_indices)
        col_index = nonzero_indices(i);
        block2x2 = Ak(col_index:col_index+1, col_index:col_index+1);
        % computing eigenvalues of the 2x2 block
        delta = (block2x2(1,1) + block2x2(2,2))^2  - 4*(block2x2(1,1)*block2x2(2,2) - block2x2(1,2)*block2x2(2,1));
        lambda1 = 0.5*(block2x2(1,1) + block2x2(2,2) - sqrt(delta));
        lambda2 = 0.5*(block2x2(1,1) + block2x2(2,2) + sqrt(delta));
        eigenvalues(col_index) = lambda1;
        eigenvalues(col_index+1) = lambda2;
    end

end
eigenvectors = zeros(n,n);
if strcmp(both, 'yes')
    % if eigenvalues are different I can apply compute_m_greatest_eigenvalues (mix of PM
    % and deflation) to Ak
    if length(eigenvalues)==length(unique(eigenvalues))
        % if Ak symmetric use shifting + PM
        if abs(Ak-Ak')<= relTol
            [~, eigenvectors_Ak] = compute_eigenpairs(Ak, maxIter, relTol);
        % otherwise use deflation + PM
        else 
            [~, eigenvectors_Ak] = compute_m_greatest_eigenvalues(Ak,n,'yes', maxIter, relTol);
        end
    else %otherwise I use the matlab eig
        [eigenvectors_Ak, ~] = eig(Ak);
    end
    eigenvectors = prodQk*eigenvectors_Ak;
end


end