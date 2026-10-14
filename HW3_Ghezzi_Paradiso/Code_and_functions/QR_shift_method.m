function [eigenvalues, kseq] = QR_shift_method(A,maxIter, relTol)
% function implementation of the QR method with shifting
%
% input: 
% A : nxn matrix of which we want to compute eigenvalues
% maxIter : maximum number of iterations
% relTol : tolerance for stopping criterion
%
% output:
% eigenvalues : n dimensional vector containing the eigenvaues of A in
% decreasing order
% kseq : sequence of number of iterations performed for every shift

% compute an upper hessemberg matrix similar to A
A = similar_hessemberg(A);
n = size(A, 1);
eigenvalues = zeros(n,1);

%CHECK IF A POSDEF
posdef =  1;
try
    chol(A);
catch
    posdef = 0;
end

% initialize sequence of Ak
Ak=A;
k=1;
kseq=zeros(n,1);
% sequence of smaller matrices
for i=n:-1:2
    % new matrix Ak
    Ak=Ak(1:i,1:i);
    % last subdiagonal element
    muk_sub = Ak(i,i-1);
    % identity matrix for next iteration
    Inmk = eye(size(Ak));

    muk = Ak(i,i);
    if eigenvalues(i)==0 && posdef == 0
        eigenvalues(i)=muk;
    end
    % end while loop when the last subdiagonal element is close to 0
    while k <= maxIter && abs(muk_sub)>relTol
        % last subdiagonal element
        muk_sub = Ak(i, i-1);
        % last diagonal element
        muk = Ak(i,i);
        % QR decomposition of the shifted Ak
        [Qk, Rk] = qr(Ak - muk*Inmk);
        % updating Ak with shift
        Ak = Rk * Qk + muk*Inmk;
        k = k+1;
    end
    kseq(n+1-i)=k-1;
    % when we exit the while loop, if Ak(i,i-1) is almost 0 
    % then Ak(i-1,i-1) = lambdai eigenvalue of A
    % otherwise we compute the eigenvalues of the 2x2 block left
    muk_sub = Ak(i, i-1);
    if abs(muk_sub) > relTol && posdef==0
        block2x2 = Ak(i-1:i, i-1:i);
        % computing eigenvalues of the 2x2 block
        delta = (block2x2(1,1) + block2x2(2,2))^2  - 4*(block2x2(1,1)*block2x2(2,2) - block2x2(1,2)*block2x2(2,1));
        lambda1 = 0.5*(block2x2(1,1) + block2x2(2,2) - sqrt(delta));
        lambda2 = 0.5*(block2x2(1,1) + block2x2(2,2) + sqrt(delta));
        eigenvalues(i) = lambda1;
        eigenvalues(i-1) = lambda2;
    elseif posdef == 1
        eigenvalues(i)=muk;
    end
end
eigenvalues(1)=Ak(1,1);
end
