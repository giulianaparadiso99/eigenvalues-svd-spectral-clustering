function [A_Hess] = similar_hessemberg(A)
% function to compute, given a matrix A, a similar matrix upper Hessemberg
% input: 
% A : nxn matrix to make upper Hessemberg
% output:
% A_Hess : nxn matrix similar to A upper Hessemberg 
n=size(A, 1);
Ai = A;

for k=1:n-2
    Ti = eye(n,n);
    % column vector corresponding to the current column of A 
    % (from the subdiagonal element on) to make upper Hess
    ai = Ai(k+1:n, k);
    e1 = zeros(size(ai, 1), 1);
    e1(1)=norm(ai, 2);
    aipluse1 = ai + e1;

    % matrix Hi of dim (n-k)x(n-k) orthogonal
    Hi = eye(n-k, n-k) - 2*(aipluse1*aipluse1')/sum((aipluse1).^2);
    % matrix Ti nxn having I kxk in the upper block and Hi in the lower
    % block
    Ti(k+1:n, k+1:n) = Hi;
    % new matrix with row k upper hessemberg
    Ai = Ti * Ai * Ti';
end
A_Hess = Ai;
end