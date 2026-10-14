function [B, P, H] = bidiagonal(A)
% function implementation to turn a rectangular (tall, m>=n) mxn matrix into
% bidiagonal form B mxn
%
% input:
% A : mxn matrix to make bidiagonal (m>=n)
%
% output:
% B : mxn matrix bidiagonal form of A st PAH = B
% P : mxm matrix applied from the left to A
% H : nxn matrix applied from the right to A

[m, n] = size(A);
% check to see if A is a tall rectngular matrix (m>=n), otherwise we
% transpose it
if m<n
 error('not a tall matrix (m>=n), transpose it');
end
B = A;
Ptilde = eye(m,m);
Htilde = eye(n,n);
P = Ptilde;
H = Htilde;
% at every step row and column k of A are transfromed to keep only their
% diagonal and upper diagonal part
for k=1:n
 Ptilde = eye(m,m);
 Htilde = eye(n,n);
 % taking column k of A from row k to end
 a = B(k:m, k); % column vector
 e1a=zeros(size(a));
 e1a(1)=1;
 % defining householder vector from a
 a = a + sign(a(1))*norm(a, 2)*e1a;
 Immkm1 = eye(m-(k-1), m-(k-1)); 
 % householder matrix part of Ptilde
 Ptilde(k:m, k:m) = Immkm1 - 2*(a*a')/sum(a.^2);
 % updating the matrix B (making column k diagonal)
 B = Ptilde*B;
 % updating matrix P
 P = Ptilde*P;
 if k<=n-2
     % taking row k of A from column k+1 to end
     b = B(k, k+1:n); % row vector
     e1b= zeros(size(b));
     e1b(1)=1;
     % defining householder vector from b
     b = b + sign(b(1))*norm(b,2)*e1b;
     Inmk = eye(n-k, n-k);
     % householder matrix part of Htilde
     Htilde(k+1:n, k+1:n) = Inmk - 2*(b'*b)/sum(b.^2);
     % updating the matrix B (making it bidiagonal until row k)
     B = B*Htilde;
     % updating matrix P
     H = H*Htilde; 
 end
end

end