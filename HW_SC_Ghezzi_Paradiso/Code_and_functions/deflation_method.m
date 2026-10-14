function [P_bar, small_A, b, PAP] = deflation_method(A, v1)
% Function implementation of deflation method :
% Given lambda_1, v1 the dominant eigenvector of a matrix A and its corresponding eigenvector,
% deflation method computes the matrix A2 having as dominant eigenvalue
% mu_2, the second greatest eigenvalue of A
%
% input: 
% A : nxn matrix having as dominant eigenvalue lambda_1 
% v1 : eigenvector associated to the dominant eigenvalue lambda_1
%
% output:
% P_bar : nxn matrix st P_bar A P_bar' = [lambda_1, b' ; 0, small_A] block matrix having as
% 1st column [lambda_1, 0, ..., 0] and 1rst row [lambda_1, --b1'--] and the rest
% defines the (n-1)x(n-1) matrix small_A
% small_A : (n-1)x(n-1) matrix having as dominant eigenvalue the 2nd greatest
% eigenvalue of A
% b : first row of P_bar from 2nd element
% PAP : product P_bar A P_bar'

In = eye(size(A)) ;
e1 = zeros(size(A, 1), 1);
e1(1)=1;
v1pluse1 = v1 + e1;

P_bar = In - 2*(v1pluse1*v1pluse1') / sum((v1pluse1).^2);
PAP = P_bar * A * P_bar';
small_A = PAP(2:end, 2:end);
b = P_bar(1, 2:end);
b = b';

end