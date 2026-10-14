function [A1] = shifting_method(A,lambda1, v1)
% function implementation of the shifting method for eigenvalues
%  
% input: 
% A : nxn symmetric matrix of which we know the 1st eigenpair and want to
% compute the 2nd
% lambda1 : dominant eigenvalue of A
% v1 : dominant eigenvector of A
% 
% output: 
% A1: nxn matrix A1 having as greatest eigenpair the second
% greatest eigenpair of A (lambda2, v2)
A1 = A - lambda1*(v1*v1');
end