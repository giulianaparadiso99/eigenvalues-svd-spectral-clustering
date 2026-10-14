function [lambda_1,v1, lambda_seq] = power_method(L, v0, maxIter, relTol)
% function implementation of the power method to compute the
% greatest eigenvalue lambda1 of the matrix L and its
% corresponding eigenvector v1
%
% input:
% L : nxn matrix of which we want to compute eigenvalues 
% NB: L must be a squared symmetric invertible matrix with greatest
% eigenvalue strictly greater than the other eigenvalues
% v0 : n dimensional column vector, initial approximation of the eigenvector
% maxIter : maximum number of iterations after which the algorithm stops
% relTol : relative tolerance used as stopping criterion for the algorithm
%
% output:
% lambda1 : greatest eigenvalue of L
% v1 : n dimensional column vector, eigenvector associated to the greatest
% eigenvalue lambda1
% lambda_seq : sequence of eigenvalues converging to lambda_1

% initializations
v0 = v0 / norm(v0, 2);
vk = v0;
k = 1;
lambda_k = 1e4;
lambda_km1 = inf;
lambda_sequence = zeros(maxIter, 1);

% loop
while k <= maxIter && abs(lambda_k - lambda_km1) >= relTol * abs(lambda_k)
    % Update mu_km1 after evaluating the stopping criterion
    lambda_km1 = lambda_k;

    % Perform one step of the power method
    v_tilde = L * vk;  
    lambda_k = vk' * v_tilde;  % Compute Rayleigh quotient
    vk = v_tilde / norm(v_tilde, 2);  % Normalize eigenvector
    lambda_sequence(k)=lambda_k;

    % Increment iteration counter
    k = k + 1;
end
lambda_seq = lambda_sequence(1:k-1);
lambda_1 = lambda_k;
v1 = vk;
end