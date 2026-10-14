function [lambda_n,vn, mu] = inverse_power_method(L, v0, maxIter, relTol )
% function implementation of the inverse power method to compute the
% smallest eigenvalue lambdan of the matrix L and its
% corresponding eigenvector vn
%
% input:
% L : nxn matrix of which we want to compute the smallest eigenvalue
% NB: L must be a squared symmetric invertible matrix with smallest
% eigenvalue strictly smaller than the other eigenvalues
% v0 : n dimensional column vector, initial approximation of the eigenvector
% maxIter : maximum number of iterations after which the algorithm stops
% relTol : relative tolerance used as stopping criterion for the algorithm
%
% output:
% lambdan : smallest eigenvalue of L
% vn : n dimensional column vector, eigenvector associated to the smallest eigenvalue lambdan
% mu : largest eigenvalue of the inverse of L L^(-1)

% initializations
v0 = v0 / norm(v0, 2);
vk = v0;
k = 0;
mu_k = 1e4;
mu_km1 = inf;

% loop
while k <= maxIter && abs(mu_k - mu_km1) >= relTol * abs(mu_k)
    % Update mu_km1 after evaluating the stopping criterion
    mu_km1 = mu_k;

    % Perform one step of the inverse power method
    v_tilde = L \ vk;  % Solve L * v_tilde = vk
    mu_k = vk' * v_tilde;  % Compute Rayleigh quotient
    vk = v_tilde / norm(v_tilde, 2);  % Normalize eigenvector

    % Increment iteration counter
    k = k + 1;
end
lambda_n = 1/mu_k;
mu=mu_k;
vn = vk;
k=k-1;
end