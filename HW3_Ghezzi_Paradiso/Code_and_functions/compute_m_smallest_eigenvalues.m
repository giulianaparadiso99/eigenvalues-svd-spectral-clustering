function [eigenvalues_A, eigenvectors_A, mu_seq] = compute_m_smallest_eigenvalues(A,m,both, maxIter, relTol)
% function to compute the m smallest eigenvalues of the nxn matrix A
% alternating between inverse power method (IPM) and deflation method; in
% particular, starting from the original matrix A, it
% 1) computes the smallest eigenvalue (and corr eigenvector) of the nxn
% matrix  A using IPM
% 3) iteratively (for i=2:m) creates a smaller (n-(i-1) x n-(i-1)) matrix whose smallest
% eigenvalue is the ith smallest eigenvalue of A
% 4) applies IPM to find the smallest eigenvalue (and corr eigenvector) of
% the smaller matrix (the eigenvalue is also the ith smaller of A)
% 5) retrieves the eigenvector of A associated the just found ith
% eigenvalue
%
% input:
% A : nxn matrix of which we want to compute eigenvalues
% m : number of smallest eigenvalues of A to compute
% both : string specifying whether to compute also the eigenvectors; if
% 'yes' the function computes both eigenvalues and eigenvectors, if 'no'
% it only computes eigenvalues
% maxIter : maximum number of iterations after which IPM stops
% relTol : relative tolerance used as stopping criterion for IPM
%
% output:
% eigenvalues_A : m dimensional vector containing the m smallest
% eigenvalues of A
% eigenvectors_A : nxm matrix having as columns the eigenvectors
% corresponding to the m smallest eigenvalues
% mu_seq : eigenvalues of A^(-1)

%%% setting default value for both
if nargin < 5
    both = 'yes';
end

%%% initializing return vector of eigenvalues and matrix of eigenvectors
eigenvalues_A = zeros(m, 1); %column vector to be filled with the m smallest eigenvalues of A in increasing order
eigenvectors_A = zeros(size(A, 1), m); %matrix nxm to be filled, for every column, with one eigenvector of A

%%% setting parameters (could be done by user)
maxIter = maxIter; %3000; 
relTol = relTol; %1e-9;
v0 = ones(size(A, 1), 1); %starting vector for the IPM

%%% compute the smallest eigenvalue of A and corresponding eigenvector using IPM 
[lambda_n, vn, ~] = inverse_power_method(A, v0, maxIter, relTol);
eigenvalues_A(1) = lambda_n;
eigenvectors_A(:, 1) = vn;

%%% successive eigenvalues and eigenvector computation

% store last computed eigenvector of the reduced size Ai;
v_last_smallA =vn;

% store last computed small_AI
small_Ai_last = A;

% initialize dictionaries containing the already computed matrices and
% vectors, necessary to make the computations
P_dict = containers.Map('KeyType', 'int32', 'ValueType', 'any');
P_bar_dict = containers.Map('KeyType', 'int32', 'ValueType', 'any');
bi_dict = containers.Map('KeyType', 'int32', 'ValueType', 'any');

% define the "eigenvalue to fill"
eig_to_fill = zeros(size(A,1), 1);

% for loop to find eigenvalues 
for i=2:m
    %deflation method to find: small_Ai (having as smallest eigenvalue the
    %next eigenvalue we look for) , Pi_bar = I-2(xi +e1i)(xi+e1i)' / || xi + e1i ||^2
    % nb: xi = v_last_smallA last computed eigenvalue of small_Ai (same dimension)
    % PAPi = [lasteigenvalueA, bi; 0, small_Ai] (block matrix)
    [Pi_bar, small_Ai, bi, ~] = deflation_method(small_Ai_last, v_last_smallA);

    %save bi
    bi_dict(i-1) = bi;

    % save Pi_bar
    P_bar_dict(i-1) = Pi_bar; %at i=2 (looking for 2nd smallest eigenvalue) I find P1, b1
    if i==2
        % at 1st iteration Pi_bar = Pi
        P_dict(i-1) = Pi_bar;
    else
        % compute Pi (nxn) and save it
        Pi = eye(size(A,1), size(A, 1));
        Pi(i-1:end, i-1:end) = Pi_bar;
        P_dict(i-1) = Pi;
    end

    % IPM to compute new smallest eigenvalue of small_Ai (eigenvalue of A) and its
    % corresponding eigenvector (not eigenvector of A being it of smaller size wrt A)
    v02 = ones(size(small_Ai, 1), 1);
    [eigenvalues_A(i), v_last_smallA, mu_seq] = inverse_power_method(small_Ai, v02, maxIter, relTol); %+update v_last_smallA
    % if both='no' eigenvectors are not computed, otherwise:
    if strcmp(both, 'yes')
        % for loop to compute the eigenvetor of A starting from v_last_smallA
        eig_to_fill(i:end)=v_last_smallA;
        mod_v_last_smallA = v_last_smallA;
        for j=1:i-1
            % compute the i-1 values to fill the eigenvector of Bi starting from
            % the last empty one
            if j==1
                %compute the last value
                alphaimj = - (bi_dict(i-j)' * mod_v_last_smallA) / (eigenvalues_A(i-j) - eigenvalues_A(i));
                %put it inside the filled eigenvector
                eig_to_fill(i-j)=alphaimj;
                %update small eigenvector
                mod_v_last_smallA = [alphaimj; mod_v_last_smallA]; 
                %NB consider preallocating a matrix of zeros to fill rather
                %than changing size
            else
                %compute the last empty value
                alphaimj = - (bi_dict(i-j)' * P_bar_dict(i-(j-1)) * mod_v_last_smallA ) / (eigenvalues_A(i-j) - eigenvalues_A(i));
                %put it inside the eigenvector that we are filling
                eig_to_fill(i-j)=alphaimj;
                %update small eigenvector
                mod_v_last_smallA = [alphaimj; mod_v_last_smallA]; 
            end
        end

        %now multiply the eigenvector of B1 times Pi st the eigenvector of A
        %corresponding to the ith eigenvalue is xi=P1*P2*...*Pi-1*eig_to_fill
        for k=1:i-1
            eig_to_fill = P_dict(i-k)*eig_to_fill;
        end
        eigenvectors_A(:,i) = eig_to_fill/norm(eig_to_fill, 2);

    end

    %update small_Ai_last for new iteration
    small_Ai_last =small_Ai;
end
if m~=2
    mu_seq = '';
end

end