function [u_new, pi_new, lambda_new, err, iter] = UzawaSignoriniSolver(p, f, K, r, g, ibcs)
%--------------------------------------------------------------------------
% UzawaSignoriniSolver 
%
% Solves a variational inequality with Signorini boundary conditions 
% using a regularized Uzawa iterative method.
%
% Syntax:
%   [u_new, pi_new, lambda_new, err, iter] = UzawaSignoriniSolver(p, f, K, r, g, ibcs)
%
% Inputs:
%   p      - Node coordinates (np x 2)
%   f      - Right-hand side vector (np x 1)
%   K      - System stiffness matrix (sparse np x np)
%   r      - Regularization parameter (scalar)
%   g      - Gap function at boundary nodes (vector)
%   ibcs   - Indices of boundary nodes where Signorini condition applies
%
% Outputs:
%   u_new      - Solution vector
%   pi_new     - Projected boundary values (Signorini projection)
%   lambda_new - Lagrange multipliers (reaction forces at boundary nodes)
%   err        - Final relative error
%   iter       - Number of iterations performed
%
%--------------------------------------------------------------------------
% (c) O.Elamraoui, E-H.Essoufi, A.Zafrar 2025. Created for the fractional temporal SIER model
%     with Signorini boundary conditions. All rights reserved.
%--------------------------------------------------------------------------

    %% Parameters
    tol = 1e-8;
    max_iter = 200;
    err = 1;
    iter = 1;
    
    en = length(ibcs);
    enum = [ibcs(1:end-1), ibcs(2:end)];
    
    lambda_prev = 0.01 * ones(en, 1);
    pi_prev = g - 0.1 * rand(en, 1);
    
    %% Initial step
    SF = kpde2drhsn(p, enum, r * pi_prev - lambda_prev) + f;
    u_prev = K \ SF;
    
    pi_prev = u_prev(ibcs) + (1 / r) * (lambda_prev - max(0, lambda_prev + r * (u_prev(ibcs) - g)));
    lambda_prev = lambda_prev + r * (u_prev(ibcs) - pi_prev);

    %% Iterative Uzawa loop
    while err > tol && iter < max_iter
        SF = kpde2drhsn(p, enum, r * pi_prev - lambda_prev) + f;
        u_new = K \ SF;
        
        pi_new = u_new(ibcs) + (1 / r) * (lambda_prev - max(0, lambda_prev + r * (u_new(ibcs) - g)));
        lambda_new = lambda_prev + r * (u_new(ibcs) - pi_new);
        
        err = (norm(u_new - u_prev)^2 + norm(pi_new - pi_prev)^2) / ...
              (norm(u_new)^2 + norm(pi_new)^2);
        
        % Update for next iteration
        u_prev = u_new;
        pi_prev = pi_new;
        lambda_prev = lambda_new;
        
        iter = iter + 1;
    end

end
