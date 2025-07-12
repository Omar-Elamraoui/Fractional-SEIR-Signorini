function [S, E, I, R] = run_SEIR_simulation(p, t, ibcs, tt, alpha, ba, gm, d, di, b, sig, r, g, dt)
%--------------------------------------------------------------------------
% RUN_SEIR_SIMULATION
%
% Solves a fractional time-space SEIR model with Signorini boundary 
% conditions using a regularized Uzawa iterative method.
%
% INPUTS:
%   p, t   - mesh nodes and connectivity
%   ibcs   - indices of boundary nodes
%   tt     - time steps
%   alpha  - fractional order
%   ba     - transmission rate
%   gm     - recovery rate
%   d, di  - natural and disease-induced death rates
%   b      - birth rate
%   sig    - incubation rate
%   r, g   - penalty and regularization parameters
%   dt     - time step size
%
% OUTPUTS:
%   S, E, I, R - matrices of solutions at each time step
%
%--------------------------------------------------------------------------
% (c) O. Elamraoui, E-H. Essoufi, A. Zafrar 2025. 
%     All rights reserved.
%--------------------------------------------------------------------------
    np = size(p, 1);
    nx = length(unique(p(:,1)));
    ny = length(unique(p(:,2)));

    %% Initial Conditions (Gaussian)
    S0 = 0.7 * ones(nx * ny, 1);
    E0 = zeros(nx * ny, 1);
    I0 = zeros(nx * ny, 1);
    R0 = zeros(nx * ny, 1);
    center_x = 0.5;
    center_y = 0.5;
    radius = 0.4;

    for i = 1:nx
        for j = 1:ny
            idx = (i - 1) * ny + j;
            x = p(idx, 1);
            y = p(idx, 2);
            if (x - center_x)^2 + (y - center_y)^2 <= radius^2
                E0(idx) = 0.5;
                I0(idx) = 0.3;
            end
        end
    end

    %% Matrices
    A = kpde2dstf(p, t, gamma(2 - alpha));
    M = full(kpde2dmss(p, t, dt^(-alpha)));
    K1 = M + A;
    Ms = zeros(np);
    Ms(ibcs, ibcs) = M(ibcs, ibcs);
    K = K1 + r * Ms;

    %% Time Loop Initialization
    S(:, 1) = S0;
    E(:, 1) = E0;
    I(:, 1) = I0;
    R(:, 1) = R0;
    aver_iter = 0;

    %% Time Integration Loop
    for n = 2:length(tt)
        bS = zeros(np, 1);
        bE = zeros(np, 1);
        bI = zeros(np, 1);
        bR = zeros(np, 1);

        if n == 2
            bS = -dt^(-alpha) * S0;
            bE = -dt^(-alpha) * E0;
            bI = -dt^(-alpha) * I0;
            bR = -dt^(-alpha) * R0;
        else
            w = zeros(1, n - 1);
            c = zeros(1, n - 1);

            for k = 1:n-1
                w(k) = (tt(n) - tt(k))^(1 - alpha) - (tt(n) - tt(k + 1))^(1 - alpha);
            end

            c(1) = ((tt(n) - tt(n-2))^(1 - alpha) - 2 * (tt(n) - tt(n-1))^(1 - alpha)) / dt;
            for k = 2:n-2
                c(k) = ((tt(n) - tt(n-k-1))^(1 - alpha) - 2 * (tt(n) - tt(n-k))^(1 - alpha) + ...
                        (tt(n) - tt(n-k+1))^(1 - alpha)) / dt;
            end
            c(n-1) = -((tt(n) - tt(1))^(1 - alpha) - (tt(n) - tt(2))^(1 - alpha)) / dt;

            for j = 1:n-1
                bS = bS + c(j) * S(:, n-j);
                bE = bE + c(j) * E(:, n-j) + w(j) * E(:, j);
                bI = bI + c(j) * I(:, n-j) + w(j) * I(:, j);
                bR = bR + c(j) * R(:, n-j);
            end
        end

        %% Solve for S, E, I, R
        FS = kpde2drhs(p, t, -gamma(2 - alpha) * ba * S0 .* I0 - bS + gamma(2 - alpha) * (b - d * S0));
        S1 = K1 \ FS;

        FE = kpde2drhs(p, t, gamma(2 - alpha) * ba * S0 .* I0 - gamma(2 - alpha) * (sig + d) * E0 - bE);
        E1 = K1 \ FE;

        FI = kpde2drhs(p, t, gamma(2 - alpha) * sig * E0 - gamma(2 - alpha) * (gm + d + di) * I0 - bI);
        [I1, ~, ~, ~, iterI] = UzawaSignoriniSolver(p, FI, K, r, g, ibcs);
        aver_iter = aver_iter + iterI;

        FR = kpde2drhs(p, t, gamma(2 - alpha) * gm * I0 - gamma(2 - alpha) * d * R0 - bR);
        R1 = K1 \ FR;

        %% Update for next step
        S0 = S1;
        E0 = E1;
        I0 = I1;
        R0 = R1;

        S(:, n) = S0;
        E(:, n) = E0;
        I(:, n) = I0;
        R(:, n) = R0;
    end
end
