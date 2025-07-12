%----------------------------------------------------------------------------------------
% Program: Fractional Temporal SIER Model with Signorini Boundary Conditions
% Article: Switching cases for fractional time SEIR Model with memory and space diffusion 
% Authors: Omar Elamraoui, EL Hassan Essoufi, Abderrahim Zafrar
% Created by: Omar Elamraoui
% Description: This program Showing effect of  varing fractional-order alpha on
%                               SIER epidemic model
%
% Copyright: © O.Elamraoui, E-H.Essoufi, A.Zafrar, 2025. All rights reserved.
%----------------------------------------------------------------------------------------
%% Parameters
nx = 20; ny = 20;           % Number of points in each dimension
dt = 0.01; Nt = 5000;       % Time step and number of time steps
tt = 0:dt:Nt*dt;            % Time vector

r = 2; g = 0.2;             % Penalty and regularization parameters
ba = 0.2;                   % Infection rate
gm = 0.1;                   % Recovery rate
d = 0.25;                   % Death rate
di = 0.4;                   % Death rate due to illness
b = 0.525;                  % Birth rate
sig = 0.35;                 % Incubation rate

%% Mesh generation
[p, t, pbx, pby] = kpde2dumsh(0, 1, 0, 1, nx, ny);
np = size(p, 1);
nt = size(t, 1);

%% Boundary conditions
ibcs = pbx(:, 1);
en = length(ibcs);
enum = [ibcs(1:end-1), ibcs(2:end)];

%% Varying alpha values
alpha_values = [0.5, 0.75, 0.9, 1.0];
S_all = cell(length(alpha_values), 1);
E_all = cell(length(alpha_values), 1);
I_all = cell(length(alpha_values), 1);
R_all = cell(length(alpha_values), 1);

%% Loop over alpha values
for alpha_idx = 1:length(alpha_values)
    alpha = alpha_values(alpha_idx);

    %% Initial conditions (Gaussian region)
    center_x = 0.5; center_y = 0.5;
    radius = 0.4;

    S0 = 0.7 * ones(nx * ny, 1);
    E0 = zeros(nx * ny, 1);
    I0 = zeros(nx * ny, 1);
    R0 = zeros(nx * ny, 1);

    for i = 1:nx
        for j = 1:ny
            idx = (i-1)*ny + j;
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

    %% Initialize solution arrays
    S(:, 1) = S0;
    E(:, 1) = E0;
    I(:, 1) = I0;
    R(:, 1) = R0;

    %% Time stepping
    aver_iter = 0;

    for n = 2:length(tt)
        % Initialize RHS vectors
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
            w_coeffs = zeros(1, n-1);
            c_coeffs = zeros(1, n-1);

            for k = 1:n-1
                w_coeffs(k) = (tt(n) - tt(k))^(1 - alpha) - (tt(n) - tt(k + 1))^(1 - alpha);
            end

            c_coeffs(1) = ((tt(n) - tt(n-2))^(1 - alpha) - 2 * (tt(n) - tt(n-1))^(1 - alpha)) / dt;

            for k = 2:n-2
                c_coeffs(k) = ((tt(n) - tt(n - k - 1))^(1 - alpha) - ...
                              2 * (tt(n) - tt(n - k))^(1 - alpha) + ...
                              (tt(n) - tt(n - k + 1))^(1 - alpha)) / dt;
            end

            c_coeffs(n-1) = -((tt(n) - tt(1))^(1 - alpha) - (tt(n) - tt(2))^(1 - alpha)) / dt;

            for j = 1:n-1
                bS = bS + c_coeffs(j) * S(:, n - j);
                bE = bE + c_coeffs(j) * E(:, n - j) + w_coeffs(j) * E(:, j);
                bI = bI + c_coeffs(j) * I(:, n - j) + w_coeffs(j) * I(:, j);
                bR = bR + c_coeffs(j) * R(:, n - j);
            end
        end

        % Update S
        FS = kpde2drhs(p, t, -gamma(2 - alpha) * ba * S0 .* I0 - bS + gamma(2 - alpha) * (b - d * S0));
        S1 = K1 \ FS;

        % Update E
        FE = kpde2drhs(p, t, gamma(2 - alpha) * ba * S0 .* I0 - gamma(2 - alpha) * (sig + d) * E0 - bE);
        E1 = K1 \ FE;

        % Update I
        FI = kpde2drhs(p, t, gamma(2 - alpha) * sig * E0 - gamma(2 - alpha) * (gm + d + di) * I0 - bI);
        [I1, ~, ~, ~, iterI] = UzawaSignoriniSolver(p, FI, K, r, g, ibcs);
        aver_iter = aver_iter + iterI;

        % Update R
        FR = kpde2drhs(p, t, gamma(2 - alpha) * gm * I0 - gamma(2 - alpha) * d * R0 - bR);
        R1 = K1 \ FR;

        % Store new values
        S0 = S1; E0 = E1; I0 = I1; R0 = R1;
        S(:, n) = S0;
        E(:, n) = E0;
        I(:, n) = I0;
        R(:, n) = R0;
    end

    % Store results for current alpha
    S_all{alpha_idx} = S;
    E_all{alpha_idx} = E;
    I_all{alpha_idx} = I;
    R_all{alpha_idx} = R;
end

