%----------------------------------------------------------------------------------------
% Program: Fractional Temporal SIER Model with Signorini Boundary Conditions
% Article: Switching cases for fractional time SEIR Model with memory and space diffusion 
% Authors: Omar Elamraoui, EL Hassan Essoufi, Abderrahim Zafrar
% Created by: Omar Elamraoui
% Description: This program simulates a fractional-order SIER epidemic model
%              using finite element methods with Signorini boundary conditions.
%
% Copyright: © O.Elamraoui, E-H.Essoufi, A.Zafrar, 2025. All rights reserved.
%----------------------------------------------------------------------------------------

%% Parameters
nx = 20; ny = 20;             % Grid size
r = 2;                      % Regularization parameter
g = 0.2;                    % Regularization parameter for UBR22
dt = 0.01;                  % Time step
Nt = 500;                   % Number of time steps
tt = 0:dt:Nt*dt;            % Time vector
% Model parameters
ba = 0.2;                  % Infection rate 
gm = 0.1;                  % Recovery rate
d = 0.25;                  % Natural death rate
di = 0.4;                  % Disease-related death rate
b = 0.525;                 % Birth rate
sig = 0.35;                % Incubation rate
alpha=0.5;                 % Fractional order
num_iter=0;
%% Mesh generation
[p, t, pbx, pby] = kpde2dumsh(0, 1, 0, 1, nx, ny);
np = size(p, 1);

% Boundary indices
ibcs = pbx(:, 1); 
en = length(ibcs);
enum = [ibcs(1:end-1), ibcs(2:end)];

%% Initial conditions
S0 = 0.7 * ones(np, 1);
E0 = zeros(np, 1);
I0 = zeros(np, 1);
R0 = zeros(np, 1);

center_x = 0.5;
center_y = 0.5;
radius = 0.4;

% Apply initial Gaussian-like infection distribution
for i = 1:np
    x = p(i,1);
    y = p(i,2);
    if (x - center_x)^2 + (y - center_y)^2 <= radius^2
        E0(i) = 0.5;
        I0(i) = 0.3;
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
S(:,1) = S0;
E(:,1) = E0;
I(:,1) = I0;
R(:,1) = R0;
%% Time loop
for n = 2:length(tt)
    bS = zeros(np,1); bE = zeros(np,1); bI = zeros(np,1); bR = zeros(np,1);
    
    if n == 2
        bS = -dt^(-alpha) * S0;
        bE = -dt^(-alpha) * E0;
        bI = -dt^(-alpha) * I0;
        bR = -dt^(-alpha) * R0;
    else
        w = zeros(1,n-1);
        for k = 1:n-1
            w(k) = (tt(n)-tt(k))^(1-alpha) - (tt(n)-tt(k+1))^(1-alpha);
        end
        
        c = zeros(1,n-1);
        c(1) = ((tt(n)-tt(n-2))^(1-alpha) - 2*(tt(n)-tt(n-1))^(1-alpha)) / dt;
        for k = 2:n-2
            c(k) = ((tt(n)-tt(n-k-1))^(1-alpha) - 2*(tt(n)-tt(n-k))^(1-alpha) + ...
                   (tt(n)-tt(n-k+1))^(1-alpha)) / dt;
        end
        c(n-1) = -((tt(n)-tt(1))^(1-alpha) - (tt(n)-tt(2))^(1-alpha)) / dt;
        
        for j = 1:n-1
            bS = bS + c(j) * S(:,n-j);
            bE = bE + c(j) * E(:,n-j) + w(j) * E(:,j);
            bI = bI + c(j) * I(:,n-j) + w(j) * I(:,j);
            bR = bR + c(j) * R(:,n-j);
        end
    end

    % Update fields
    FS = kpde2drhs(p, t, -gamma(2-alpha)*ba*S0.*I0 - bS + gamma(2-alpha)*(b - d*S0));
    S1 = K1 \ FS;
    
    FE = kpde2drhs(p, t, gamma(2-alpha)*ba*S0.*I0 - gamma(2-alpha)*(sig + d)*E0 - bE);
    E1 = K1 \ FE;
    
    FI = kpde2drhs(p, t, gamma(2-alpha)*sig*E0 - gamma(2-alpha)*(gm + d + di)*I0 - bI);
    [I1, pI1, lambdaI1, errI, iterI] = UzawaSignoriniSolver(p, FI, K, r, g, ibcs);
    num_iter = num_iter + iterI;

    FR = kpde2drhs(p, t, gamma(2-alpha)*gm*I0 - gamma(2-alpha)*d*R0 - bR);
    R1 = K1 \ FR;

    % Save new values
    S0 = S1; E0 = E1; I0 = I1; R0 = R1;
    S(:,n) = S0; E(:,n) = E0; I(:,n) = I0; R(:,n) = R0;
    opS(n-1) = sum(S0);opE(n-1) = sum(E0);
    opI(n-1) = sum(I0);opR(n-1) = sum(R0);
end 


