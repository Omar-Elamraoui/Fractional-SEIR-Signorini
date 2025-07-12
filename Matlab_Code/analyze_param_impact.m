function [S_all, E_all, I_all, R_all] = analyze_param_impact(param_name, param_values)
%--------------------------------------------------------------------------
% ANALYZE_PARAM_IMPACT
%
% Runs the SEIR simulation for different values of a specified parameter 
% to study its impact on the solution.
%
% INPUTS:
%   param_name   - name of the parameter (string) e.g. '\beta', '\gamma', '\alpha'
%   param_values - vector of values to test
%
% OUTPUTS:
%   S_all, E_all, I_all, R_all - cell arrays containing solution matrices 
%                                for S, E, I, R at each parameter value
%
%--------------------------------------------------------------------------
% (c) O. Elamraoui, E-H. Essoufi, A. Zafrar 2025. All rights reserved.
%--------------------------------------------------------------------------

    %% Time and mesh setup
    dt = 0.01; 
    Nt = 5000; 
    tt = 0:dt:Nt*dt;
    nx = 20; ny = 20;
    [p, t, pbx, ~] = kpde2dumsh(0, 1, 0, 1, nx, ny);
    ibcs = pbx(:, 1);

    %% Default model parameters
    ba = 0.2;    % transmission rate
    gm = 0.1;    % recovery rate
    d = 0.25;    % natural death rate
    di = 0.4;    % disease-induced death rate
    b = 0.525;   % birth rate
    sig = 0.35;  % incubation rate
    r = 2.0;     % penalty parameter
    g = 0.2;     % regularization parameter
    alpha = 1;   % fractional order

    %% Preallocate result storage
    num_vals = length(param_values);
    S_all = cell(num_vals, 1);
    E_all = cell(num_vals, 1);
    I_all = cell(num_vals, 1);
    R_all = cell(num_vals, 1);

    %% Loop over parameter values
    for i = 1:num_vals
        val = param_values(i);

        % Update specified parameter
        switch param_name
            case {'\beta', 'ba'},    ba = val;
            case {'\gamma', 'gm'},   gm = val;
            case {'\sigma', 'sig'},  sig = val;
            case {'\alpha', 'alpha'},alpha = val;
            otherwise
                error('Unknown parameter name: %s', param_name);
        end

        % Run simulation
        [S, E, I, R] = run_SEIR_simulation(p, t, ibcs, tt, alpha, ba, gm, d, di, b, sig, r, g, dt);

        % Store results
        S_all{i} = S;
        E_all{i} = E;
        I_all{i} = I;
        R_all{i} = R;
    end
end
