function sens_indices = sensitivity_index_all(param_name, p_val)
%--------------------------------------------------------------------------
% SENSITIVITY_INDEX_ALL
%
% Computes sensitivity indices for a fractional time-space SEIR model
% with Signorini boundary conditions.
%
% INPUTS:
%   param_name - string, name of the parameter to perturb (e.g., 'beta', 'gamma')
%   p_val      - double, nominal value of the parameter
%
% OUTPUT:
%   sens_indices - struct containing normalized sensitivity indices for
%                  S, E, I, R at final simulation time
%
%--------------------------------------------------------------------------
% (c) O. Elamraoui, E-H. Essoufi, A. Zafrar 2025.
%     Created for the fractional temporal SEIR model with Signorini
%     boundary conditions. All rights reserved.
%--------------------------------------------------------------------------

    % Define perturbation (1% of parameter value)
    delta = 0.01 * p_val;

    % Baseline simulation
    [S0, E0, I0, R0] = analyze_param_impact(param_name, p_val);
    S_base = sum(S0{1}(:, end));
    E_base = sum(E0{1}(:, end));
    I_base = sum(I0{1}(:, end));
    R_base = sum(R0{1}(:, end));

    % Perturbed simulation
    [S1, E1, I1, R1] = analyze_param_impact(param_name, p_val + delta);
    S_pert = sum(S1{1}(:, end));
    E_pert = sum(E1{1}(:, end));
    I_pert = sum(I1{1}(:, end));
    R_pert = sum(R1{1}(:, end));

    % Finite difference derivatives
    dS_dp = (S_pert - S_base) / delta;
    dE_dp = (E_pert - E_base) / delta;
    dI_dp = (I_pert - I_base) / delta;
    dR_dp = (R_pert - R_base) / delta;

    % Compute normalized sensitivity indices
    sens_indices.S = dS_dp * (p_val / S_base);
    sens_indices.E = dE_dp * (p_val / E_base);
    sens_indices.I = dI_dp * (p_val / I_base);
    sens_indices.R = dR_dp * (p_val / R_base);

end
