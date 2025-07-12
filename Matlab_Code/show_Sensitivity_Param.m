%----------------------------------------------------------------------------------------
% Program: Parameter Sensitivity Analysis for Fractional Temporal SEIR Model
% Article: Switching cases for fractional time SEIR Model with memory and space diffusion 
% Authors: Omar Elamraoui, EL Hassan Essoufi, Abderrahim Zafrar
% Created by: Omar Elamraoui
% Description: This program performs sensitivity analysis on the fractional-order SEIR
%              model with Signorini boundary conditions. It varies a specified parameter,
%              runs simulations, and generates plots of population dynamics over time.
%
% Copyright: © O.Elamraoui, E-H.Essoufi, A.Zafrar, 2025. All rights reserved.
%----------------------------------------------------------------------------------------


clear; clc; close all;
%% Parameter configuration
param_name = '\beta';
param_values = [0.1, 0.2,0.3,0.4];

%% Run simulations
[S_all, E_all, I_all, R_all] = analyze_param_impact(param_name, param_values);

%% Time vector (must match analyze_param_impact)
dt = 0.01; 
Nt = 5000; 
tt = 0:dt:Nt*dt;

%% Plot settings
fontSize = 12;
lineWidth = 2;
colors = lines(length(param_values));

component_names = {'Susceptible', 'Exposed', 'Infected', 'Recovered'};
component_data = {S_all, E_all, I_all, R_all};

%% Plot results
for k = 1:4
    figure(k); clf; hold on;
    
    for i = 1:length(param_values)
        % Sum over spatial nodes to get total population
        total_comp = sum(component_data{k}{i}, 1);
        % Ensure correct time mapping (may already match tt)
        t_plot = linspace(tt(1), tt(end), length(total_comp));
        
        plot(t_plot, total_comp, 'LineWidth', lineWidth, 'Color', colors(i, :));
    end
    
    xlabel('Time', 'FontSize', fontSize);
    ylabel(['Number of ', component_names{k}], 'FontSize', fontSize);
    
    title(sprintf('%s Population for Varying $%s$', ...
        component_names{k}, param_name), ...
        'Interpreter', 'latex', 'FontSize', fontSize);
    
    legend(arrayfun(@(v) sprintf('$%s = %.2f$', param_name, v), param_values, ...
        'UniformOutput', false), ...
        'Location', 'best', 'Interpreter', 'latex');
    
    grid on;
    set(gca, 'FontSize', fontSize);
end

%% Map param_name to save identifier
switch param_name
    case {'\beta', 'ba'},    name = 'ba';
    case {'\gamma', 'gm'},   name = 'gm';
    case {'\sigma', 'sig'},  name = 'sig';
    case {'\alpha', 'alpha'},name = 'alpha';
    otherwise
        error('Unknown parameter name.');
end

%% Save figures
for k = 1:4
    fig_base = sprintf('%s_%s', name, component_names{k}(1));
    saveas(figure(k), [fig_base, '.jpg']);
    saveas(figure(k), [fig_base, '.eps']);
end
