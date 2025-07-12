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
clear all; clc; close all;

%% Run SEIR model with varying alpha
SEIR_var_alpha_fun
%% Set common plot properties
lineWidth = 2;
fontSize = 12;
alpha_labels = {'\alpha = 0.5', '\alpha = 0.75', '\alpha = 0.9', '\alpha = 1.0'};
colors = {'r', 'b', 'g', 'c'};

%% Global SEIR dynamics for α = 0.5
f5 = figure(5); hold on;
plot(tt, sum(S_all{1}, 1), 'r', 'LineWidth', lineWidth);
plot(tt, sum(E_all{1}, 1), 'm', 'LineWidth', lineWidth);
plot(tt, sum(I_all{1}, 1), 'b', 'LineWidth', lineWidth);
plot(tt, sum(R_all{1}, 1), 'g', 'LineWidth', lineWidth);
xlabel('Time (days)', 'FontSize', fontSize);
ylabel('Number of individuals', 'FontSize', fontSize);
title('SEIR Model Dynamics Over Time for \alpha = 0.5', 'FontSize', fontSize);
legend('S: Susceptible', 'E: Exposed', 'I: Infected', 'R: Recovered', 'Location', 'best');
grid on;
set(gca, 'FontSize', fontSize);

%% Susceptible population for different α values
f1 = figure(1); hold on;
for k = 1:4
    plot(tt, sum(S_all{k}, 1), colors{k}, 'LineWidth', lineWidth);
end
xlabel('Time', 'FontSize', fontSize);
ylabel('Number of Susceptible Individuals', 'FontSize', fontSize);
title('Susceptible Population Over Time for Different \alpha Values', 'FontSize', fontSize);
legend(alpha_labels, 'Location', 'best');
grid on;
set(gca, 'FontSize', fontSize);

%% Exposed population for different α values
f2 = figure(2); hold on;
for k = 1:4
    plot(tt, sum(E_all{k}, 1), colors{k}, 'LineWidth', lineWidth);
end
xlabel('Time', 'FontSize', fontSize);
ylabel('Number of Exposed Individuals', 'FontSize', fontSize);
title('Exposed Population Over Time for Different \alpha Values', 'FontSize', fontSize);
legend(alpha_labels, 'Location', 'best');
grid on;
set(gca, 'FontSize', fontSize);

%% Infected population for different α values
f3 = figure(3); hold on;
for k = 1:4
    plot(tt, sum(I_all{k}, 1), colors{k}, 'LineWidth', lineWidth);
end
xlabel('Time', 'FontSize', fontSize);
ylabel('Number of Infected Individuals', 'FontSize', fontSize);
title('Infected Population Over Time for Different \alpha Values', 'FontSize', fontSize);
legend(alpha_labels, 'Location', 'best');
grid on;
set(gca, 'FontSize', fontSize);

%% Recovered population for different α values
f4 = figure(4); hold on;
for k = 1:4
    plot(tt, sum(R_all{k}, 1), colors{k}, 'LineWidth', lineWidth);
end
xlabel('Time', 'FontSize', fontSize);
ylabel('Number of Recovered Individuals', 'FontSize', fontSize);
title('Recovered Population Over Time for Different \alpha Values', 'FontSize', fontSize);
legend(alpha_labels, 'Location', 'best');
grid on;
set(gca, 'FontSize', fontSize);

%% Save figures
saveas(f1, 'S_varying_alpha.jpg');
saveas(f2, 'E_varying_alpha.jpg');
saveas(f3, 'I_varying_alpha.jpg');
saveas(f4, 'R_varying_alpha.jpg');
saveas(f5, 'SEIR_varying_alpha.jpg');