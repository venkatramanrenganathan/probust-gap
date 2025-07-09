%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code simulates gap metric between state-space systems under 
% random parameter θ with simpler nominal model
%
% Copyrights Authors: Venkatraman Renganathan 
%                     Cranfield University, United Kingdom.
%
% Emails: v.renganathan@cranfield.ac.uk
%
% Date last updated: 9 July, 2025.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a fresh start
clear; close all; clc;

% set properties for plotting
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
addpath(genpath('src'));

% Nominal system: simple stable first-order system
A = -1;
B = 1;
C = 1;
D = 0;

% Form transfer function from state space of nominal model
P_nom = tf(ss(A, B, C, D));

% Simulation parameters
p = 1;                          % One-dimensional θ
N = 1000;                       % Number of samples
mu_theta = 0.0;                 % Mean of θ
sigma_theta = 0.25;             % Standard Deviation of θ

% Place holder for storing gap and theta values
gap_values = zeros(1, N);
theta_norms = zeros(1, N);

% Compute gap for all N samples
for i = 1:N

    % Sample a θ from Gaussian distribution
    theta = normrnd(mu_theta, sigma_theta);

    % Perturb the A matrix (ensure stability)
    A_pert = A + theta;        % Small perturbation around stable pole
    if A_pert >= 0             % skip unstable systems
        gap_values(i) = NaN;
        continue;
    end

    % Form the tranfer function from state space of the perturbed model
    P_theta = tf(ss(A_pert, B, C, D));

    try
        % Compute gap with Matlab inbuilt gapmetric command
        [gap_values(i),~] = gapmetric(P_nom, P_theta);
    catch
        gap_values(i) = NaN;
    end
    % Store the theta values
    theta_norms(i) = norm(theta);

end

% Remove NaNs
gap_values = gap_values(~isnan(gap_values));
N_valid = length(gap_values);
theta_norms = theta_norms(1:N_valid);


% ------------------------
% Empirical and theoretical stats
% ------------------------
expected_gap = mean(gap_values);
% Estimate Lipschitz constant of gap as Gap/theta_norm
L_gap_est = mean(gap_values ./ theta_norms);

% sub-Gaussian variance prixy of Gap(θ)
sigma_gap = L_gap_est * sigma_theta;
% Controller stabilizing nominal system (already stable)
C = 1; 
% Form Sensitivity & Complementary sensitivity transfer functions
sen = 1/(1+P_nom*C);
b_PC = norm(1-sen,inf);

% Compute the tolerance
epsilon = b_PC - expected_gap;
% Compute Empirical P[Gap < b_PC]
empiricalSatisfactionProbability = mean(gap_values < b_PC);
% Compute Theoretical Lower Bound on P[Gap < b_PC]
theoreticalLowerBound = 1 - exp(-epsilon^2 / (2 * sigma_gap^2));

% Display summary
fprintf('\n--- Summary of Simulation Results ---\n');
fprintf('Valid Samples:                    %d\n', N_valid);
fprintf('Expected[Gap]:                    %.4f\n', expected_gap);
fprintf('Robust Stability Margin b_PC:     %.4f\n', b_PC);
fprintf('Empirical P[Gap < b_PC]:          %.4f\n', empiricalSatisfactionProbability);
fprintf('Theoretical Lower Bound:          %.4f\n', theoreticalLowerBound);

% ------------------------
%% Plot
% ------------------------
figure;
histogram(gap_values, 40, 'Normalization', 'probability', 'FaceAlpha', 0.6);
hold on;
xline(b_PC, 'r--', 'LineWidth', 5);
xlabel('Gap($\theta$)');
ylabel('Probability Density');
legend('Empirical $\mathbf{f}_{\mathrm{gap}}$', '$b_{\bar{\Sigma}, \bar{C}}$');
grid on;
xlim([0, 1]);
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 5);
set(a, 'linewidth', 5);
set(a, 'FontSize', 50);
set(gca,'fontweight','bold');