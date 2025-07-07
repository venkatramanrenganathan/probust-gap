%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code simulates Scenario-Based Gap Robustness Numerical Illustration
%
% Copyrights Authors: Venkatraman Renganathan 
%                     Cranfield University, United Kingdom.
%
% Emails: v.renganathan@cranfield.ac.uk
%
% Date last updated: 7 July, 2025.
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

% Nominal system P
s = tf('s');
P = 1/(s+1);
C = 1; % Controller stabilizing nominal system
sen = 1/(1+P*C);
CompSenNorm = norm(1-sen,inf);

% Gap robustness margin (approximate)
b_PC = norm(1-sen,inf); % assumed known from design 0.6

% Parameters
d = 2; % dimension
theta0 = zeros(d,1);
L = 0.5; % gap Lipschitz constant
sigma = 0.1;
N = 1000; % Number of scenarios
epsilon = 0.05; % violation probability
beta = 0.01; % confidence level

% Generate scenarios
thetas = mvnrnd(theta0, sigma^2*eye(d), N);

gap_vals = zeros(N,1);
for i=1:N
    theta = thetas(i,:)';
    % Random perturbed system
    delta_pole = 0.5*sum(theta);
    P_theta = 1/(s+1+delta_pole);
    
    % Simple gap estimate using H-inf norm of difference
    gap_vals(i) = norm(P-P_theta,inf)/norm(P,inf);
end

% Scenario-based gap threshold
alpha_hat = max(gap_vals);

% Theoretical sample size check (Calafiore & Dabbene)
N_required = ceil(log(1/beta)/log(1/(1-epsilon)));
fprintf('Scenarios taken: %d, Required: %d\n', N, N_required);

% Probabilistic stabilization check
if alpha_hat < b_PC
    fprintf('Probabilistic robustness guarantee met: alpha_hat=%.4f < b_PC=%.4f\n', alpha_hat, b_PC);
else
    fprintf('Probabilistic robustness guarantee NOT met: alpha_hat=%.4f >= b_PC=%.4f\n', alpha_hat, b_PC);
end

% Plot scenario gap values
figure;
histogram(gap_vals,30,'Normalization','probability');
xline(alpha_hat,'r','LineWidth',5);
xline(b_PC,'k','LineWidth',5);
xlabel('Gap');
ylabel('Probability Density');
xlim([0, 1]);
legend('Scenario $\mathbf{f}_{gap}$','$\hat{\alpha}_{N}$','$b_{\bar{\Sigma}, \bar{C}}$');
grid on;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 5);
set(a, 'linewidth', 5);
set(a, 'FontSize', 50);
set(gca,'fontweight','bold');