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
% Date last updated: 12 July, 2025.
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
C = 0.5*(s+2)/(s+3); %0.5; % Controller stabilizing nominal system

% Gap robustness margin
b_PC = compute_bpc(P, C);

% Simulation Parameters
N = 10000;            % Number of scenarios
beta = 0.01;          % confidence level
epsilon = 0.05;       % violation probability
muTheta = 0;          % Mean of theta
sigmaTheta = 0.25;    % Standard deviation of theta

% Generate N scenarios of theta
thetaValues = normrnd(muTheta, sigmaTheta, [1, N]);

% Place holder to store gap values
gapValues = zeros(N,1);

% Run Monte-carlo simulation with N scenarios
for i = 1:N
    % Generate the pole perturbation as theta
    deltaPole = thetaValues(i);
    
    % Generate the random perturbed system
    P_theta = 1/(s+1+deltaPole);
    
    % Compare with Matlab inbuilt gapmetric command
    [gapValues(i), ~] = gapmetric(P, P_theta);
end

% Scenario-based gap threshold
alphaHatN = max(gapValues);

% Theoretical sample size check (Calafiore & Dabbene)
requiredSamples = ceil(log(1/beta)/log(1/(1-epsilon)));
fprintf('Scenarios taken: %d, Required: %d\n', N, requiredSamples);

% Probabilistic stabilization check
if alphaHatN < b_PC
    fprintf('Probabilistic robustness guarantee met: alpha_hat=%.4f < b_PC=%.4f\n', alphaHatN, b_PC);
else
    fprintf('Probabilistic robustness guarantee NOT met: alpha_hat=%.4f >= b_PC=%.4f\n', alphaHatN, b_PC);
end

%% Plot scenario gap values
figure;
histogram(gapValues, 30,'Normalization','probability');
xline(alphaHatN,'r','LineWidth',5);
xline(b_PC,'k','LineWidth',5);
xlabel('$\mathrm{Gap}(\theta)$');
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