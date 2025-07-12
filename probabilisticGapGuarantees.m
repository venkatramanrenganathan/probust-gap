%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code simulates Probabilistic Robust Control Using Gap Metric
%
% Copyrights @ 2025
% 
% Authors: Venkatraman Renganathan 
%          Cranfield University, United Kingdom.
%
% Email: v.renganathan@cranfield.ac.uk
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

%% Nominal SISO System \bar{\Sigma}
A1 = -1;
B1 = 1;
C1 = 1;
D1 = 0;
nominalSystem = tf(ss(A1, B1, C1, D1));

% Set the desired closed loop pole location
desiredPoleLocation = -2;

% Place the pole at desiredPoleLocation & get nominal controller gain
nominalController = place(A1, B1, desiredPoleLocation);

% Get the state space and the transfer function
nominalControllerStateSpace = ss(nominalController);

% Compute the performance measure for nominal plant and controller
bPC_Nominal = compute_bpc(nominalSystem, nominalControllerStateSpace); % 0.8;

%% Parameters for Simulation

% Set the number of samples of \theta for Monte-Carlo simulation
numSamples = 10000;

% Mean of \theta
meanTheta = [0.1; -0.05];

% Infer the dimension of \theta parameter
thetaDimension = size(meanTheta, 1);

% Set the standard deviation for \theta parameter
standardDeviationTheta = 0.5;

% Form the covariance matrix of \theta parameter
covarianceTheta = standardDeviationTheta^(2)*eye(thetaDimension);

% Sample \theta from Gaussian distribution
thetaSamples = mvnrnd(meanTheta', covarianceTheta, numSamples);

% Place holders for storing computed values
gapValues = zeros(numSamples,1);
thetaNorms = zeros(numSamples,1);
hInftyNormValues = zeros(numSamples,1);

% Do Monte-Carlo Simulation with all \theta samples
for i = 1:numSamples
    
    % Get the ith \theta sample
    theta = thetaSamples(i,:)';

    % Form ith Perturbed System \tilde{\Sigma}(theta) using \theta
    A2 = A1 + theta(1);
    C2 = C1 + theta(2);
    perturbedSystem = tf(ss(A2, B1, C2, D1));
    
    % Compute gap with Matlab inbuilt command
    [gapValues(i), ~] = gapmetric(nominalSystem, perturbedSystem);
    
    % Compute the norm of theta variations from mean
    thetaNorms(i) = norm(theta - meanTheta);
    
    % Form the closed loop state space
    closedLoopSystem = feedback(series(nominalController, perturbedSystem), 1);   
    
    % Compute H-infinity norm of Closed-loop TF
    hInftyNormValues(i) = hinfnorm(closedLoopSystem);
end

% Remove NaN & Infinity entries
gapValues = gapValues(~isnan(gapValues));
thetaNorms = thetaNorms(1:length(gapValues));
hInftyNormValues = hInftyNormValues(~isinf(hInftyNormValues));

% Computing Bound on Probability of H-infinity satisfaction
% Set the gamma range
gammaValues = linspace(1.01, 10.0, 100);

% Estimate Lipschitz constant of gap as Gap/theta_norm
gapLispchitz = mean(gapValues ./ thetaNorms);

% Placeholder for Upper bound of Prob(\norm{T_zw}_{\infty} \leq \gamma)
upperBoundProbability = zeros(size(gammaValues));

% For every gamma, find the upperBoundProbability
for i = 1:length(gammaValues)
    % Get gammaBar
    gammaBar = (gammaValues(i) - bPC_Nominal) / (1 + gammaValues(i));
    
    % Find Upper bound for Prob(\norm{T_zw}_{\infty} \leq \gamma)
    upperBoundProbability(i) = 1 - exp(-(gammaBar^2)/(2*gapLispchitz^2*standardDeviationTheta^2));
end

%% Print the Summary of Results
fprintf('--- Printing the Summary of Results ---\n');

% Report the Lipschitz constant of gap
fprintf('Estimated L_gap: %.4f\n', gapLispchitz);

% Compute & report the expected gap
expectedGap = mean(gapValues);
fprintf('E[Gap]: %.4f\n', expectedGap);

% Compute & report the upper bound for expected gap
expectedGapUpperBound = gapLispchitz * mean(thetaNorms);
fprintf('L_gap * E[theta_norms]: %.4f\n', expectedGapUpperBound);

% Check if expectedGap exceeds the upperBound
if(expectedGap <= expectedGapUpperBound)
    fprintf('E[Gap] <= L_gap*E[theta_norms] \n');
else
    fprintf('E[Gap] > L_gap*E[theta_norms] \n');
end

% Compute & report the expected Hinf norm of T_zw transfer function
expectedTzwHinfNorm = mean(hInftyNormValues);
fprintf('E[||T_{zw}(P2, C)||_inf]: %.4f\n', expectedTzwHinfNorm);

% Compute & report the upper bound on expected Hinf norm of T_zw transfer function
CinvGap = (1/(1-expectedGapUpperBound))*(1+expectedGapUpperBound + 8 * standardDeviationTheta^(2) * gapLispchitz^(2) * exp(-(1-expectedGapUpperBound)^(2)/(8 * standardDeviationTheta^(2) * gapLispchitz^(2) )));
expectedTzwHinfNormUpperBound = (bPC_Nominal + 1) * CinvGap;
fprintf('Upper Bound on E[||T_{zw}(P2, C)||_inf]: %.4f\n', expectedTzwHinfNormUpperBound);


%% Plotting Code
% Plot Gap vs Theta Norm
figure;
scatter(thetaNorms, gapValues, 150, 'filled');
xlabel('$||\theta - \mu_{\theta}||$', 'Interpreter', 'latex');
ylabel('$\mathrm{Gap}(\theta)$', 'Interpreter', 'latex');
grid on;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 5);
set(a, 'linewidth', 5);
set(a, 'FontSize', 70);
set(gca,'fontweight','bold');

% Plot gamma vs upperBoundProbability
figure;
plot(gammaValues, upperBoundProbability, 'LineWidth', 2);
xlabel('$\gamma$', 'Interpreter', 'latex');
ylabel('$P(||T_{zw}(\tilde{\Sigma}(\theta), \bar{C})||_{\mathcal{H}_{\infty}} \leq \gamma)$', 'Interpreter', 'latex');
grid on;
xticks(1:1:10);
xlim([1, 10]);
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 5);
set(a, 'linewidth', 5);
set(a, 'FontSize', 50);
set(gca,'fontweight','bold');
