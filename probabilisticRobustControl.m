%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code simulates Probabilistic Robust Control Using Gap Metric
%
% Copyrights Authors: Venkatraman Renganathan 
%                     Cranfield University, United Kingdom.
%
% Emails: v.renganathan@cranfield.ac.uk
%
% Date last updated: 27 June, 2025.
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
sys1 = ss(A1, B1, C1, D1);
[num1, den1] = tfdata(tf(sys1), 'v');

% Prepare Nominal "Graph" operator foor normalized RCF
G1 = [num1, den1]';
G1 = G1 / norm(G1, 'fro');
% Pi_G1 = G1 * ((G1' * G1) \ G1');   % Projection matrix
Pi_G1 = G1 * G1';   % Projection matrix
Pi_G1_perp = eye(size(Pi_G1)) - Pi_G1;

%% Parameters for Simulation
% Set the performance measure for nominal plant and controller
b_PC = 0.8;
% Set the standard deviation for \theta parameter
sigma = 0.5;
% Set the number of samples of \theta for Monte-Carlo simulation
numSamples = 1000;
% Mean of \theta
barTheta = [0.1; -0.05];
% Sample \theta from Gaussian distribution
thetaSamples = mvnrnd(barTheta', sigma^(2)*eye(2), numSamples);

% Place holders 
gapValues = zeros(numSamples,1);
thetaNorms = zeros(numSamples,1);
TzwValues = zeros(numSamples,1);

% Do Monte-Carlo Simulation with all \theta samples
for i = 1:numSamples
    
    % Get the ith \theta sample
    theta = thetaSamples(i,:)';

    % Form ith Perturbed System \tilde{\Sigma}(theta) using \theta
    A2 = -1 + theta(1);
    C2 = 1 + theta(2);
    sys2 = ss(A2, B1, C2, D1);
    [num2, den2] = tfdata(tf(sys2), 'v');
    G2 = [num2, den2]';
    % If Frobenius norm of G2 is 0, then Gap is undefined
    if norm(G2, 'fro') == 0
        gapValues(i) = NaN;
        continue;
    end
    % Normalize it
    G2 = G2 / norm(G2, 'fro');

    % Compute gap as \norm{Pi_G1_perp * G2}
    % gapValues(i) = norm(Pi_G1_perp * G2, 2);
    % Compare with Matlab inbuilt gapmetric command
    [gapValues(i),~] = gapmetric(tf(sys1),tf(sys2));
    % Compute the norm of theta variations from mean
    thetaNorms(i) = norm(theta - barTheta);

    % Compute H-infinity norm of Closed-loop TF
    % Set the desired closed loop pole location
    desiredPoleLocation = -2;
    % Place the pole at desiredPoleLocation & get nominal controller gain
    barC = place(A1, B1, desiredPoleLocation);
    % Form the closed loop matrix A - BK
    Acl = A2 - B1 * barC;
    % Form the closed loop state space
    sys_cl = ss(Acl, 1, 1, 0);
    % Compute the Hinfty norm
    TzwValues(i) = norm(sys_cl, Inf);
end

% Remove NaN entries
gapValues = gapValues(~isnan(gapValues));
thetaNorms = thetaNorms(1:length(gapValues));
TzwValues = TzwValues(1:length(gapValues));

% Computing Bound on Probability of H-infinity satisfaction
% Set the gamma range
gammaValues = linspace(1.01, 3.0, 100);
% Estimate Lipschitz constant of gap as Gap/theta_norm
L_gap = mean(gapValues ./ thetaNorms);
% Placeholder for Upper bound of Prob(\norm{T_zw}_{\infty} \leq \gamma)
upperBoundProbability = zeros(size(gammaValues));
% For every gamma, find the upperBoundProbability
for i = 1:length(gammaValues)
    gammaBar = (gammaValues(i) - b_PC) / (1 + gammaValues(i));
    % Find Upper bound for Prob(\norm{T_zw}_{\infty} \leq \gamma)
    upperBoundProbability(i) = 1 - exp(-(gammaBar^2)/(2*L_gap^2*sigma^2));
end

%% Summary
fprintf('--- Summary ---\n');
fprintf('Estimated L_gap: %.4f\n', L_gap);
% Compute & report the expected gap
expectedGap = mean(gapValues);
fprintf('E[Gap]: %.4f\n', expectedGap);
% Compute & report the upper bound for expected gap
expectedGapUpperBound = L_gap * mean(thetaNorms);
fprintf('L_gap * E[theta_norms]: %.4f\n', expectedGapUpperBound);
if(expectedGap <= expectedGapUpperBound)
    fprintf('E[Gap] <= L_gap*E[theta_norms]. So, Lemma 4 is satisfied \n');
else
    fprintf('E[Gap] > L_gap*E[theta_norms]. So, Lemma 4 not satisfied \n');
end
% Compute & report the expected Hinf norm of T_zw transfer function
expectedTzwHinfNorm = mean(TzwValues);
fprintf('E[||T_{zw}(P2, C)||_inf]: %.4f\n', expectedTzwHinfNorm);

%% Plotting Code
% Plot Gap vs Theta Norm
figure;
scatter(thetaNorms, gapValues, 80, 'filled');
xlabel('$||\theta - \theta_0||$', 'Interpreter', 'latex');
ylabel('$\delta_{g}(\bar{\Sigma}, \tilde{\Sigma}(\theta))$', 'Interpreter', 'latex');
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
ylabel('$P(||T_{zw}(\tilde{\Sigma}(\theta), C)||_{\mathcal{H}_{\infty}} \leq \gamma)$', 'Interpreter', 'latex');
grid on;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 5);
set(a, 'linewidth', 5);
set(a, 'FontSize', 70);
set(gca,'fontweight','bold');
