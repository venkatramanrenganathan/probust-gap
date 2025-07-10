%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code simulates to obtain probability of robust stability using
% gap metric between LTI dynamical systems under random parameter θ
%
% Copyrights @ 2025
% 
% Authors: Venkatraman Renganathan 
%          Cranfield University, United Kingdom.
%
% Email: v.renganathan@cranfield.ac.uk
%
% Date last updated: 10 July, 2025.
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
nominalSystem = tf(ss(A, B, C, D));

% Set the parameters for simulation
p = 2;                          % Dimension of θ 
numSamples = 10000;             % Number of θ samples
meanTheta = [0.50 -0.25];       % Mean of θ
standardDeviationTheta = 0.25;  % Standard Deviation of θ

% Set covariance of theta
covarianceTheta = standardDeviationTheta^(2)*eye(p); 

% Place holder for storing gap and theta values
gapValues = zeros(1, numSamples);
thetaNorms = zeros(1, numSamples);

% Set the stability flag to false
stableFlag = 0;

% Compute gap for all numSamples of theta
for i = 1:numSamples
    
    % % Loop until you set the stable flag = 1 (find a stable system)
    % while(stableFlag == 0)
    % 
    %     % Sample θ from Gaussian distribution
    %     theta = mvnrnd(meanTheta, covarianceTheta);
    % 
    %     % Perturb the A matrix - Small perturbation around stable pole
    %     A_perturb = A + theta(1);   
    % 
    %     % If perturbed system is stable, exit the while loop
    %     if (A_perturb < 0) 
    %         % Perturb C matrix
    %         C_perturb = C + theta(2);
    %         % Set stableFlag as 1 (we have found a stable perturbedSystem)
    %         stableFlag = 1;
    %         break;
    %     end
    % end

    % Sample θ from Gaussian distribution
    theta = mvnrnd(meanTheta, covarianceTheta);

    % Perturb the A matrix - Small perturbation around stable pole
    A_perturb = A + theta(1);   

    % Perturb C matrix
    C_perturb = C + theta(2);

    % Form the tranfer function from state space of the perturbed model
    perturbedSystem = tf(ss(A_perturb, B, C_perturb, D));

    % Compute gap with Matlab inbuilt gapmetric command
    [gapValues(i),~] = gapmetric(nominalSystem, perturbedSystem);
    
    % Store the theta values
    thetaNorms(i) = norm(theta - meanTheta);

end


%% Compute & compare the empirical and theoretical stats

% Find the empirical expected gapvalue
expectedGap = mean(gapValues);

% Estimate Lipschitz constant of gap as Gap/theta_norm
gapLipschitzConstant = mean(gapValues ./ thetaNorms);

% sub-Gaussian variance proxy of Gap(θ)
sigmaGap = gapLipschitzConstant * standardDeviationTheta;

% Set the desired closed loop pole location
desiredPoleLocation = -2;

% Place the pole at desiredPoleLocation & get nominal controller gain
nominalController = place(A, B, desiredPoleLocation);

% Get the state space and the transfer function
nominalControllerStateSpace = ss(nominalController);

% Compute the performance measure for nominal plant and controller
b_PC = compute_bpc(nominalSystem, nominalControllerStateSpace);

% Compute the tolerance (Theorem 6 in Paper)
epsilonTolerance = b_PC - expectedGap;

% Compute Empirical probability of robust stability P[Gap < b_PC]
empiricalRobustStabilityProbability = sum(gapValues < b_PC)/numSamples;

% Compute Theoretical Lower Bound on probability of robust stability P[Gap < b_PC]
lowerBoundRobustStabilityProbability = 1 - exp(-epsilonTolerance^(2) / (2 * sigmaGap^(2)));

%% Display summary
fprintf('\n--- Summary of Simulation Results ---\n');
fprintf('Expected[Gap]:                    %.4f\n', expectedGap);
fprintf('Robust Stability Margin b_PC:     %.4f\n', b_PC);
fprintf('Empirical Probability of Robust Stability P[Gap < b_PC]: %.4f\n', empiricalRobustStabilityProbability);
fprintf('Lower Bound on Probability of Robust Stability: %.4f\n', lowerBoundRobustStabilityProbability);


%% Plot the results
figure;
histogram(gapValues, 40, 'Normalization', 'probability', 'FaceAlpha', 0.6);
hold on;
xline(b_PC, 'r--', 'LineWidth', 5);
xline(expectedGap, 'k--', 'LineWidth', 5);
xlabel('Gap($\theta$)');
ylabel('Probability Density');
legend('Empirical $\mathbf{f}_{\mathrm{gap}}$', '$b_{\bar{\Sigma}, \bar{C}}$', 'Expected $\mathrm{Gap}(\theta)$');
grid on;
xlim([0, 1]);
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 5);
set(a, 'linewidth', 5);
set(a, 'FontSize', 50);
set(gca,'fontweight','bold');