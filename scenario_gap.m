% Scenario-Based Gap Robustness Numerical Illustration
clear; clc;

% Nominal system P
s = tf('s');
P = 1/(s+1);
C = 1; % Controller stabilizing nominal system

% Gap robustness margin (approximate)
b_PC = 0.6; % assumed known from design

% Parameters
d = 2; % dimension
theta0 = zeros(d,1);
L = 1.0; % gap Lipschitz constant
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
histogram(gap_vals,30,'Normalization','pdf');
xline(alpha_hat,'r','LineWidth',2,'Label','Scenario max gap');
xline(b_PC,'g','LineWidth',2,'Label','Robustness margin b_{PC}');
title('Scenario-Based Gap Robustness');
xlabel('Gap Metric Values');
ylabel('Probability Density');
legend('Scenario Gap Distribution','Scenario Max Gap','Robustness Margin');