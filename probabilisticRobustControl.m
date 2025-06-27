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
% Date last updated: 26 June, 2025.
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

%% Nominal SISO System P1
A1 = [-1];
B1 = [1];
C1 = [1];
D1 = 0;
sys1 = ss(A1, B1, C1, D1);
[num1, den1] = tfdata(tf(sys1), 'v');

%% Parameters
d = 2;
N_samples = 1000;
sigma = 0.5;
theta0 = [0.1; -0.05];
b_PC = 0.8;
gamma_vals = linspace(1.01, 3.0, 100);

theta_samples = mvnrnd(theta0', sigma^2 * eye(d), N_samples);

gap_vals = zeros(N_samples,1);
theta_norms = zeros(N_samples,1);
Tzw_vals = zeros(N_samples,1);

% Nominal "graph" operator
G1 = [num1, den1]';
G1 = G1 / norm(G1, 'fro');
P_G1 = G1 * ((G1' * G1) \ G1');   % Projection matrix
P_G1_perp = eye(size(P_G1)) - P_G1;

for i = 1:N_samples
    theta = theta_samples(i,:)';

    %% Perturbed System P2(theta)
    A2 = -1 + theta(1);
    C2 = 1 + theta(2);
    sys2 = ss(A2, B1, C2, D1);
    [num2, den2] = tfdata(tf(sys2), 'v');

    G2 = [num2, den2]';
    if norm(G2, 'fro') == 0
        gap_vals(i) = NaN;
        continue;
    end
    G2 = G2 / norm(G2, 'fro');

    %% Gap approximation
    gap_vals(i) = norm(P_G1_perp * G2, 2);
    theta_norms(i) = norm(theta - theta0);

    %% Closed-loop H-infinity norm
    Cfb = place(A1, B1, [-2]);
    Acl = A2 - B1 * Cfb;
    sys_cl = ss(Acl, 1, 1, 0);
    Tzw_vals(i) = norm(sys_cl, Inf);
end

%% Remove NaN entries
gap_vals = gap_vals(~isnan(gap_vals));
theta_norms = theta_norms(1:length(gap_vals));
Tzw_vals = Tzw_vals(1:length(gap_vals));

%% Plot Gap vs Theta Norm
figure;
scatter(theta_norms, gap_vals, 20, 'filled');
xlabel('$||\theta - \theta_0||$', 'Interpreter', 'latex');
ylabel('$\delta_{g}(\bar{\Sigma}, \tilde{\Sigma}(\theta))$', 'Interpreter', 'latex');
grid on;
grid on;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 5);
set(a, 'linewidth', 5);
set(a, 'FontSize', 50);
set(gca,'fontweight','bold');

%% Probability Curve for H-infinity satisfaction
C_est = mean(gap_vals ./ theta_norms);
probs = zeros(size(gamma_vals));
for i = 1:length(gamma_vals)
    delta_bar = (gamma_vals(i) - b_PC) / (1 + gamma_vals(i));
    probs(i) = 1 - exp(-(delta_bar^2)/(2*C_est^2*sigma^2));
end

%% Summary
fprintf('--- Summary ---\n');
fprintf('Estimated C: %.4f\n', C_est);
fprintf('E[Gap]: %.4f\n', mean(gap_vals));
fprintf('C * E[norm]: %.4f\n', C_est * mean(theta_norms));
fprintf('E[||T_{zw}(P2, C)||_inf]: %.4f\n', mean(Tzw_vals));


figure;
plot(gamma_vals, probs, 'LineWidth', 2);
xlabel('$\gamma$', 'Interpreter', 'latex');
ylabel('$P(||T_{zw}(\tilde{\Sigma}(\theta), C)||_{\mathcal{H}_{\infty}} \leq \gamma)$', 'Interpreter', 'latex');
grid on;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 5);
set(a, 'linewidth', 5);
set(a, 'FontSize', 50);
set(gca,'fontweight','bold');