% set properties for plotting
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
addpath(genpath('src'));

% Parameters
b_PC = 0.8;           % nominal performance
C_delta = 0.5;        % Lipschitz constant
gamma = linspace(b_PC + 0.01, 10.0, 300);  % performance threshold

% Compute delta_bar and probability
delta_bar = (gamma - b_PC) ./ (1 + gamma);
probability = 1 - exp(-(delta_bar.^2) ./ (2 * C_delta^2));

% Plot
figure;
plot(gamma, probability, 'LineWidth', 5);
xlabel('$\gamma$ ', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$P(\|\mathbf{T}(\tilde{\Sigma}, \bar{C})\|_{\mathcal{H}_{\infty}} \leq \gamma)$', ...
       'Interpreter', 'latex', 'FontSize', 20);
grid on;
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 5);
set(a, 'linewidth', 5);
set(a, 'FontSize', 50);
set(gca,'fontweight','bold');