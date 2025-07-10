
% coprime_gap_vs_stability.m
% Plot coprime perturbation norm vs closed-loop stability indicator

clear; clc;

% Nominal plant P(s) = 1 / (s + 1)
P = tf(1, [1 1]);
% [N, M] = coprime(P);
[~,M,N] = rncf(P);

% Controller C(s)
s = tf('s');
C = 5 * (s + 2)/(s + 10);

% Compute b_PC using full formula
G = [P; 1];                        % [P; I]
L = eye(1) - C*P;
invL = inv(L);
W = [-C 1];
T = G * invL * W;
b_PC = 1 / hinfnorm(T);

% Simulation parameters
N_samp = 1000;
sigma = 0.2;
L_lip = 2;

coprime_radius = zeros(1, N_samp);
is_stable = false(1, N_samp);

for i = 1:N_samp
    eps = sigma * randn(2,1);
    dN = tf(eps(1), [10 1]);   % stable perturbations
    dM = tf(eps(2), [10 1]);
    Np = N + dN;
    Mp = M + dM;

    try
        Ptheta = Np / Mp;
        CL = feedback(Ptheta * C, 1);
        poles = pole(CL);
        is_stable(i) = all(real(poles) < 0);
        d_coprime = [dN; dM];
        coprime_radius(i) = norm(d_coprime, inf);
    catch
        coprime_radius(i) = NaN;
        is_stable(i) = false;
    end
end

% Remove NaNs
valid = ~isnan(coprime_radius);
coprime_radius = coprime_radius(valid);
is_stable = is_stable(valid);

% Plot
figure;
scatter(coprime_radius, is_stable, 30, is_stable, 'filled');
xlabel('Coprime perturbation norm (H_\infty)');
ylabel('Stabilized (1 = yes, 0 = no)');
title('Stability vs Coprime Perturbation Radius');
grid on;
