%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is useful while simulating Probabilistic Robust Control 
% using random gap metric formulation
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

function b_PC = compute_bpc(P, C)

% Function compute_bpc computes the robust stability margin b_{P,C} using
% b_{P,C} = 1 / || [P; I](I - CP)^{-1}[-C I] ||_∞
%
% Inputs:
%   P - Plant as transfer function or ss system
%   C - Controller as transfer function or ss system
%
% Output:
%   b_PC - Robust stability margin from gang of four matrices

% Check dimensions
if size(P,2) ~= size(C,1)
    error('Inner dimensions of P and C do not match');
end

% Step 1: Compute CP = C * P
CP = series(C, P);

% Step 2: Compute (I - CP)^{-1}
I_sys = eye(size(CP));

try
    inv_I_minus_CP = feedback(I_sys, CP);  % inv(I - CP)
catch
    error('(I - CP) is not invertible (possibly unstable closed-loop)');
end

% Step 3: Build [P; I]
I_tf = eye(size(P));           % identity of same size as P
PI_block = [P; I_tf];          % vertical concatenation

% Step 4: Build [-C I]
minusC_I_block = [-C, I_tf];   % horizontal concatenation

% Step 5: Form transfer H_PC = [P; I]*(I - CP)^{-1} * [-C I]
H_PC = PI_block * inv_I_minus_CP * minusC_I_block;

% Step 6: Use hinfnorm (accurate H-infinity norm)
[gamma, ~] = hinfnorm(H_PC);

% Step 7: Compute b_{P,C}
if gamma > 0
    b_PC = 1 / gamma;
else
    b_PC = 0;
end

end