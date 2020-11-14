clear all
close all
clc

%% INPUT PARAMETERS 
% Market parameteres
r = 0.2;
S0 = 100;

% Contract
T = 1;
Stike = 100;
payoff = @(ST) max(ST-Strike,0);
% Model parameters
sigma = 0.4;

% Domain Boundaries 
Smax = S0*exp( (r - sigma^2/2)*T + 6*sigma*sqrt(T));

% Discretization parameters 
M = 50; dt = T/M; t_grid = 0:dt:M;
N = 30; ds = Smax/N; s_grid = 0:ds:N;

% Resolving matrix              
low_price_coef = -r*s_grid/ds/2 + sigma^2*s_grid/2/(ds^2);
mid_price_coef = -1/dt -s_grid.^2*sigma^2/(ds^2)-r;
high_price_coef = r*s_grid*2/ds + sigma^2/2/(ds^2) * s_grid.^2; 
A = spdiags([low_price_coef mid_price_coef high_price_coef], -1:1, N+1, N+1);
b = -1/dt * payoff(s_grid)

for j = (M-1):-1:1
    %boundary adjustment
    b(1) = b(1) - 0; % bottom bdn adjustment
    b(end) = b(end) - F(end)*(r*Smax/2/ds + sigma^2*Smax^2/2/(ds^2)); % upper bnd adjustment
    %compute prices at time tj
    F = A\b;
end

callPrice = interp1(s_grid, F, S0)
