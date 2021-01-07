%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options under Bates model
% using euler scheme approx of SDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.001;               % riskfree interest rate 
S0 = 100;            % spot price
% Model parameters
epsilon = 0.2;% vol-of-vol
k = 0.01;     % mean reversion speed
rho = -0.2;   % correlation
theta = 0.1;  % long-term mean
V0 = 0.1;     % initial var-of-vol
kZ = 0;       % log(1+kZ) = mean of jump size
lambda = 0;   % intensity of jump time
sigma = 0;    % vol of jump size

% Contract parameters
T = 1;      % maturity
K = S0;     % strike
M = 52;     % number of time discretization steps

par = struct('S0',S0,'r',r,'TTM',T,'epsilon',epsilon, ...
'kappa',k,'rho',rho,'theta',theta,'V0',V0,'lambda',lambda,'kappaZ',kZ,'sigma',sigma);
% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

%% Simulate Underlying Asset
% Simulate the underlying at time T
rng(0)
[~,ST] = Bates_simulate_asset(par,Nsim,M);
%% Compute the discounted payoff
DiscPayoff = exp(-r * T) * max(ST - K, 0);

%% Compute call price and asymptotic CI 
disp("Bates Model - EU Vanilla options price via euler simulation:")
[call_price, ~, call_CI_price] = normfit(DiscPayoff)

%% Compute put price and asymptotic CI via Call-Put Parity
put_parity = @(call_p) call_p - S0 + K*exp(-r*T);

put_price = put_parity(call_price)
put_CI_price = put_parity(call_CI_price)

CI_len = call_CI_price(2)-call_CI_price(1)