%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options under Heston model
% using Andersen QE techinque
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
theta = 0.1;  % mean 
V0 = 0.1; 
% Contract parameters
T = 1;      % maturity
K = S0;     % strike
M = 52;     % number of time discretization steps

par = struct('S0',S0,'r',r,'TTM',T,'epsilon',epsilon,'kappa',k,'rho',rho,'theta',theta,'V0',V0);
% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

%% Simulate Underlying Asset
% Simulate the underlying at time T
rng(0)
[~,ST] = Heston_simulate_asset_QE(par,Nsim,M);
%% Compute the discounted payoff
DiscPayoff = exp(-r * T) * max(ST - K, 0);

%% Compute call price and asymptotic CI 
disp("Heston Model - EU Vanilla options price via QE simulation:")
[call_price, ~, call_CI_price] = normfit(DiscPayoff)

%% Compute put price and asymptotic CI via Call-Put Parity
put_parity = @(call_p) call_p - S0 + K*exp(-r*T);

put_price = put_parity(call_price)
put_CI_price = put_parity(call_CI_price)

CI_len = call_CI_price(2)-call_CI_price(1)