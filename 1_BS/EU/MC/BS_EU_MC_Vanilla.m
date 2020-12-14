%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put and Options in the B&S model using MC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.02;               % riskfree interest rate 
S0 = 218.75;            % spot price
% Model parameters
sigma = 0.2516;         % standard deviation 
% Contract parameters
T = 1;                  % maturity
K = S0;                 % strike
% Discretization parameter
Nsim = 1e8;             % number of MC simulations 

%% Simulate Underlying Asset
% Simulate the underlying at time T
ST = S0 * exp((r - sigma^2/2) * T + sigma * sqrt(T) * randn(Nsim, 1));

%% Compute the discounted payoff
DiscPayoff = exp(-r * T) * max(ST - K, 0);

%% Compute call price and asymptotic CI 
disp("EU Vanilla options price via plain MC:")
[call_price, ~, call_CI_price] = normfit(DiscPayoff)

%% Compute put price and asymptotic CI via Call-Put Parity
put_parity = @(call_p) call_p - S0 + K*exp(-r*T);

put_price = put_parity(call_price)
put_CI_price = put_parity(call_CI_price)
