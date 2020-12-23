%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options under Merton model using MC
% with antithetic variance reduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.0367;             % riskfree interest rate 
S0 = 100;               % spot price
% Model parameters
sigma = 0.17801;        % BM vol
delta = sqrt(0.4);      % jump vol
mu    = 0.01;           % jump drift
lambdaK = 0.2;          % n jumps intensity
% Contract parameters
T = 1;                  % maturity
K = S0;                 % strike
% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'mu',mu,'delta',delta,'lambdaK',lambdaK);

%% Simulate Underlying Asset
% Simulate the underlying at time T
[~,~,ST,STAV] = Merton_simulate_asset_AV(par,Nsim,1);
%% Compute the discounted payoff
DiscPayoff = exp(-r * T) * max(ST - K, 0);
DiscPayoffAV = exp(-r * T) * max(STAV - K, 0);

%% Compute call price and asymptotic CI 
disp("Merton Model - EU Vanilla options price via plain MC:")
[call_price, ~, call_CI_price] = normfit( (DiscPayoff+DiscPayoffAV)/2 )

%% Compute put price and asymptotic CI via Call-Put Parity
put_parity = @(call_p) call_p - S0 + K*exp(-r*T);

put_price = put_parity(call_price)
put_CI_price = put_parity(call_CI_price)
