%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options under Kou model using MC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.001;              % riskfree interest rate 
S0 = 1;                 % spot price
% Model parameters
sigma = 0.4;            % BM vol
p = 0.5;                % jump vol
lambdap = 5;            % pos jump intensity
lambdam = 6;            % neg jump intensity
lambdaK = 3;            % n jumps intensity
% Contract parameters
T = 1;                  % maturity
K = S0;                 % strike
% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'p',p,'lambdap',lambdap,'lambdam',lambdam,'lambdaK',lambdaK);

%% Simulate Underlying Asset
% Simulate the underlying at time T
[~,ST] = Kou_simulate_asset(par,Nsim,1);
%% Compute the discounted payoff
DiscPayoff = exp(-r * T) * max(ST - K, 0);

%% Compute call price and asymptotic CI 
disp("Kou Model - EU Vanilla options price via plain MC:")
[call_price, ~, call_CI_price] = normfit(DiscPayoff)

%% Compute put price and asymptotic CI via Call-Put Parity
put_parity = @(call_p) call_p - S0 + K*exp(-r*T);

put_price = put_parity(call_price)
put_CI_price = put_parity(call_CI_price)
