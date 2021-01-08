%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options under NIG model using MC
% with antithetic variance reduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.01;            % riskfree interest rate 
S0 = 100;            % spot price
% Model parameters
sigma = 0.6;
theta = 0.05;
kVG   = 0.2;
% Contract parameters
T = 1;   % maturity
K = S0;  % strike

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'theta',theta,'kNIG',kNIG);
% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

%% Simulate Underlying Asset
% Simulate the underlying at time T
[~,~,ST,STAV] = NIG_simulate_assetAV(par,Nsim,1);
%% Compute the discounted payoff
DiscPayoff = exp(-r * T) * max(ST - K, 0);
DiscPayoffAV = exp(-r * T) * max(STAV - K, 0);

%% Compute call price and asymptotic CI 
disp("NIG Model - EU Vanilla options price via MC with antithetic variance reduction:")
[call_price, ~, call_CI_price] = normfit( (DiscPayoff+DiscPayoffAV)/2 )

%% Compute put price and asymptotic CI via Call-Put Parity
put_parity = @(call_p) call_p - S0 + K*exp(-r*T);

put_price = put_parity(call_price)
put_CI_price = put_parity(call_CI_price)


CI_len = call_CI_price(2)-call_CI_price(1)