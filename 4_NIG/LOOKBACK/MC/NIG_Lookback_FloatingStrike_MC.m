%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price Lookback floating strike call option under VG model using MC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
% Market parameters
r = 0.01;            % riskfree interest rate 
S0 = 100;            % spot price
% Model parameters
sigma = 0.6;
theta = 0.05;
kNIG   = 0.2;
% Contract parameters
T = 1;                  % maturity
M = round(12*T);        % monthly monitoring

disc_payoff_fun = @(S) exp(-r*T).*max(S(:,end)-min(S,[],2),0);     % disc payoff function. S(i,:) = i-th simulation of an underlying PATH 
par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'theta',theta,'kNIG',kNIG);
% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

%% Simulate Underlying Asset
% Simulate the underlying at time T
S = NIG_simulate_asset(par,Nsim,M);
%% Compute the discounted payoff
DiscPayoff = disc_payoff_fun(S);

%% Compute call price and asymptotic CI 
disp("VG Model - Lookback floating strike call option price via plain MC:")
[lb_floatStrike_call_price, ~, lb_floatStrike_call_price_CI] = normfit(DiscPayoff)

CI_len = lb_floatStrike_call_price_CI(2)-lb_floatStrike_call_price_CI(1)