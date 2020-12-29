%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price a Barrier Down&Out call option under Merton model 
% using MC with antithetic variance reduction 
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
D = 90;                 % barrier
M = 52*T;               % weekly monitoring
disc_payoff_fun = @(S) exp(-r*T).*max(S(:,end)-K,0).*(min(S,[],2)>D);     % disc payoff function. S(i,:) = i-th simulation of an underlying PATH 
% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'mu',mu,'delta',delta,'lambdaK',lambdaK);

%% Simulate Underlying Asset
% Simulate the underlying at time T
[S,SAV] = Merton_simulate_asset_AV(par,Nsim,M);
%% Compute the discounted payoff
DiscPayoff = disc_payoff_fun(S);
DiscPayoffAV = disc_payoff_fun(SAV);

%% Compute call price and asymptotic CI 
disp("Merton model - Barrier DO option price via plain MC:")
[barrier_DO_price, ~, barrier_DO_CI_price] = normfit( (DiscPayoff+DiscPayoffAV)/2 )