%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price Barrier D&O Call option under Kou model using MC
% with Antithetic Variance Reduction
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
L = 0.6;                % barrier
M = round(12*T);        % monthly monitoring
% Discretization parameter
Nsim = 1e6;             % number of MC simulations 

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma,'p',p,'lambdap',lambdap,'lambdam',lambdam,'lambdaK',lambdaK);

%% Simulate Underlying Asset
% Simulate the underlying at time T
[S,SAV,ST,STAV] = Kou_simulate_asset_AV(par,Nsim,M);
%% Compute the discounted payoff
DiscPayoff = exp(-r * T) * max(ST - K, 0).*(min(S,[],2)>L);
DiscPayoffAV = exp(-r * T) * max(STAV - K, 0).*(min(SAV,[],2)>L);

%% Compute call price and asymptotic CI 
disp("Kou Model - Barrier Down and Out Call option price via MC with Antithetic Variance Reduction:")
[DO_call_price_AV, ~, DO_call_CI_price_AV] = normfit( (DiscPayoff+DiscPayoffAV)/2 )
