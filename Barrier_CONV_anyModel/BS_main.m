addpath(genpath('.'))
%% Parameters
S_0 = 218.17;
% Market parameters
param.rf = 0.02;             % riskfree interest rate 
% Model parameters
param.sigma = 0.2516;        % BM vol
% Contract parameters
param.T = 1;                  % maturity
K = S_0;                 % strike
M = round(12*param.T);        % monthly monitoring
param.q= 0; %dividend
param.dt= param.T/M;

param.distr=1; % BS
Barrier=0.98*S_0; N=2^12;
[S,v] = CONV( S_0, K, M, N, Barrier, param);
price=interp1(S,v,S_0,'spline')


close all
