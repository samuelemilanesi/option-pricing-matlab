%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calibrate Merton model parameters on a set of vanilla options 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

%% Marked to market EU Call options data
data=[0.5  0.3  13.0391
      0.5  0.6  13.2413
      0.5  0.9  13.4410
      0.5  1    13.5070
      0.5  1.5  13.8325
      0.5  2    14.1501
      0.6  0.6  10.7861
      0.7  0.6  8.36675
      0.8  0.6  6.01538
      0.9  0.6  3.79936
      1    0.6  1.87872
      1.1  0.6  0.61598
      1.3  0.6  0.05320
      1.4  0.6  0.01986
      1.5  0.6  0.00832];

S0 = 25.67;    % spot price of the underlying asset 
r=0.05;         % risk free rate 
% vector of market prices
mkt_prices = data(:,3);
% vector of strikes of each traded option
Strikes = S0*data(:,1);
% vector of time to maturity of each traded option
TTMs = data(:,2);

% Market parameters and initial ansaz of model parameters 
par = struct('S0',S0,'r',r,'TTM',0,'sigma',0,'mu',0,'delta',0,'lambdaK',1);
t=tic();

%% RSS loss minimization 
gaps = @(x) fun_Merton_model_mkt_gap(par,x,Strikes, TTMs, mkt_prices);

% Params optimization
options = optimoptions('lsqnonlin','FunctionTolerance',1e-8,'OptimalityTolerance',1e-8,'StepTolerance',1e-8);
opt_params = lsqnonlin(@(x) gaps(x), [0.1 0.5 0.5 10],[0 -10 0 1], [1 1 1 1000], options)    

toc(t)
%% Error and plot
figure
plot(Strikes, mkt_prices,'+'); hold on
[err_vect, model_prices] = fun_Merton_model_mkt_gap(par,opt_params,Strikes, TTMs,mkt_prices);
Error = norm(err_vect)
plot(Strikes, model_prices,'s'); 
legend('Market price', 'Merton price')