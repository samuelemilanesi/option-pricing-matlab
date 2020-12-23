%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calibrate BS model parameter on a set of vanilla options 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

%% Marked to market EU Call options data

% vector of strikes of each traded option
Strikes = [200 225 210 180 215 230 220 175 150 170 160 195 230 210 220 225 200 215 205 235 ]';
% vector of time to maturity of each traded option
TTMs = [0.3333 0.3333 0.3333 0.3333 0.3333 0.3333 0.3333 0.3333 0.3333 0.3333 0.3333 0.3333 0.0159 0.0159 0.0159 0.0159 0.0159 0.0159 0.0159 0.0159 ]';
% vector of market prices
mkt_prices = [25.3000 10.0000 18.2000 41.8500 15.2000 7.8500 12.3600 46.0500 70.7000 50.8500 60.4000 29.4000 0.2400 9.2500 2.2600 0.7600 19.0300 5.2000 14.0000 0.1000]';
S0 = 218.75;    % spot price of the underlying asset 
r=0.02;         % risk free rate 

%% RSS loss minimization  
gaps = @(vol) BS_model_mkt_gap(S0, Strikes, r, TTMs, vol, mkt_prices);

sigma = lsqnonlin(@(vol) gaps(vol), 0.2, 0.01, 0.8)    % 0.2 = ansaz; 0.01 = LB; 0.8 = UB 


figure
plot(Strikes, mkt_prices,'+'); hold on
model_prices = blsprice(S0, Strikes, r, TTMs, sigma);
plot(Strikes, model_prices,'s'); 
legend('Market price', 'BS price')