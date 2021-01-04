%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options in under Heston model using 
% Carr and Madan Formula via FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
% Market parameters
r = 0.001;               % riskfree interest rate 
S0 = 100;            % spot price
% Model parameters
epsilon = 0.2;% vol-of-vol
k = 0.01;     % mean reversion speed
rho = -0.2;   % correlation
theta = 0.1;  % mean 
V0 = 0.1; 
% Contract parameters
T = 1;                  % maturity
Strike = [S0; S0*0.9];                 % strike

par = struct('S0',S0,'r',r,'TTM',T,'epsilon',epsilon,'kappa',k,'rho',rho,'theta',theta,'V0',V0);

% Discretization parameter
Npow = 15; A = 600; 

% Chararacteristic function
CharFunc = @(v) char_fun_Heston(v,T,par);

% Carr&Madan method
disp("Vanilla prices using Carr-Madan under Heston model:")
[call_price, put_price] = fun_FFT_CM(par, Strike, CharFunc, Npow,A)