%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price European Call and Put Options in the B&S model using 
% Carr and Madan Formula via FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
% Market parameters
r = 0.02;               % riskfree interest rate 
S0 = 218.75;            % spot price
% Model parameters
sigma = 0.2516;         % standard deviation 
% Contract parameters% logstrike grid
T = 1;                  % maturity
Strike = [S0,S0*0.9];            % strike in which to evaluate the option

par = struct('S0',S0,'r',r,'TTM',T,'sigma',sigma);

% Discretization parameter
Npow = 20; A = 1200;

% Chararacteristic function
CharFunc = @(v) exp(T* char_exponent_BS(v,par) );

% Carr&Madan method
disp("Vanilla prices using Carr-Madan under BS model:")
[call_price, put_price] = fun_FFT_CM(par, Strike, CharFunc, Npow,A)


disp('----BS closed formula----:')
[bs_call, bs_put] = blsprice(S0,Strike,r,T,sigma)