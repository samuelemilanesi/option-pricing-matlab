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
Strike = S0;            % strike in which to evaluate the option
% Discretization parameter
Npow = 20; N = 2^Npow; A = 1200;

% v-> compute integral as a summation
eta = A / N; v = [0:eta:A * (N - 1) / N]; v(1) = 1e-22;
% lambda-> compute summation via FFT
lambda = 2 * pi / (N * eta); k = -lambda * N / 2 + lambda * (0:N - 1);              % logstrike grid

% Fourier transform of z_k
CharFunc = @(v) exp(T * char_exponent_BS(v,sigma));                                 % Characteristic function of X_T
Z_k = exp(1i * r * v * T) .* (CharFunc(v - 1i) - 1) ./ (1i * v .* (1i * v + 1));    % Carr-Madan Formula: F-transform of tempered price call

%% Option Price
w = ones(1, N); w(1) = 0.5; w(end) = 0.5;           % Grid in Fourier Space
x = w .* eta .* Z_k .* exp(1i * pi * (0:N - 1));    % Grid in logprices space
z_k = real(fft(x) / pi);                            % Tempered price call value on logstrike grid 
C = S0 * (z_k + max(1 - exp(k - r * T), 0));        % Call price value on logstrikes grid 
K = S0 * exp(k);                                    % Stike grid

% Remove uninteresting strikes and plot
index = find(K > 0.1 * S0 & K < 3 * S0);
C = C(index); K = K(index);
plot(K, C)
title('Option Price');
xlabel('Strike');

disp("Vanilla prices using Carr-Madan under BS model:")
call_price = interp1(K, C, Strike, 'spline')

%% Compute put price and asymptotic CI via Call-Put Parity
put_parity = @(call_p) call_p - S0 + Strike*exp(-r*T);

put_price = put_parity(call_price)
