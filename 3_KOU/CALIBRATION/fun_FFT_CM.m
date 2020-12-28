function [call_prices, put_prices] = fun_FFT_CM(par, Strikes, CharFunc, Npow, A)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Computes EU call and put prices on a grid of stikes under any model with known characteristic
    %   function, using Carr&Madan procedure via FFT
    %
    %   INPUT: - par = parameters struct, different for every model but with at least{
    %                  - par.S0 = underlying  spot price
    %                  - par.r = risk free rate
    %                  - par.TTM = time to maturity
    %                  }
    %           - Strikes = vector of strikes on which compute options' prices
    %           - CharFunct = handle function with the characteristic function of model logprice
    %           - Npow = log_2(number of discretization points in frequencies' domain)
    %
    %   OUTPUT: - call_prices: vector of EU call option prices with strike from Strikes input
    %           - put_prices: vector of EU put option prices with strike from Strikes input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% de-struct params 
    S0 = par.S0;
    r = par.r;
    T = par.TTM;
    
    % Discretization parameter
    N = 2^Npow;

    % v-> compute integral as a summation
    eta = A / N; v = [0:eta:A * (N - 1) / N]; v(1) = 1e-22;
    % lambda-> compute summation via FFT
    lambda = 2 * pi / (N * eta); k = -lambda * N / 2 + lambda * (0:N - 1);              % logstrike grid

    % Fourier transform of z_k
    Z_k = exp(1i * r * v * T) .* (CharFunc(v - 1i) - 1) ./ (1i * v .* (1i * v + 1));    % Carr-Madan Formula: F-transform of tempered price call

    %% Option Price
    w = ones(1, N); w(1) = 0.5; w(end) = 0.5;           % Grid in Fourier Space
    x = w .* eta .* Z_k .* exp(1i * pi * (0:N - 1));    % Grid in logprices space
    z_k = real(fft(x) / pi);                            % Tempered price call value on logstrike grid 
    C = S0 * (z_k + max(1 - exp(k - r * T), 0));        % Call price value on logstrikes grid 
    K = S0 * exp(k);                                    % Stike grid

    % Remove uninteresting strikes
    index = find(K > 0.1 * S0 & K < 3 * S0);
    C = C(index); K = K(index);
    call_prices = interp1(K, C, Strikes, 'spline');

    %% Compute put price and asymptotic CI via Call-Put Parity
    put_parity = @(call_p) call_p - S0 + Strikes.*exp(-r*T);

    put_prices = put_parity(call_prices);
end