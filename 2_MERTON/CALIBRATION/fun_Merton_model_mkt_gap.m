function [gaps, model_prices] =  fun_Merton_model_mkt_gap(par, x, Strikes, TTMs, mkt_prices);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Returns a vector of gaps between Merton model prices and mkt prices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Discretization parameters
    Npow = 15; A = 1200;
    % par struct
    par.sigma = x(1);
    par.mu = x(2);
    par.delta = x(3);
    par.lambdaK = x(4);
    % calibration on call options
    model_prices = zeros(size(TTMs));
    for T = unique(TTMs)'
        idx = TTMs==T;
        par.TTM = T; 
        % Chararacteristic function
        CharFunc = @(v) exp(T*char_exponent_Merton(v,par));
        % model prices evaluation: EU CALL options 
        model_prices(idx) = fun_FFT_CM(par, Strikes(idx), CharFunc, Npow, A);

        % To calibrate using EU PUT OPTIONS USE THIS
        % [~, put_model_prices] = fun_FFT_CM(par, Strikes(idx), CharFunc, Npow, A);
        % model_prices(idx) = put_model_prices
        
    end
    gaps = abs(model_prices - mkt_prices);
end