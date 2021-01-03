function V = char_exponent_BS(v, par)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns the characteristic exponent of logprices in BS model
    % INPUT: v = point of evaluation; 
    %        par = parameter struct{
    %              - sigma = volatility of the model
    % OUTPUT: characteristic exponent for the BS model under risk neutral measure 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % de-struct params
    sigma = par.sigma;
    r = par.rf;
    % evaluate char exp under risk-neutral measure 
    V = @(v) - sigma^2/2 * v.^2; % without drift
    drift_rn = -V(-1i); % Drift Risk_neutral
    V = drift_rn * 1i * v + V(v);
end