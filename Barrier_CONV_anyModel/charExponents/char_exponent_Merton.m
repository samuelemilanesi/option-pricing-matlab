function V = char_exponent_Merton(v, par)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns the characteristic exponent of logprices in Merton model
    %   INPUT: - par = parameters structure{
    %                  - par.r = risk free rate
    %                  - par.sigma = BM vol 
    %                  - par.mu = jump drift
    %                  - par.delta = jump vol
    %                  - par.lambdaK = jump time intensity
    %                  }
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% de-struct params 
    sigma = par.sigma;
    mu = par.mu;
    delta = par.delta;
    lambdaK  = par.lambdaK;
    
    % evaluate char exp under risk-neutral measure 
    V=@(v)-sigma^2/2*v.^2+lambdaK.*(exp(-delta^2/2*v.^2+mu*1i*v)-1); % without drift
    drift_rn = -V(-1i); % Drift Risk_neutral
    V = drift_rn*1i*v+V(v);
end