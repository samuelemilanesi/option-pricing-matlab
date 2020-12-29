function V = char_exponent_Kou(v, par)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns the characteristic exponent of logprices in Kou model
    %   INPUT: - par = parameters structure{
    %                  - par.sigma = BM vol 
    %                  - par.p = prob of positive jump
    %                  - par.lambdap = pos jump intensity
    %                  - par.lambdam = neg jump intensity
    %                  - par.lambdaK = jump time intensity
    %                  }
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% de-struct params 
    sigma = par.sigma;
    p = par.p;
    lambdap = par.lambdap;
    lambdam = par.lambdam;
    lambdaK  = par.lambdaK;
    
    % evaluate char exp under risk-neutral measure 
    V = @(v)-sigma^2/2*v.^2+1i*v.*lambdaK.*(p./(lambdap-1i*v)-(1-p)./(lambdam+1i*v)); % without drift
    drift_rn = -V(-1i); % Drift Risk_neutral
    V = drift_rn*1i*v+V(v);
end