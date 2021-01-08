function V = char_exponent_NIG(v, par)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns the characteristic exponent of logprices in VG model
    %   INPUT: - par = parameters structure{
    %                  - par.sigma 
    %                  - par.theta 
    %                  - par.k
    %                  }
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% de-struct params 
    sigma = par.sigma;
    theta = par.theta;
    k     = par.kNIG;
    
    % evaluate char exp under risk-neutral measure 
    
    V=@(v) (1-sqrt(1+v.^2.*sigma^2*k-1i*theta*v*k))/k; 
    drift_rn=-V(-1i); % Drift Risk_neutral
    V=drift_rn*1i*v+V(v);
end