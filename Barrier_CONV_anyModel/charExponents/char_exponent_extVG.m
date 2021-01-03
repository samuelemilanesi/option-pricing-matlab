function V = char_exponent_extVG(v, par)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns the characteristic exponent of logprices in VG model
    %   INPUT: - par = parameters structure{
    %                  - par.sigma
    %                  - par.sigmaGBM 
    %                  - par.theta 
    %                  - par.k
    %                  }
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% de-struct params 
    sigma = par.sigma;
    sigmaGBM = par.sigmaGBM;
    theta = par.theta;
    k     = par.kVG;
    
    % evaluate char exp under risk-neutral measure 
    V=@(v) -sigmaGBM^2/2 *v.^2 -log(1+v.^2.*sigma^2*k/2-1i*theta*k*v)/k; % without drift
    drift_rn=-V(-1i); % Drift Risk_neutral
    V=drift_rn*1i*v+V(v);
end