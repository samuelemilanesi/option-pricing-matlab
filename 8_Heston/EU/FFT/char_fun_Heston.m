function phi = char_fun_Heston(u, T, par)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns the characteristic function of logprices in Heston model
    %   INPUT:  - v   = evaluation point 
    %           - T   = evaluation time
    %           - par = parameters structure{
    %                  - par.epsilon = vol-of-vol
    %                  - par.kappa   = mean reversion speed
    %                  - par.rho     = correlation
    %                  - par.V0      = initial var
    %                  - par.theta   = mean
    %                  }
    %           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% de-struct params 
    epsilon = par.epsilon; % vol-of-vol
    kappa = par.kappa; % mean reversion speed
    rho = par.rho; % correlation
    theta = par.theta; % mean 
    V0 = par.V0;
    
    alfa = -.5*(u.*u + u*1i);
    beta = kappa - rho*epsilon*u*1i;
    gamma = .5 * epsilon^2;
    
    D = sqrt(beta .* beta - 4.0 * alfa .* gamma);
    
    bD = beta - D;
    eDt = exp(- D * T);
    
    G = bD ./ (beta + D);
    B = (bD ./ epsilon^2) .* ((1.0 - eDt) ./ (1.0 - G .* eDt));
    psi = (G .* eDt - 1.0) ./(G - 1.0);
    A = ((kappa * theta) / (epsilon^2)) * (bD * T - 2.0 * log(psi));
    
    y = A + B*V0;

    phi = exp(y);
end
