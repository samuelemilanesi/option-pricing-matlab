function [S,ST]=Heston_simulate_asset(par,Nsim,M)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   INPUT: - par = parameters structure{
    %                  - par.S0 = underlying  spot price
    %                  - par.r = risk free rate
    %                  - par.TTM = time to maturity
    %                  - par.epsilon = vol-of-vol
    %                  - par.kappa   = mean reversion speed
    %                  - par.rho     = correlation
    %                  - par.V0      = initial var
    %                  - par.theta   = mean
    %                  }
    %           - Nsim = number of simulations
    %           - M = number of time steps
    %   OUTPUT: - S = matrix of Nsim underlying paths (one for each row)
    %           - ST = vector of Nsim underlying sim at time T 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% de-struct params 
    epsilon = par.epsilon; % vol-of-vol
    k = par.kappa; % mean reversion speed
    rho = par.rho; % correlation
    theta = par.theta; % mean 
    V0 = par.V0;
    T = par.TTM; 
    S0 = par.S0;
    r = par.r;
    dt=T/M; 
    V = zeros(M+1,Nsim); S = zeros(M+1,Nsim);
    
    % Random changes for volatility are sampled from one of two distributions,
    % depending on the ratio Psi = s^2/m^2, where m & s are mean and variance
    % of next volatility value, conditioned to current one.
    % Scheme 1 (exponential approx.) selected when Psi>Psi_cutoff
    % Scheme 2 (quadratic approx.) selected when 0<Psi<Psi_cutoff
    
    % choose Psi_cutoff = 1.5 as suggested in article
    Psi_cutoff = 1.5;
    
    %--- DISCRETIZE V ---------------------------------------------------------
    V(1,:)=V0;
    for i=1:M
        % STEP 1 & 2
        % Calculate m,s,Psi
        m = theta+(V(i,:)-theta)*exp(-k*dt);
        m2 = m.^2;
        s2 = V(i,:)*epsilon^2*exp(-k*dt)*(1-exp(-k*dt))/k+...
            theta*epsilon^2*(1-exp(-k*dt))^2/(2*k);
        s = sqrt(s2);
        
        Psi = (s2)./(m2);
        
        % STEP 3,4,5
        % Depending on Psi, use exp or quad scheme to calculate next V
        index=find(Psi>Psi_cutoff);
        % Exponential approximation, used for Psi>Psi_cutoff)
        % PDF of V(t+dt) is p*delta(0) + (1-p) * (1-exp(-beta x )
        % thus a prob. mass in 0 & an exp tail after that
        p_exp = (Psi(index)-1)./(Psi(index)+1);	% Prob. mass in 0, Eq 29
        beta_exp = (1-p_exp)./m(index);		% exponent of exp. density tail, Eq 30
        % gets x from inverse CDF applied to uniform U
        U = rand(size(index));
        V(i+1,index) = (log((1-p_exp)./(1-U))./beta_exp).*(U>p_exp);
        index=find(Psi<=Psi_cutoff);
        % Quadratic approx, used for 0<Psi<Psi_cutoff
        % V(t+dt) = a(b+Zv)^2, Zv~N(0,1)
        invPsi = 1./Psi(index);
        b2_quad = 2*invPsi-1+sqrt(2*invPsi).*sqrt(2*invPsi-1); % Eq 27
        a_quad = m(index)./(1+b2_quad);	% Eq 28
        V(i+1,index) = a_quad.*(sqrt(b2_quad)+randn(size(index))).^2; % Eq.23
        
    end
    
    %-- DISCRETIZE S----------------------------------------------------------
    S(1,:)=S0;
    gamma1=0.5; gamma2=0.5; %central discretization scheme
    k0=r*dt-rho*k*theta*dt/epsilon;
    k1=gamma1*dt*(k*rho/epsilon-0.5)-rho/epsilon;
    k2=gamma2*dt*(k*rho/epsilon-0.5)+rho/epsilon;
    k3=gamma1*dt*(1-rho^2);
    k4=gamma2*dt*(1-rho^2);
    for i=1:M
        S(i+1,:)=exp(log(S(i,:))+k0+k1*V(i,:)+k2*V(i+1,:)+...
            sqrt(k3*V(i,:)+k4*V(i+1,:)).*randn(1,Nsim));
    end

    S = S';
    ST = S(:,end);
end

    