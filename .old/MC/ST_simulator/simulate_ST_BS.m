function ST = simulate_ST_BS(r, S0, T, Strike, sigma, Nsim, seed)
    % Simulates the final value ST of the underlying under the B&S model 
    % Returns a vector ST of Nsim simulations
    
    if nargin < 7
	seed = 0
    end

    rng(seed);
    Z = randn(Nsim,1);
    ST = S0*exp( (r - sigma^2/2)*T + sigma*sqrt(T)*Z );

end
