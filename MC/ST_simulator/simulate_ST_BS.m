function ST = simulate_ST_BS(r, S0, T, Strike, sigma, Nsim, seed=0)
    % Simulates the final value ST of the underlying under the B&S model 
    % Returns a vector ST of Nsim simulations
                             
    Z = randn(Nsim,1);
    ST = S0*exp( (r - sigma^2/2)*T + sigma*sqrt(T)*Z )

    return ST;
end
