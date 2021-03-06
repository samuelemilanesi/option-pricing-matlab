function S = BS_simulate(S0, r, sigma, T, Nsim, M)
    % Simulate the evolution of a process S
    % dS = (r) S dt + sigma S dW(t), S(0)=S0.

    S = zeros(Nsim, M + 1); S(:, 1) = S0;
    dt = T / M;
    for i = 1:M
        S(:, i + 1) = S(:, i) .* exp((r - sigma^2/2) * dt + sigma * sqrt(dt) * randn(Nsim, 1));
    end
end