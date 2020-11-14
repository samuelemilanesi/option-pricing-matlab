function payoff = payoff_EU_Call(ST, Nsim_alpha, Strike)

     
    

    payoff = max(ST - Strike, 0);
end
