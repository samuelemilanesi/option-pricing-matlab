function payoff = payoff_EU_Call(ST, Strike)
    payoff = max(ST - Strike, 0);
end
