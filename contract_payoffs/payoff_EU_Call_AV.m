function payoff = payoff_EU_Call_AV(ST, ST_AV, Strike)

    payoff = max( (ST+ST_AV)/2 - Strike, 0);

end
