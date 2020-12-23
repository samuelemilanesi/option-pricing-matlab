function gaps =  BS_model_mkt_gap(S0, Strikes, r, TTMs, vol, mkt_prices)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Returns a vector of gaps between BS model prices and mkt prices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model_prices = blsprice(S0,Strikes,r,TTMs,vol);
    gaps = model_prices - mkt_prices;
end