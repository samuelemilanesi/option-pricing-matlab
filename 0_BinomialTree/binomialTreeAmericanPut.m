
%% CRR (1979) binomial tree

clear all
close all
clc

%% American Put Option

%1. parameters
rf = 0.01; S0 = 100; sigma = 0.4; %market parameter
T = 1; K = 100; %contract parameter
M = 1000; %discr. parameter

%2. Tree generation
dt = T/M;
u = exp(sigma*sqrt(dt)); d = 1/u;
q = (exp(rf*dt)-d) / (u-d);

finalPayoffs = zeros(M+1,1);
finalPayoffs(:,1) = max(K - S0 * d.^(0:M) .* u.^(M:-1:0), 0);

Tree = finalPayoffs;

for i = M:-1:1
    earlyExerc = max(K - S0 * d.^(0:i-1) .* u.^(i-1:-1:0), 0)';
    continueContr = exp(-rf*dt) * (q * Tree(1:i) + (1-q) * Tree(2:i+1));
    Tree = max(earlyExerc, continueContr);
end

PriceAm = Tree(1,1)
PriceEU = blsprice(S0,K,rf,T,sigma,0) - S0 + K * exp(-rf*T)


