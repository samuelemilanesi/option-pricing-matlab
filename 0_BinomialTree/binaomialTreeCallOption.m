
%% CRR (1979) binomial tree

clear all
close all
clc
%% European Call Options

%1. parameters
rf = 0.01; S0 = 100; sigma = 0.4; %market parameter
T = 1; K = 100; %contract parameter
M = 100; %discr. parameter

%2. Tree generation
dt = T/M;
u = exp(sigma*sqrt(dt)); d = 1/u;
q = (exp(rf*dt)-d) / (u-d);
Tree = zeros(M+1,1);

finalPayoffs = zeros(M+1,1);
finalPayoffs(:,1) = max(S0 * d.^(0:M) .* u.^(M:-1:0) - K, 0);

Tree = finalPayoffs;

for i = M:-1:1
    Tree = exp(-rf*dt) * (q * Tree(1:i) + (1-q) * Tree(2:i+1));
end

Price = Tree(1,1)
Price_BS = blsprice(S0,K,rf,T,sigma,0)

