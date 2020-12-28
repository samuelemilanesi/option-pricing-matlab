clear
close
% PRICE AN EUROPEAN CALL OPTION
% FD SCHEME - IMPLICIT EULER
% LOG-PRICE TRANSFORM
% FINITE ACTIVITY LEVY - KOU MODEL

%% Input
T=1; K=100; S0=100; r=0.01; sigma=0.3;
p=0.5; lambda=2; lambdap=15; lambdam=18;
nu=@(y) lambda*(p*lambdap*exp(-lambdap*y).*(y>0)+...
               (1-p)*lambdam*exp(-lambdam*abs(y)).*(y<0));
M=75; N=2000;

%% Grids
dt=T/M; 
% xmin=1.5*((r-sigma^2/2)*T-6*sigma*sqrt(T));
% xmax=1.5*((r-sigma^2/2)*T+6*sigma*sqrt(T));
Smin=0.2*S0; Smax=3*S0; xmin=log(Smin/S0); xmax=log(Smax/S0);
dx=(xmax-xmin)/N;
x=linspace(xmin, xmax, N+1)';

%% Computing alpha, lambda, and truncating the integral
% Truncate the integral
tol=1e-14;
ymin=-0.1;
while nu(ymin)>tol
    ymin=ymin-0.5;
end
ymax=0.1;
while nu(ymax)>tol
    ymax=ymax+0.5;
end
ynodes=linspace(ymin,ymax,2*N);
figure
plot(ynodes,nu(ynodes));
alpha=trapz(ynodes, (exp(ynodes)-1).*nu(ynodes))
lambdaNum=trapz(ynodes, nu(ynodes))

%% Construct the matrix
A=-(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2);
B=-1/dt-sigma^2/(dx^2)-(r+lambda);
C=(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2);
% 1nd way
Mhat=spdiags( [A*ones(N-1,1) B*ones(N-1,1) C*ones(N-1,1)],[-1 0 1],N-1,N-1);
% 2nd way
Mat=sparse(N+1,N+1); Mat(2:end-1,2:end-1)=Mhat; clear Mhat
Mat(1,1)=1; Mat(2,1)=A; Mat(end,end)=1; Mat(end-1,end)=C;

%% Backward-in-time procedure
V=max( S0*exp(x)-K,0); 
for j=M-1:-1:0
    I=integral(nu,x,V,ynodes,S0,K*exp(-r*(T-(j+1)*dt)));
    rhs=-1/dt*V-I; rhs(1)=0; rhs(end)=S0*exp(xmax)-K*exp(-r*(T-j*dt));
    V=Mat\rhs;
end
figure
plot(S0*exp(x),V); title('Price'); xlabel('S - spot price');
Price=interp1( S0*exp(x),V, S0,'spline')
Price_CM=FFT_CM_Call

function I=integral(nu,x,V,ynodes,S0,discK)
% trapezoidal formula
dy=ynodes(2)-ynodes(1);
w=ones(size(ynodes)); w(1)=1/2; w(end)=1/2;
I=zeros(size(x));
for i=2:length(I)-1
    I(i)=sum( w.*valueV(x(i)+ynodes,x,V,S0,discK).*nu(ynodes) )*dy;
end
end

function f=valueV(y,x,V,S0,discK)
f=zeros(size(y));
index=find( (y>=x(1)) .* (y<=x(end)) );
f(index)=interp1(x,V,y(index));
index=find( (y>x(end)) );
f(index)=S0*exp(y(index))-discK;
end    
    

    
    
    
    
    
    
    
    
    
    


