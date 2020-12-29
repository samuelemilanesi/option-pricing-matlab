function I=JD_integral(nu,x,V,ynodes,lb_cond_t,ub_cond_t)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Integral part of the numerical approximation of option pricess via PIDE on
    %   logprices in Lévy framework for Jump Diffusion processes
    %   INPUT:
    %           - nu = handle, Lévy measure 
    %           - x  = logprices grid 
    %           - V  = last evaluation of the option price
    %           - ynodes = integration nodes
    %           - S0 = spot price
    %           - discK = discounted strike 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dy=ynodes(2)-ynodes(1);
    w=ones(size(ynodes)); w(1)=1/2; w(end)=1/2; % trapez schema
    I=zeros(size(x));
    for i=2:length(I)-1
        I(i)=sum( w.*valueV(x(i)+ynodes,x,V,lb_cond_t,ub_cond_t).*nu(ynodes) )*dy;
    end
end

function f=valueV(y,x,V,lb_cond_t,ub_cond_t)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Helper: computes the value of V on the grid y adjusting boundary conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f=zeros(size(y));

    index=find( (y>=x(1)) .* (y<=x(end)) );
    f(index)=interp1(x,V,y(index));

    % upper boundary condition for x>xmax 
    index=find( (y>x(end)) );     
    f(index) = ub_cond_t(y(index)); 
    % lower boundary condition for x<xmin
    index=find( (y<x(1)) );     
    f(index) = lb_cond_t(y(index)); 
end    
       
