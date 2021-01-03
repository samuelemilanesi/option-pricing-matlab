function [S,v] = CONV( S_0, K, Ndate, N, Barriera, param)

    b=2.5;
    [x,~,~,H] = kernel(N,-b,b,param,0); % logprice grid x, conjugate of the characteristic function
    S = S_0*exp(x);
    % Payoff e trasformata di Fourier della density
    v = max(S-K,0).*(S>Barriera); %DOWN AND OUT CALL
    H=ifftshift(H);
    for j = 1:Ndate
        %v(S<=Barriera) = 0; % se in t=0 la barriera NON ha effetto
        %v=real((fft(ifft(v).*H)))*exp(-param.rf*param.dt);
        v=real(fftshift(fft( ...
                              ifft(ifftshift(v)) ...
                                .*H ...
                            ) ...
                ) ...
            )*exp(-param.rf*param.dt);
        v(S<=Barriera) = 0; % se in t=0 la barriera ha effetto
    end
    index=find( (S>0.1*S_0).*(S<3*S_0));
    S=S(index); v=v(index);
    figure
    plot(S,v)

end

