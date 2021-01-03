function F = charfunction(u,parameters,flag)
    % flag=0 --> funzione caratteristica per problema backward 
    % flag=1 --> funzione caratteristica per problema forward
    if nargin==2
        flag=0;
    end
    
    dt = parameters.dt;
    
    switch parameters.distr
    case 1 % Normal
        F = exp(dt*char_exponent_BS(u,parameters));
    case 2 % Merton
        F = exp(dt*char_exponent_Merton(u,parameters));
    case 3 % Kou
        F = exp(dt*char_exponent_Kou(u,parameters));
    case 4 % NIG
        error('TO BE IMPLEMENTED')
    case 5 % VG 
        F = exp(dt*char_exponent_VG(u,parameters));
    case 6 % ExtVG
        F = exp(dt*char_exponent_extVG(u,parameters));
    case 7 % ExtNIG
        error('TO BE IMPLEMENTED')
    case 8 % Heston 
        error('TO BE IMPLEMENTED')
    end
    
    if flag==0
        F = F.*exp(parameters.rf.*1i.*u*dt);
        F=conj(F);
    end
   
    
    