function codec=codec(M,N,params)  % M,N reversed from convention M is signal size, N is number of measurements

    %y= Phi * x 
    %y= Phi * Psi * s 
    %y= Phi * Psi * s 
    %y= Theta * s 
    %Theta = Phi * Psi
    %Theta' = Psi'*Phi' 
    % x= Psi * s 

    Psi = speye(M);
    Af=@(x) FastCSOperator(1,N,M,Psi*x,true(1,M),M);
    At=@(y) Psi'*FastCSOperator(2,N,M,y,true(1,M),M);
    Theta = opFunction(N,M,{ Af, At });
    
    codec.encode=@encode;
    codec.decode=@decode;
    
    function y=encode(s)
       % x is in the image domain (not wavelet domain)
        
        y = Theta*s;
    end
    function [s_hat,results]=decode(y)
    
        
        x0 = At(y);
        %x_hat = tveq_logbarrier(x0, Af, At, y, 1e-3, 5, 1e-8, 200);
        x_hat = tv_l1magic.optimization.tveq_logbarrier(x0, Af, At, y, 1e-1, 2, 1e-8, 600);
%         x_hat = tv_l1magic.optimization.l1qc_logbarrier(x0, A, [], y, epsilon, 1e-3);
        s_hat = Psi' * x_hat;
        results=[];%samples;
    end

end
