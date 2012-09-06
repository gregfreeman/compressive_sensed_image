function encoder_decoder=codec(M,N,params)  % M,N reversed from convention M is signal size, N is number of measurements
   

    Phi = randn(N,M);
    Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[N,1]);	
    % Phi = opFunction(N,M,{  @(x) FastCSOperator(1,N,M,x,true(1,M),M), @(x) FastCSOperator(2,N,M,x,true(1,M),M) });

    encoder_decoder.encode=@encode;
    encoder_decoder.decode=@decode;
    
    function y=encode(x)
        y = Phi*x;
    end
    function [x_hat,results]=decode(y)
        initsigma2 = std(y)^2/1e2;
        eta=1e-8;
        [weights,index,sigma2,errbars,basis] = bcs_rvm.BCS_fast_rvm(Phi,y,initsigma2,eta);
        x_hat=zeros(M,1);
        x_hat(index)=weights;
        results.errbars=errbars;
        results.sigma2=sigma2;
    end

end