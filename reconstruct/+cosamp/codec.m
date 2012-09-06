function encoder_decoder=codec(M,N,params)  % M,N reversed from convention M is signal size, N is number of measurements
   
    if ~isfield(params,'cosamp_maxiterations')
        params.cosamp_maxiterations=50;
    end
    if isfield(params,'cosamp_K_factor')
        params.cosamp_K=floor(M*params.cosamp_K_factor);
    end
    if ~isfield(params,'cosamp_K')
        params.cosamp_K=floor(M*.10);
    end
    if ~isfield(params,'cosamp_tol')
        params.cosamp_tol=1e-2;
    end
    if ~isfield(params,'Pstate')
        params.Pstate = 4972169;
    end
    if ~isfield(params,'Qstate')
        params.Qstate = 7256157;
    end
    
    global Pstate Qstate
    Pstate = params.Pstate;
    Qstate = params.Qstate;
%     Phi = randn(N,M);
%     Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[N,1]);	
    Phi = opFunction(N,M,{  @(x) FastCSOperator(1,N,M,x,true(1,M),M), @(x) FastCSOperator(2,N,M,x,true(1,M),M) }); 
   
    encoder_decoder.encode=@encode;
    encoder_decoder.decode=@decode;
    
    function y=encode(x)
        y = Phi*x;
    end
    function [x_hat,results]=decode(y)
        maxiterations=params.cosamp_maxiterations;
        K=params.cosamp_K;
        tol=params.cosamp_tol;
        [x_hat,iterations,error] = cosamp.cosamp(Phi,y,K,tol,maxiterations);
        results.iterations=iterations;
        results.error=error;
    end

end