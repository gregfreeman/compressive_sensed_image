function encoder_decoder=codec_spgl1(N,M,params)

    global Pstate Qstate
    if ~isfield(params,'lasso_lambda')
        params.lasso_lambda=0.1;
    end
    if ~isfield(params,'optTol')
        params.optTol=1e-4;
    end
    if ~isfield(params,'residual_constraint')
        params.residual_constraint=sqrt(M)*0.15;
    end
    if ~isfield(params,'Pstate')
        params.Pstate = 4972169;
    end
    if ~isfield(params,'Qstate')
        params.Qstate = 7256157;
    end

    Pstate = params.Pstate;
    Qstate = params.Qstate;
    if mod(N,2)==1  % force samples to be even
        Nfinal=N;
        N=N+1;
        pad=[0];
    else
        Nfinal=N;
        N=N;
        pad=[];
    end
    A_measure = opFunction(M,N,{  @(x) FastCSOperator(1,M,N,x,true(1,N),N), @(x) FastCSOperator(2,M,N,x,true(1,N),N) });
        
    encoder_decoder.encode=@encode;
    encoder_decoder.decode=@decode;
    encoder_decoder.Phi=A_measure;
    
    function y=encode(x)
    y  = A_measure*[x;pad];
    end
    function [x_hat,results]=decode(y)
        opts = spgSetParms('optTol', params.optTol, 'verbosity', 1);
    
        [x_hat,r,g,info] = spg_bpdn(A_measure,y,params.residual_constraint,opts);
        results.r=r;
        results.g=g;
        results.info=info;
        x_hat=x_hat(1:Nfinal);
    end

end