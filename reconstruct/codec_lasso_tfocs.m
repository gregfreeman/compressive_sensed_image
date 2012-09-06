function encoder_decoder=codec_lasso_tfocs(N,M,params)

    global Pstate Qstate
    if ~isfield(params,'lasso_lambda')
        params.lasso_lambda=0.1;
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
    

    lambda =params.lasso_lambda;

    
    encoder_decoder.encode=@encode;
    encoder_decoder.decode=@decode;
    encoder_decoder.Phi=A_measure;
    
    function y=encode(x)
    y  = A_measure*[x;pad];
    end
    function [x_hat,results]=decode(y)
    
        %% Call the TFOCS solver

        % mu              = 5;
        % er              = @(s_hat) norm(A_measure*s_hat-y) ;
        % opts = [];
        % opts.errFcn     = @(f,dual,primal) er(primal);
        % opts.maxIts     = 100;
        % opts.printEvery = 20;
        % opts.tol        = 1e-4;


        % build operators:
        %A           = linop_handles([N,N], @(x)x, @(x) x );


        %normA2      = 1;
        %W_measure   = linop_spot( A_measure );
        % W_wavelet   = linop_spot( opWavelet(n1,n2,'Daubechies') );
        % W_curvelet  = linop_spot( opCurvelet(n1,n2) );
        % W_tv        = linop_TV( [n1,n2] );
        % normWavelet      = linop_normest( W_wavelet );
        % normCurvelet     = linop_normest( W_curvelet );
        % normTV           = linop_TV( [n1,n2], [], 'norm' );

        % contOpts            = [];
        % contOpts.maxIts     = 4;


        % [ x, odata, opts ] = solver_L1RLS( A, b, lambda, x0, opts )
        %    Solves the l1-regularized least squares problem
        %        minimize (1/2)*norm( A * x - b )^2 + lambda * norm( x, 1 )
        %    using the Auslender/Teboulle variant with restart. A must be a matrix

        opts = [];
        opts.restart = 100; 
        opts.printEvery = 20;
        opts.tol        = 1e-4;
   
        x0=zeros(N,1);
%         x0=A_measure'*y;
        % [x,odata,opts] = tfocs( smooth_quad, { A, -b }, prox_l1( lambda ), x0, opts );
        [x_hat,odata,opts] = tfocs( smooth_quad, { linop_spot( A_measure ), -y }, prox_l1( lambda ), x0, opts );
        results.odata=odata;
        results.opt=opts;
        x_hat=x_hat(1:Nfinal);

    end

end