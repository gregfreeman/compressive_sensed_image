function encoder_decoder=codec_tv_tfocs(N,M,params)

    global Pstate Qstate
    if ~isfield(params,'tv_lambda')
        params.tv_lambda=0.1;
    end
    if ~isfield(params,'Pstate')
        params.Pstate = 4972169;
    end
    if ~isfield(params,'Qstate')
        params.Qstate = 7256157;
    end

    n1=params.imageSide;
    n2=params.imageSide;
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
   
    lambda=params.tv_lambda;
    vec = @(x) x(:);
    encoder_decoder.encode=@encode;
    encoder_decoder.decode=@decode;

    
    function y=encode(x)
    y  = A_measure*[x;pad];
    end
    function [x_hat,results]=decode(y,coarse)
        % y are the fine compressive measurement 
        % coarse are the coarse direct measurements
        
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

%         Af=@(x) FastCSOperator(1,N,M,Psi*x,true(1,M),M);
%         At=@(y) Psi'*FastCSOperator(2,N,M,y,true(1,M),M);
        normTV           = linop_TV( [n1,n2], [], 'norm' );
        W_tv   = linop_TV( [n1,n2] );
        
        opts = [];
        opts.restart = 100; 
        opts.printEvery = 20;
        opts.tol        = 1e-4;
        opts.normW2     = normTV^2;
        
        contOpts            = [];
        contOpts.maxIts     = 4;
        contOpts.betaTol    = 2;
        x0=A_measure'*y;
        EPS= 1e-2; %.9*norm(b-b_original);
        z0  = [];
        mu  = 1e-3*255;
        % [x,odata,opts] = tfocs( smooth_quad, { A, -b }, prox_l1( lambda ), x0, opts );
        solver = @(mu,x0,z0,opts) solver_sBPDN_W( A_measure, W_tv, y, EPS, mu, x0, z0, opts);
        [ x, odata, opts ] = continuation(solver,5*mu,x0,z0,opts, contOpts);

        results.odata=odata;
        results.opt=opts;
        x_hat=x(1:Nfinal);

    end

end