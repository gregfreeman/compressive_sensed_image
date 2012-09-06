function encoder_decoder=codec_stomp(N,M,params)


    if ~isfield(params,'lambda')
        params.lambda=0.1;
    end
    if ~isfield(params,'Pstate')
        params.Pstate = 4972169;
    end
    if ~isfield(params,'Qstate')
        params.Qstate = 7256157;
    end
    params.N=N;
    params.M=M;
    if mod(N,2)==1  % force samples to be even
        params.Nfinal=N;
        params.N=N+1;
        params.pad=[0];
    else
        params.Nfinal=N;
        params.N=N;
        params.pad=[];
    end
    

    global Pstate Qstate
    Pstate = params.Pstate;
    Qstate = params.Qstate;
    
%     A_measure = opFunction(M,N,{  @(x) FastCSOperator(1,M,N,x,true(1,N),N), @(x) FastCSOperator(2,M,N,x,true(1,N),N) });
    
   
    encoder_decoder.encode=@encode;
    encoder_decoder.decode=@decode;
    
    function samples=encode(input)

        input = [input(:);params.pad];
% generate the vector b
        samples = FastCSOperator(1,params.M,params.N,input,true(size(input)),params.N);

%         y  = A_measure*x;
    end
    function [output,results]=decode(samples)
   
        stages=10;
        m=length(samples);
        if isfield(params,'StompTheshold') && strcmp(params.StompTheshold,'FDR')  
            if(isfield(params,'StompQ') )  
                q=params.StompQ;
            else
                q=0.5;
            end
            [output, results.iters,results.nActiveSet,results.normRes] = SolveStOMP('FastCSOperator', samples, params.N, 'FDR', q, stages, 1);
            samples_hat = FastCSOperatorWrapper(1,params.M,params.N,output);  
            results.normRes2= norm(samples_hat-samples);
        elseif isfield(params,'StompTheshold') && strcmp(params.StompTheshold,'RCL')  
            [output, results.iters,results.nActiveSet,results.normRes] = SolveStOMP('FastCSOperator', samples, params.N, 'RCL', 0.7, stages, 1);

            samples_hat = FastCSOperatorWrapper(1,params.M,params.N,output);  
            results.normRes2= norm(samples_hat-samples);

        elseif isfield(params,'usepdco')
            [ output, results ] = FastCSExtract( samples, params.N, params );
            samples_hat = FastCSOperatorWrapper(1,m,N,output);  
            results.normRes= norm(samples_hat-samples);


        else
            a_0 = 0.4*m/params.N/stages;
            [output, results.iters,results.nActiveSet,results.normRes] = SolveStOMP('FastCSOperator', samples, params.N, 'FAR', a_0, stages, 1);
            samples_hat = FastCSOperatorWrapper(1,params.M,params.N,output);  
            results.normRes2= norm(samples_hat-samples);

        end
        output =output(1:params.Nfinal);
    end


end