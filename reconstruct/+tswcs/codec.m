function encoder_decoder=codec(M,N,params)  % M,N reversed from convention M is signal size, N is number of measurements
   

    Phi = randn(N,M);
    Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[N,1]);	
    % Phi = opFunction(N,M,{  @(x) FastCSOperator(1,N,M,x,true(1,M),M), @(x) FastCSOperator(2,N,M,x,true(1,M),M) });

    encoder_decoder.encode=@encode;
    encoder_decoder.decode=@decode;
    save_mcmc_iterations=0;
    if(isfield(params,'save_mcmc_iterations') && params.save_mcmc_iterations~=0)
        save_mcmc_iterations=1;
    end
    
    function y=encode(x)
        y = Phi*x;
    end
    function [x_hat,results]=decode(y)
    
        [IdxParent, IdxChildren, Ms]=tswcs.WaveRelation2D_opWavelet([params.n1 params.n2], params.levels);
        idx=1:(params.n1*params.n2);
        idxFine=idx(params.vfine);
        IdxParent=tswcs.mapIndexes(idxFine,IdxParent);
        IdxChildren2=zeros(4,length(IdxParent));
        for i=1:4
            IdxChildren2(i,:)=tswcs.mapIndexes(idxFine,IdxChildren(i,:));
        end
        % indicator of level for each wavelet coefficient
        wLevel=zeros(params.n1,params.n2);
        L=length(Ms);
        n1log2=log2(params.n1);
        for s=L:-1:1
            k=n1log2-L+s;
            wLevel(1:2^k,1:2^k)=s;
        end
        wLevel=reshape(wLevel,params.n1*params.n2,1);
        wLevel_fine=wLevel(params.vfine);
 
        [x_hat, samples,burnin_samples] = tswcs.TSWCSmcmc_NoSc(Phi, y, Ms, IdxParent', IdxChildren2',wLevel_fine);
        if save_mcmc_iterations
            results.iterations=samples;
            results.burnin_samples=burnin_samples;
        else
            results=[];
        end
    end

end