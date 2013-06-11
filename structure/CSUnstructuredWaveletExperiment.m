function [outputImage,results]=CSUnstructuredWaveletExperiment(inputImage,params)

params.imageSide=size(inputImage,1);
if ~isfield(params,'delta')
    params.delta=0.1;
end
if ~isfield(params,'coarsestScale')
    params.coarsestScale=4;
end

if ~isfield(params,'reconstruct')
    params.reconstruct = 'stomp';
end

if ~isfield(params,'qmf')
    params.qmf = MakeONFilter('Symmlet',8);
end

[n1,n2]=size(inputImage);
if(n1~=n2)
    error('expect image with equal dimensions for each side')
end
n1log2=log2(n1);
if(n1log2~=floor(n1log2))
    error('expect image a size of 2 power')
end
    
delta = params.delta;
N=n1*n2;
M=floor(delta * N);
levels = n1log2-params.coarsestScale;

vec = @(x) x(:);
mat = @(x) reshape(x,n1,n2);
W_op=opWavelet2(n1,n2,'Custom',params.qmf,levels);


quantizeOptions=params;
quantizeOptions.Mcoarse=0;
quantizer=sampleQuantizer(quantizeOptions);


switch(params.reconstruct)
case 'stomp'
   make_coder=@(N,M,parms) codec_stomp(N,M,parms);
case 'lasso_tfocs'
   make_coder=@(N,M,parms) codec_lasso_tfocs(N,M,parms);
case 'tswcs'
   params.levels=levels;
   params.n1=n1;
   params.n2=n2;
   make_coder=@(N,M,parms) tswcs.codec(N,M,parms);
case 'bcs_rvm'
   make_coder=@(N,M,parms) bcs_rvm.codec(N,M,parms);
case 'bcs_gf'
   params.vfine=vfine;
   params.levels=levels;
   params.n1=n1;
   params.n2=n2;
   make_coder=@(N,M,parms) bcs_gf.codec(N,M,parms);
case 'cosamp'
   make_coder=@(N,M,parms) cosamp.codec(N,M,parms);
case 'tv_l1magic'
   params.sparseBasis=W_op';
   make_coder=@(N,M,parms) tv_l1magic.codec(N,M,parms);   
case 'tv_tfocs'
   params.sparseBasis=W_op';
   make_coder=@(N,M,parms) codec_tv_tfocs(N,M,parms);   
case 'spgl1'
   make_coder=@(N,M,parms) codec_spgl1(N,M,parms);

end
coder=make_coder(N,M,params);   
  
samples=CompressHere(inputImage,params);
[outputImage results]=ExtractHere(samples,params);
results.nSamples=length(samples);

if(isfield(params,'save_wavelet_data') && params.save_wavelet_data~=0)
    results.W_op = W_op;
    results.mat = mat;
end
if(isfield(params,'save_coder_data') && params.save_coder_data~=0)
    results.coder = coder;
    results.samples = samples;    
end


function [samples] = CompressHere(inputImage,params)
    disp(params)
 

    theta=W_op*vec(inputImage);   
    samples= coder.encode(theta);
    samples=quantizer.encode(samples);
    if length(samples)~=M
        warning('There was an error counting the samples')
    end

end


function [outputImage,results] = ExtractHere(samples,params)

    disp(params)

    samples=quantizer.decode(samples);
    results=struct();


    if length(samples)~=M
        warning('There was an error in the number of samples when reconstructing the image')
    end

    % Solve the CS problem 
    [w_hat,results.decoder] = coder.decode(samples);   

    % Reconstruct image
    outputImage = mat(W_op'*w_hat);

end

end
