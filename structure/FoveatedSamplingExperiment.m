function [outputImage,results]=FoveatedSamplingExperiment(inputImage,params)

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

if ~isfield(params,'foveation_implementation')
    options1=struct();
    options1.viewDist=1;
    options1.numBands=20;
    options1.filterHalfSize=20;
    foveator=foveateFilterBank([n1,n2],params.fovea,options1);
else
    foveator=params.foveation_implementation;
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
j0=params.coarsestScale;

vec = @(x) x(:);
mat = @(x) reshape(x,n1,n2);
W_op=opWavelet(n1,n2,params.qmf,levels);


fine=true(n1,n2);
fine(1:2^j0,1:2^j0)=false;
vfine=vec(fine);


Ncoarse = sum(~vfine);
Nfine= sum(vfine);
Mfine = M - Ncoarse;

quantizeOptions=params;
quantizeOptions.Mcoarse=Ncoarse;
quantizeOptions.Mfine=Mfine;
quantizer=sampleQuantizer(quantizeOptions);

switch(params.reconstruct)
case 'stomp'
   make_coder=@(N,M,parms) codec_stomp(Nfine,Mfine,parms);
case 'lasso_tfocs'
   make_coder=@(N,M,parms) codec_lasso_tfocs(Nfine,Mfine,parms);
case 'tswcs'
   params.vfine=vfine;
   params.levels=levels;
   params.n1=n1;
   params.n2=n2;
   make_coder=@(N,M,parms) tswcs.codec(Nfine,Mfine,parms);
case 'bcs_rvm'
   make_coder=@(N,M,parms) bcs_rvm.codec(Nfine,Mfine,parms);
case 'bcs_gf'
   params.vfine=vfine;
   params.levels=levels;
   params.n1=n1;
   params.n2=n2;
   make_coder=@(N,M,parms) bcs_gf.codec(Nfine,Mfine,parms);
case 'cosamp'
   make_coder=@(N,M,parms) cosamp.codec(Nfine,Mfine,parms);
case 'tv_l1magic'
   params.sparseBasis=W_op';
   make_coder=@(N,M,parms) tv_l1magic.codec(Nfine,Mfine,parms);   
case 'tv_tfocs'
   params.sparseBasis=W_op';
   make_coder=@(N,M,parms) codec_tv_tfocs(Nfine,Mfine,parms);   
end
coder=make_coder(Nfine,Mfine,params);   
 

samples=CompressHere(inputImage,params);
[outputImage results]=ExtractHere(samples,params);
results.nSamples=length(samples);


function [samples] = CompressHere(inputImage,params)
    disp(params)

    foveatedImage=foveator.foveate(inputImage);

    alpha0=W_op*vec(foveatedImage);   
    samples_coarse = alpha0(~vfine);
    theta = alpha0(vfine);
    samples_fine= coder.encode(theta);
    samples=[samples_coarse; samples_fine];
    
    samples=quantizer.encode(samples);
    if length(samples)~=M
        warning('There was an error counting the samples')
    end

end


function [outputImage,results] = ExtractHere(samples,params)
    disp(params)

    samples=quantizer.decode(samples);
     
    results=struct();

    alpha_CS = zeros(N,1);
    alpha_CS(~vfine)=samples(1:Ncoarse);
    start_i=Ncoarse+1;

    fine_samples = samples(start_i:end);

    if length(fine_samples)~=Mfine
        warning('There was an error in the number of samples when reconstructing the image')
    end

    % Solve the CS problem for fine scale

    [alpha_fine,results.fine] = coder.decode(fine_samples);
    alpha_CS(vfine)=alpha_fine;

    % Reconstruct
    outputImage = mat(W_op'*alpha_CS);

end

end
