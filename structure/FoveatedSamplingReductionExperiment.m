function [outputImage,results]=FoveatedSamplingReductionExperiment(inputImage,params)

params.imageSide=size(inputImage,1);
if ~isfield(params,'delta')
    params.delta=0.1;
end
if ~isfield(params,'coarsestScale')
    params.coarsestScale=4;
end

if ~isfield(params,'Pstate')
    params.Pstate = 4972169;
end
if ~isfield(params,'Qstate')
    params.Qstate = 7256157;
end

if ~isfield(params,'qmf')
    params.qmf = MakeONFilter('Symmlet',8);
end

if ~isfield(params,'viewDist')
    params.viewDist = 3; % view distance as a multiple of screen size
end

if ~isfield(params,'reconstruct')
    params.reconstruct = 'stomp';
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
W_op=opWavelet2(n1,n2,'Custom',params.qmf,levels);

weight=foveateWaveletSensitivity([n1,n2],params.fovea, levels, params.viewDist);
mask=weight>=0.01;
fine=true(n1,n2);
fine(1:2^j0,1:2^j0)=false;
mask(1:2^j0,1:2^j0)=false;
vfine=vec(fine);
vmask=vec(mask);

Ncoarse = sum(~vfine);
Nfine= sum(vmask);
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
end
coder=make_coder(Nfine,Mfine,params);   

samples=CompressHere(inputImage,params);
[outputImage results]=ExtractHere(samples,params);
results.nSamples=length(samples);


function [samples] = CompressHere(inputImage,params)
    disp(params)

    alpha0=W_op*vec(inputImage);   
    samples_coarse = alpha0(~vfine);
    theta = alpha0(vmask);
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
    alpha_CS(vmask)=alpha_fine;

    % Reconstruct
    outputImage = mat(W_op'*alpha_CS);
    
end

end
