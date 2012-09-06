function [outputImage,results]=ImageDomainExperiment(inputImage,params)

params.imageSide=size(inputImage,1);
if ~isfield(params,'delta')
    params.delta=0.1;
end

if ~isfield(params,'reconstruct')
    params.reconstruct = 'tv_tfocs';
end

[n1,n2]=size(inputImage);
if(n1~=n2)
    error('expect image with equal dimensions for each side')
end
n1log2=log2(n1);
if(n1log2~=floor(n1log2))
    error('expect image a size of 2 power')
end
    
% Set coarsest scale
delta = params.delta;
N=n1*n2;
M=floor(delta * N);

vec = @(x) x(:);
mat = @(x) reshape(x,n1,n2);

switch(params.reconstruct)
case 'tv_l1magic'
   make_coder=@(N,M,parms) tv_l1magic.codec(N,M,parms);   
case 'tv_tfocs'
   make_coder=@(N,M,parms) codec_tv_tfocs(N,M,parms);   
end
coder=make_coder(N,M,params);   
  
samples=CompressHere(inputImage,params);
[outputImage results]=ExtractHere(samples,params);
results.nSamples=length(samples);

function [samples] = CompressHere(inputImage,params)
    disp(params)
 
    samples= coder.encode(vec(inputImage));
    if length(samples)~=M
        warning('There was an error counting the samples')
    end

end


function [outputImage,results] = ExtractHere(samples,params)
    disp(params)
    results=struct();

    if length(samples)~=M
        warning('There was an error in the number of samples when reconstructing the image')
    end

    [x,results.fine] = coder.decode(samples);

    % Reconstruct
    outputImage = mat(x);
end

end
