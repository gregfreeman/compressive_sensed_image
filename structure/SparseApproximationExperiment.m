function [outputImage,results]=SparseApproximationExperiment(inputImage,params)

params.imageSide=size(inputImage,1);
if ~isfield(params,'k_sparse')
    params.k_sparse=0.1;
end
if ~isfield(params,'qmf')
    params.qmf=MakeONFilter('Symmlet',8);
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
levels = n1log2;

vec = @(x) x(:);
mat = @(x) reshape(x,n1,n2);
W_op=opWavelet(n1,n2,'Custom',params.qmf,levels);


 
samples=CompressHere(inputImage,params);
[outputImage results]=ExtractHere(samples,params);
results.nSamples=length(samples);


function [samples] = CompressHere(inputImage,params)
    disp(params)
 
    alpha=W_op*vec(inputImage);   
    samples=zeros(size(alpha));
    % select M largest wavelet coefficeints
    [~, sIdx]=sort(abs(alpha),'descend');
    samples(sIdx(1:M))=alpha(sIdx(1:M));

end


function [outputImage,results] = ExtractHere(samples,params)

    % Reconstruct with k largest wavelet coefficeints

    outputImage = mat(W_op'*samples);
    results=[];
end

end
