function obj=sampleQuantizer(params)

    if ~isfield(params,'quantizationType')
        params.quantizationType = 'none';
    end

    if ~isfield(params,'quantizationBits1')  % coarse sample
        params.quantizationBits1 = 8;
    end
    if ~isfield(params,'quantizationBits2')  % fine sample
        params.quantizationBits2 = 8;
    end
    bits1=params.quantizationBits1;
    bits2=params.quantizationBits2;

    Mcoarse=params.Mcoarse;  % selects M1 samples

    obj.type=params.quantizationType;
    obj.bits1=bits1;
    obj.bits2=bits2;
    coarse_sample_max=1e4;
    fine_sample_max=256;
    fine_sample_std=40;
    switch(params.quantizationType)
        case 'none'
            obj.encode= @(x) x;
            obj.decode= @(x) x;
        case 'linear'
            M1encode = @(x) encodeLinear(x,0,coarse_sample_max,bits1);% coarse 
            M1decode = @(x) decodeLinear(x,0,coarse_sample_max,bits1);
            M2encode = @(x) encodeLinear(x,-fine_sample_max,2*fine_sample_max,bits2); % fine
            M2decode = @(x) decodeLinear(x,-fine_sample_max,2*fine_sample_max,bits2);
            obj.encode=@encode;
            obj.decode=@decode;

        case 'cdf'

            M1encode = @(x) encodeLinear(x,0,coarse_sample_max,bits1);% coarse 
            M1decode = @(x) decodeLinear(x,0,coarse_sample_max,bits1);
            M2encode = @(x) encodeCdf(x,0,fine_sample_std,bits2); % fine
            M2decode = @(x) decodeCdf(x,0,fine_sample_std,bits2);
            obj.encode=@encode;
            obj.decode=@decode;
    end


    function y=encode(x)

        x1=x(1:Mcoarse);
        x2=x(Mcoarse+1:end);

        y1=M1encode(x1); 
        y2=M2encode(x2); 

        y=[y1;y2];
    end

    function x=decode(y)

        y1=y(1:Mcoarse);
        y2=y(Mcoarse+1:end);

        x1=M1decode(y1); 
        x2=M2decode(y2); 

        x=[x1;x2];

    end
end

function y=encodeCdf(x,u,sigma,bits)
    sample_max=2.^bits -1;
    y=normcdf(x,u,sigma);
    %sample_max+2
    % scale [1/257 , 256/257] to [0,255] (for 8 bits)
    y=y.*(sample_max+2)-1;  
    y=round(y);
    if any(y>sample_max) || any(y<0)
        y(y>sample_max)=sample_max;
        y(y<0)=0;
    end
end

function x=decodeCdf(y,u,sigma,bits)
    sample_max=2.^bits -1;
    y=(y+1)./(sample_max+2);  
    % scale [0,255] to [1/257 , 256/257] (for 8 bits)
    x=norminv(y,u,sigma);
    % F(x)=1, x=+inf
    % F(x)=0, x=-inf
    if any(isinf(x))
        error('samples not valid after icdt transform')
    end
end


function y=encodeLinear(x,offset,range,bits)
    sample_max=2.^bits -1;
    y=(x-offset)./range;  % scale 0 to 1
    y=y.*sample_max;  % scale 0 to sample_max
    y=round(y);
    if any(y>sample_max) || any(y<0)
        warning('truncating samples in quantization')
        y(y>sample_max)=sample_max;
        y(y<0)=0;
    end
end

function x=decodeLinear(y,offset,range,bits)
    sample_max=2.^bits -1;
    y=y./sample_max; % scale 0 to 1 
    x=(y.*range)+offset;

end
