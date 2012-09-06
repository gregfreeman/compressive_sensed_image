function [IdxParent, IdxChildren,Ms]=WaveRelation2D_opWavelet(sz,levels)
% WaveRelation2D_opWavelet: parent/children relationships for coefficients
% of 2-dimensional wavelet transform
% USAGE: [IdxParent, IdxChildren]=WaveRelation2D(sz,levels)
% INPUT:  sz image size [n1,n2] for n1xn2 image 
%         levels - number of decomposition levels
% OUTPUT: IdxParent: 1 x M, parent index. IdxParent(i)=j means the parent of the ith coefficient is coefficient j 
%               M is the number of elements in the image (n1xn2)
%         IdxChildren: nChildren x M, children index. IdxChildren(i,:)=[j(1),...,j(nChildren)] means
%                      the children of the ith coefficient are coefficients
%                      j(1),...,j(nChildren). 
%              coeffients that don't have parents or children have 0 for
%              the index

%-------------------------------
% initial code by Lihan He, ECE, Duke University
% rewrote by Greg Freeman, University of texas - for use with the opWavelet
% wavelet operator.
%-------------------------------


    n1=sz(1);
    n2=sz(2);
    N=n1*n2;
    logn=log2(n1);
%     coarse_length=2^(logn-levels);
%     mat = @(x) reshape(x,n1,n2);
%     vec = @(x) x(:);
%     fine_map=true(n1,n2);
%     fine_map(1:coarse_length,1:coarse_length)=false;
%     
    parent_map=zeros(4,n1,n2);
    child_map=zeros(n1,n2);
    % ---------
    % | a | v |
    % |-------|
    % | h | d |
    % ---------
    %         a   h   v   d
    aOffset=[1 0;0 1;1 1];
    cOffset=[0 0;1 0;0 1;1 1];
    
   % coarse level coefficients have 3 children
   k=(logn-levels);
%    for x1=1:2^k
%    for x2=1:2^k
%        psub=[x1 x2];
%        for a=1:3
%             csub=psub+aOffset(a,:)*2^k;
%             child_map(csub(1), csub(2))=sub2ind(sz,psub(1),psub(2));
%             parent_map(psub(1),psub(2),a)=sub2ind(sz,csub(1), csub(2));
%        end
%    end
%    end
   Ms=[3*4^k];
   % other level coefficients have 4 children   
    for k=(logn-levels):(logn-2)
        Ms=[Ms 3*4^(k+1)];
        for x1=1:2^k
        for x2=1:2^k
            for a=1:3
                psub=[x1 x2]+aOffset(a,:)*2^k;
                for c=1:4
                    csub=psub*2+cOffset(c,:)+[-1 -1];
                    child_map(csub(1), csub(2))=sub2ind(sz,psub(1),psub(2));
                    parent_map(c,psub(1),psub(2))=sub2ind(sz,csub(1), csub(2));
                end            
            end
        end
        end
    end
    
    IdxParent=reshape(child_map,1,N);
    IdxChildren=reshape(parent_map,4,N);
    