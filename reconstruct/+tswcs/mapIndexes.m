function idx2 =mapIndexes(map,idx1)

% similar to ind2=map(ind1), but allows zero indexes in ind1;

revmap=zeros(size(map));
fidx2=1:numel(map);
revmap(map)=fidx2;
idx2=zeros(size(map));
nonzero1=idx1~=0;
nonzero2=nonzero1(map);
idx2(nonzero2)=revmap(idx1(nonzero1));