function mapIndexesTest


% similar to ind2=map(ind1), but allows zero indexes in ind1;

idx1=[0 0 0 0 2 8 0 10 0 0];
map=[2 4 5 6 8 10];
idx2=mapIndexes(map,idx1);
idx2