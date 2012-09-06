% function WaveRelation2D_opWavelet_Test

sz=[16 16];
levels=2;
[IdxParent, IdxChildren, Ms]=WaveRelation2D_opWavelet(sz,levels);
sum(IdxParent==0)
parents=reshape(IdxParent,16,16);
children=reshape(IdxChildren,4,16,16);

sz=[128 128];
levels=4;
[IdxParent, IdxChildren, Ms]=WaveRelation2D_opWavelet(sz,levels);
