% function CSUnstructuredWaveletExperimentTest


data=FoveatedImageData(7);%,'resize',[256 256]);   
img=data.image;
clear data;


params=struct();
params.delta=0.3;
params.qmf = MakeONFilter('Symmlet',8);
params.reconstruct = 'spgl1';
[outputImage,results]=CSUnstructuredWaveletExperiment(img,params);

subplot(1,2,1)
imagesc(outputImage)
colormap gray
subplot(1,2,2)
imagesc(img)
colormap gray
