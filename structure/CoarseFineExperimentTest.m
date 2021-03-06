% function CoarseFineExperimentTest


data=FoveatedImageData(7);%,'resize',[256 256]);   
img=data.image;
clear data;


params=struct();
params.delta=0.3;
params.qmf = MakeONFilter('Symmlet',8);
params.reconstruct = 'lasso_fista';
[outputImage,results]=CoarseFineExperiment(img,params);

subplot(1,2,1)
imagesc(outputImage)
colormap gray
subplot(1,2,2)
imagesc(img)
colormap gray
