%% Load first 20 frames of a movie
D = '/Volumes/smFRET/smFRET data/14Mar21/Nucs75pMSnf2h100nMATP1mMAfter_3';
allimgs = double(LoadUManagerTifsV5(D,[1 20]));
imshow(mean(allimgs,3),[])
%% If neccesary histogram all intensity values to figure out a nice scaling
temp = reshape(allimgs,1,size(allimgs,1)*size(allimgs,2)*size(allimgs,3));
hist(temp,1000)
%%
allimgs2 = mat2gray(allimgs,[0 800]);
imshow(mean(allimgs2,3))
%%
hold on
plot([10 10+3.74*10],[10 10],'-w','Linewidth',3') %Scale bar is 10 um long