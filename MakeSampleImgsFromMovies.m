%% Load first 20 frames of a movie
D = '/Users/Steph/Dropbox/Steph Dropbox/Narlikar Lab DB/smFRET data/13Sept13/BeadsMovie_1';
allimgs = LoadUManagerTifsV5(D,[1 20]);
imshow(mean(allimgs,3),[])
%% If neccesary histogram all intensity values to figure out a nice scaling
temp = reshape(allimgs,1,size(allimgs,1)*size(allimgs,2)*size(allimgs,3));
hist(temp,1000)
%%
allimgs2 = mat2gray(allimgs,[0 500]);
imshow(mean(allimgs2,3))
%%
hold on
plot([10 10+3.74*10],[10 10],'-w','Linewidth',3')