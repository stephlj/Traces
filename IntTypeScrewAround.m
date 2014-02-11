D = '/Users/Steph/Dropbox/Steph Dropbox/Narlikar Lab DB/smFRET data/13Sept13/Nucleosomes4000xdil_3';
moviebit = LoadUManagerTifsV5(D,[1 20]);
moviebitReshaped = reshape(moviebit,512*20*512,1);
hist(double(moviebitReshaped),100)
set(gca,'YScale','log')
ylabel('Counts','Fontsize',14)
xlabel('Raw intensity (a.u.)','Fontsize',14)
set(gca,'Fontsize',14)

rescaled = (2^16-1).*mat2gray(moviebit);
rescaledReshaped = reshape(rescaled,512*20*512,1);
figure, hist(double(rescaledReshaped),100)
set(gca,'YScale','log')
ylabel('Counts','Fontsize',14)
xlabel('Scaled intensity (a.u.)','Fontsize',14)
set(gca,'Fontsize',14)

figure, imshow(mean(moviebit,3),[])
figure, imshow(mean(rescaled,3),[0 (2^16-1)])
%The above are different!  Why?  because of the order of mat2gray (which is
%implicit in imshow(mean(moviebit,3),[]))
figure, imshow(mean(mat2gray(moviebit),3)) %Note though that this is no longer scaled between 0 and 1!
figure, imshow(mat2gray(mean(moviebit,3)))
%Note further that mat2gray(mean(moviebit,3)) = mat2gray(mean(mat2gray(moviebit),3))
%So if you want the mean scaled between 0 and 1, use mean first then
%mat2gray

rescaled2 = mat2gray(moviebit);
rescaled2Reshaped = reshape(rescaled2,512*20*512,1);
figure, hist(double(rescaled2Reshaped),100)
set(gca,'YScale','log')
ylabel('Counts','Fontsize',14)
xlabel('Scaled intensity (a.u.)','Fontsize',14)
set(gca,'Fontsize',14)
%So there's no loss of precision if I use mat2gray alone vs. times 2^16-1

%So, what I want to do is:
%To scale an entire movie, first run a for-loop to load in each 100 frame
%increment, and calculate the min and the max over the whole movie.
%Then run the same for-loop but this time scale each chunk of 100 frames
%according to
%moviebit = mat2gray(moviebit,[overallMin overallMax]);
%Note that this will also convert everything to doubles.
%Note to run the first for-loop, initialize min to be the maximum possible
%value (2^16-1), and max to be the minimum possible value (0).  Then at
%each step, overallMin = min(overallmin,min(moviebit)); etc.