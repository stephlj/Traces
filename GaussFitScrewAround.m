rootname =  'Nucleosomes4000xdil';

smFRETsetup;
params = load('AnalysisParameters.mat');

prevmapdir = load('PathToRecentMap.mat');
D_Beads = uigetdir(prevmapdir.MostRecentMapDir,'Select directory with old map');

Map = load(fullfile(D_Beads,'ChannelMapping.mat'));
A = Map.A;
b = Map.b;
Amatlab = Map.Amatlab;
bmatlab = Map.bmatlab;
clear Map prevmapdir

D_Data = uigetdir(D_Beads,'Select directory with data to analyze');
ToAnalyze = dir(fullfile(D_Data,strcat(rootname,'_*')));

i=1;

% Figure out the total number of image files in this movie:
alltifs = dir(fullfile(D_Data,ToAnalyze.name,'img*.tif'));
num_imgs = length(alltifs);

TotImg = LoadUManagerTifsV5(fullfile(D_Data,ToAnalyze(i).name),[1 params.FramesToAvg]);

if size(TotImg,3) > 1
    TotImg = mean(TotImg,3);
end

TotImg = mat2gray(TotImg);

[imgRed,imgGreen] = SplitImg(TotImg,params);

composite = CalcCombinedImage(A,b,imgGreen,imgRed);
composite = mat2gray(composite);

[spotsR,~,~] = FindSpotsV5(composite,'ShowResults',1,'ImgTitle','Composite Image',...
    'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize);
pause
close

%%
%Given how GetGaussParams works, I wonder if I should scale imgR and imgG
%separately?
[RefinedCenters,Vars] = GetGaussParams(spotsR,composite,imgGreen,imgRed,A,b, params.DNASize,1);

%%

figure('Position',[200,0,900,700])
subplot(2,3,1)
title('Xcen','Fontsize',14)
hold on
subplot(2,3,2)
title('Ycen','Fontsize',14)
hold on
subplot(2,3,3)
title('Amplitude','Fontsize',14)
hold on
subplot(2,3,4)
title('Xvar','Fontsize',14)
hold on
subplot(2,3,5)
title('Yvar','Fontsize',14)
hold on
subplot(2,3,6)
title('Background','Fontsize',14)
hold on

for j = 1:params.FramesToAvg:num_imgs
    j
    TotImg = LoadUManagerTifsV5(fullfile(D_Data,ToAnalyze(i).name),[j j-1+params.FramesToAvg]);

    if size(TotImg,3) > 1
        TotImg = mean(TotImg,3);
    end

    TotImg = mat2gray(TotImg);

    [imgRed,imgGreen] = SplitImg(TotImg,params);
    
    composite = CalcCombinedImage(Amatlab,bmatlab,imgGreen,imgRed);
    composite = mat2gray(composite);

    if j==1
        [spotsR,~,~] = FindSpotsV5(composite,'ShowResults',1,'ImgTitle','Composite Image',...
         'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize);
        pause
        close
        spotGaussParams = zeros(6,size(spotsR,2));
        spotGaussParamsR = zeros(6,size(spotsR,2));
        spotGaussParamsG = zeros(6,size(spotsR,2));
        Gspots = CalcSpotTransform([],spotsR,Amatlab,bmatlab);
        cumdiff = zeros(6,size(spotsR,2),num_imgs/params.FramesToAvg);
        cumdiffR = zeros(6,size(spotsR,2),num_imgs/params.FramesToAvg);
        cumdiffG = zeros(6,size(spotsR,2),num_imgs/params.FramesToAvg);
    end

    for ss = 1:size(spotsR,2)
       %Get ROI in composite image
       spotimg = ExtractROI(composite,params.DNASize,spotsR(:,ss));
       %Get ROI in red channel
       spotimgR = ExtractROI(imgRed,params.DNASize,spotsR(:,ss));
       %Get ROI in green channel:
       %Get coordinates of this spot in the other channel:
       %Gspots(:,ss) = CalcSpotTransform([],spotsR(:,ss),Amatlab,bmatlab);
       spotimgG = ExtractROI(imgGreen,params.DNASize,Gspots(:,ss));

       %Fit Gaussian in composite channel:
       %disp(strcat('Composite for spot number',int2str(ss)))
       [Xcen, Ycen, Xvar, Yvar, bkgnd, Amp] = Fit2DGaussToSpot(spotimg);%,'Debug',1);
       if j==1
           spotGaussParams(:,ss) = [Xcen-floor(params.DNASize/2)-1+spotsR(1,ss),...
               Ycen-floor(params.DNASize/2)-1+spotsR(2,ss),Xvar, Yvar, bkgnd, Amp];
       else
           cumdiff(:,ss,(j-1)/10)= transpose([Xcen-floor(params.DNASize/2)-1+spotsR(1,ss),...
               Ycen-floor(params.DNASize/2)-1+spotsR(2,ss),Xvar, Yvar,...
               bkgnd, Amp]) - spotGaussParams(:,ss);
       end
       clear Xcen Ycen Xvar Yvar bkgnd Amp

       %Fit Gaussian in red channel:
       %disp(strcat('Red for spot number',int2str(ss)))
       [Xcen, Ycen, Xvar, Yvar, bkgnd, Amp] = Fit2DGaussToSpot(spotimgR);%,'Debug',1,...
           %'StartParams',spotGaussParams(:,ss));
       if j==1
           spotGaussParamsR(:,ss) = [Xcen-floor(params.DNASize/2)-1+spotsR(1,ss),...
               Ycen-floor(params.DNASize/2)-1+spotsR(2,ss), Xvar, Yvar, bkgnd, Amp];
       else
           cumdiffR(:,ss,(j-1)/10)= transpose([Xcen-floor(params.DNASize/2)-1+spotsR(1,ss),...
               Ycen-floor(params.DNASize/2)-1+spotsR(2,ss),Xvar, Yvar,...
               bkgnd, Amp]) - spotGaussParamsR(:,ss);
       end
       clear Xcen Ycen Xvar Yvar bkgnd Amp

       %Fit Gaussian in green channel:
       %disp(strcat('Green for spot number',int2str(ss)))
       %newcen = CalcSpotTransform([],[spotGaussParams(1,ss),spotGaussParams(2,ss)],A,b);
       [Xcen, Ycen, Xvar, Yvar, bkgnd, Amp] = Fit2DGaussToSpot(spotimgG);%,'Debug',1,...
           %'StartParams',[newcen(1), newcen(2), spotGaussParams(3,ss),...
           %spotGaussParams(4,ss), spotGaussParams(5,ss), spotGaussParams(6,ss)]);
       if j==1
           spotGaussParamsG(:,ss) = [Xcen-floor(params.DNASize/2)-1+Gspots(1,ss),...
                Ycen-floor(params.DNASize/2)-1+Gspots(2,ss), Xvar, Yvar, bkgnd, Amp];
       else
           cumdiffG(:,ss,(j-1)/10)= transpose([Xcen-floor(params.DNASize/2)-1+Gspots(1,ss),...
                Ycen-floor(params.DNASize/2)-1+Gspots(2,ss), Xvar, Yvar,...
                bkgnd, Amp]) - spotGaussParamsG(:,ss);
       end
       clear Xcen Ycen Xvar Yvar bkgnd Amp
    end
    
end

subplot(2,3,1)
    hist(sum(cumdiff(1,:,:),3),100)
subplot(2,3,2)
    hist(sum(cumdiff(2,:,:),3),100)
subplot(2,3,3)
    hist(sum(cumdiff(6,:,:),3),100)
subplot(2,3,4)
    hist(sum(cumdiff(3,:,:),3),100)
subplot(2,3,5)
    hist(sum(cumdiff(4,:,:),3),100)
subplot(2,3,6)
    hist(sum(cumdiff(5,:,:),3),100)