% function [RedI,GrI] = CalcIntensitiesV2(PathToMovie, Rspots, spotVars,...
%       tform,params)
%
% Calculates intensities for all spots in a movie.
%
% Inputs:
% PathToMovie: a path to a directory with the movie
% Rspots: locations of the spots in the acceptor channel
% spotVars: x,y variances to use for each spot for the Gaussian weighting
%   performed by CalcSpotIntensityV4
% tform: mapping information for finding spots in the donor channel
% params: file saved by smFRETsetup
%
% Outputs:
% RedI, GrI: Intensity-vs-time information for each spot
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [RedI, GrI] = CalcIntensitiesV2(PathToMovie, Rspots, spotVars, ...
    tform,params)

% Figure out the total number of image files in this movie:
alltifs = dir(fullfile(PathToMovie,'img*.tif'));

RedI = zeros(size(Rspots,2),length(alltifs));
GrI = zeros(size(Rspots,2),length(alltifs));
% Find the spots in the coordinate system of the other (green) channel:
Gspots = transpose(transformPointsInverse(tform,Rspots'));
    
% Updated 2/2014: scaling the entire movie such that the brightest pixel is
% 1 and the dimmest is 0.  This is equivalent to mat2gray(), but I can't
% load all frames into memory at once.  So, ScaleMovie() loads in 100
% frames at a time, calculates the min and max for the whole movie, then
% below each 100-frame chunk is again loaded and scaled to the min and max
% of the whole movie.

[overallMin,overallMax] = ScaleMovie(PathToMovie,length(alltifs));

for jj = 1:100:length(alltifs)
    moviebit = LoadUManagerTifsV5(PathToMovie,[jj jj+99]);
    % As noted in ScaleMovie, it's not clear to me if I should be scaling
    % each channel separately, or as one image.  Right now ScaleMovie
    % returns only one min and one max so scaling the image globally:
    moviebit = mat2gray(moviebit,[overallMin overallMax]); % This also converts it to double precision
    [imgRraw,imgGraw] = SplitImg(moviebit,params);
    
    % Subtract time-local background:
    % For debugging:
    %[imgRbkgnd,imgGbkgnd,imgR,imgG] = SubBkgnd(mean(imgRraw,3),mean(imgGraw,3),params,1);
%     [imgRbkgnd,imgGbkgnd,~,~] = SubBkgnd(mean(imgRraw,3),mean(imgGraw,3),params);
%     imgR = imgRraw-repmat(imgRbkgnd,1,1,size(imgRraw,3));
%     imgG = imgGraw-repmat(imgGbkgnd,1,1,size(imgRraw,3));
    imgR = imgRraw;
    imgG = imgGraw;
    
    % For debugging:
%     temp = reshape(imgRraw,1,size(imgRraw,1)*size(imgRraw,2)*size(imgRraw,3));
%     hist(temp,1000)
%     xlim([0 1])
%     title('Red, Raw')
%     clear temp
%     temp = reshape(imgRbkgnd,1,size(imgRbkgnd,1)*size(imgRbkgnd,2));
%     figure,hist(temp,1000)
%     xlim([0 1])
%     title('Red, Bkgnd')
%     clear temp
%     temp = reshape(imgR,1,size(imgR,1)*size(imgR,2)*size(imgR,3));
%     figure,hist(temp,1000)
%     xlim([0 1])
%     title('Red, Minus bkgnd')
%     clear temp
%     pause 
%     close all

    % Save the background-subtracted, scaled image for later use in
    % UserSpotSelection:
    save(fullfile(PathToMovie,strcat('ScaledMovieFrames',int2str(jj),...
        'to',int2str(jj+99),'.mat')),'imgR','imgG')
    
    for kk = 1:size(Rspots,2)
        % Get ROI in red channel
       spotimgR = ExtractROI(imgR,params.DNASize,Rspots(:,kk));
       localcenR = Rspots(:,kk)-(round(Rspots(:,kk))-[floor(params.DNASize)/2; floor(params.DNASize)/2]);
       % Get ROI in green channel:
       spotimgG = ExtractROI(imgG,params.DNASize,Gspots(:,kk));
       localcenG = Gspots(:,kk)-(round(Gspots(:,kk))-[floor(params.DNASize)/2; floor(params.DNASize)/2]);
       if params.IntensityGaussWeight==1
            RedI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityV4(spotimgR,localcenR,spotVars(:,kk));
            GrI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityV4(spotimgG,localcenG,spotVars(:,kk));
       else
           RedI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityNoGauss(imgR,Rspots(:,kk));
           GrI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityNoGauss(imgG,Gspots(:,kk));
       end
       clear spotimgG spotimgR
    end
   clear imgR imgG moviebit
   disp(sprintf('Calculated intensity for frames %d to %d of %d', jj, jj+99,length(alltifs)))
end