% function [RedI,GrI] = CalcIntensities(PathToMovie, Rspots, A, b)
%
% Calculates intensities for all spots in a movie.
%
% Inputs:
% PathToMovie: a path to a directory with the movie
% numframes: the number of frames in the movie
% Rspots: locations of the spots in the acceptor channel
% A, b: mapping information for finding spots in the donor channel
% params: file saved by smFRETsetup
%
% Outputs:
% RedI, GrI: Intensity-vs-time information for each spot
% imgRinit, imgGinit: Averaged and scaled images of the first few frames, for use
%   in display in UserSpotSelection, which calls this function
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [allRedI, allGrI, Gspots, imgRinit, imgGinit] = CalcIntensities(PathToMovie, Rspots, A, b,params)

% Figure out the total number of image files in this movie:
alltifs = dir(fullfile(PathToMovie,'img*.tif'));
allRedI = zeros(size(Rspots,2),length(alltifs));
allGrI = zeros(size(Rspots,2),length(alltifs));
Gspots = zeros(size(Rspots));
    
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
    moviebit = mat2gray(moviebit,[overallMin overallMax]); %This also converts it to double precision
    [imgR,imgG] = SplitImg(moviebit,params);
    
    if jj==1
        if size(imgR,3)>=params.FramesToAvg
            imgRinit = mat2gray(mean(imgR(:,:,1:params.FramesToAvg),3));
            imgGinit = mat2gray(mean(imgG(:,:,1:params.FramesToAvg),3));
        else
            imgRinit = mat2gray(mean(imgR,3));
            imgGinit = mat2gray(mean(imgG,3));
        end
    end
    
    for kk = 1:size(Rspots,2)
       allRedI(kk,jj:jj+99) = CalcSpotIntensityV2(imgR,Rspots(:,kk));
       % Find matching green spot:
       Gspots(:,kk) = inv(A)*(Rspots(:,kk)-repmat(b,1,size(Rspots(:,kk),2)));
       allGrI(kk,jj:jj+99) = CalcSpotIntensityV2(imgG,Gspots(:,kk));
    end
   clear imgR imgG moviebit
   disp(sprintf('Calculated intensity for frames %d to %d of %d', jj, jj+99,length(alltifs)))
end