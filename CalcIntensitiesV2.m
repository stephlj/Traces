% function [RedI,GrI] = CalcIntensitiesV2(PathToMovie, Rspots, spotVars,...
%       tform,params,varargin)
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
% Steph 2/2014, updated 5/2014 to calculate a background value for each
% frame
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [RedI, GrI] = CalcIntensitiesV2(PathToMovie, Rspots, spotVars, ...
    tform,params)

% Figure out the total number of image files in this movie:
alltifs = dir(fullfile(PathToMovie,'img*.tif'));

RedI = zeros(size(Rspots,2),length(alltifs));
GrI = zeros(size(Rspots,2),length(alltifs));
% Find the spots in the coordinate system of the other (green) channel:
Gspots = tform.FRETmapInv(Rspots);

for jj = 1:100:length(alltifs)

    temp = load(fullfile(PathToMovie,strcat('ScaledMovieFrames',int2str(jj),...
        'to',int2str(jj+99),'.mat')));
    imgR = temp.imgR;
    imgG = temp.imgG;
    clear temp
    
    for kk = 1:size(Rspots,2)
        % Get ROI in red channel
       [spotimgR,localcenR] = ExtractROI(imgR,params.DNASize,Rspots(:,kk));
       % localcenR = Rspots(:,kk)-(round(Rspots(:,kk))-[floor(params.DNASize)/2; floor(params.DNASize)/2]);
       % Get ROI in green channel:
       [spotimgG,localcenG] = ExtractROI(imgG,params.DNASize,Gspots(:,kk));
       % localcenG = Gspots(:,kk)-(round(Gspots(:,kk))-[floor(params.DNASize)/2; floor(params.DNASize)/2]);
       if params.IntensityGaussWeight==1
           % Update 5/2014: First fit a Gaussian, with the only parameters
           % that can vary the background and the amplitude, to get the
           % local background in this frame.
           % TODO: Should I take an average over some number of frames?
           [~,~,~,~,bkgndR,~] = Fit2DGaussToSpot(spotimgR,'Background',...
               'StartParams',[localcenR(1,kk),localcenR(2,kk),spotVars(1,kk),spotVars(2,kk),...
               min(spotimgR(:)),max(spotimgR(:))],'symGauss',params.UseSymGauss);
           % Subtract this background value from every pixel in the ROI
           % containing the spot:
            RedI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityV4(spotimgR-bkgndR,...
                localcenR,spotVars(:,kk));
            [~,~,~,~,bkgndG,~] = Fit2DGaussToSpot(spotimgG,'Background',...
               'StartParams',[localcenG(1,kk),localcenG(2,kk),spotVars(1,kk),spotVars(2,kk),...
               min(spotimgG(:)),max(spotimgG(:))],'symGauss',params.UseSymGauss);
            GrI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityV4(spotimgG-bkgndG,...
                localcenG,spotVars(:,kk));
            clear bkgndR bkgndG
       else
           RedI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityNoGauss(imgR,Rspots(:,kk));
           GrI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityNoGauss(imgG,Gspots(:,kk));
       end
       clear spotimgG spotimgR
    end
    
   clear imgR imgG 
   disp(sprintf('Calculated intensity for frames %d to %d of %d', jj, jj+99,length(alltifs)))
end