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
% It can happen that a red channel spot, when transformed to the green
% channel, is too close to the edge to be useable.  Remove any such spots:
[imgtogetsize,~] = LoadScaledMovie(PathToMovie,[1 1]);
if length(find(Gspots(1,:)>=1+floor(params.DNASize/2)))~=length(Gspots(1,:))
    oldGspots = Gspots;
    clear Gspots
    Gspots = oldGspots(:,Gspots(1,:)>=1+floor(params.DNASize/2));
end
if length(find(Gspots(2,:)>=1+floor(params.DNASize/2)))~=length(Gspots(2,:))
    oldGspots = Gspots;
    clear Gspots
    Gspots = oldGspots(:,oldGspots(2,:)>=1+floor(params.DNASize/2));   
end
if length(find(Gspots(1,:)<=size(imgtogetsize,1)+floor(params.DNASize/2)))~=length(Gspots(1,:))
    oldGspots = Gspots;
    clear Gspots
    Gspots = oldGspots(:,Gspots(1,:)<=size(imgtogetsize,1)+floor(params.DNASize/2));
end
if length(find(Gspots(2,:)<=size(imgtogetsize,2)+floor(params.DNASize/2)))~=length(Gspots(2,:))
    oldGspots = Gspots;
    clear Gspots
    Gspots = oldGspots(:,Gspots(2,:)<=size(imgtogetsize,2)+floor(params.DNASize/2));    
end
clear imgtogetsize

for jj = 1:100:length(alltifs)
tic
    [imgR,imgG] = LoadScaledMovie(PathToMovie,[jj jj+99]);
%     temp = load(fullfile(PathToMovie,strcat('ScaledMovieFrames',int2str(jj),...
%         'to',int2str(jj+99),'.mat')));
%     imgR = temp.imgR;
%     imgG = temp.imgG;
%     clear temp
    
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
           bkgndR = zeros(1,size(spotimgR,3));
           bkgndG = zeros(1,size(spotimgR,3));
           for bb = 1:size(spotimgR,3)
               [~,~,~,~,bkgndR(bb),~] = Fit2DGaussToSpot(spotimgR(:,:,bb),'Background',...
                   'StartParams',[localcenR(1),localcenR(2),spotVars(1,kk),spotVars(2,kk),...
                   min(min(spotimgR(:,:,bb))),max(max(spotimgR(:,:,bb)))],'symGauss',params.UseSymGauss);
               [~,~,~,~,bkgndG(bb),~] = Fit2DGaussToSpot(spotimgG(:,:,bb),'Background',...
                   'StartParams',[localcenG(1),localcenG(2),spotVars(1,kk),spotVars(2,kk),...
                   min(min(spotimgG(:,:,bb))),max(max(spotimgG(:,:,bb)))],'symGauss',params.UseSymGauss);
           end
           % Subtract this background value from every pixel in the ROI
           % containing the spot:
           bkgndR = repmat(bkgndR,size(spotimgR,1)*size(spotimgR,2),1,1);
           bkgndR = reshape(bkgndR,size(spotimgR,1),size(spotimgR,2),size(spotimgR,3));
           bkgndG = repmat(bkgndG,size(spotimgG,1)*size(spotimgG,2),1,1);
           bkgndG = reshape(bkgndG,size(spotimgG,1),size(spotimgG,2),size(spotimgG,3));
            RedI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityV4(spotimgR-bkgndR,...
                localcenR,spotVars(:,kk));
            
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
   toc
end