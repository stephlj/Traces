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
% Steph 6/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [RedI, GrI] = CalcIntensitiesV3(PathToMovie, Rspots, spotVars, ...
    tform,params)

% Figure out the total number of image files in this movie:
alltifs = dir(fullfile(PathToMovie,'img*.tif'));

RedI = zeros(size(Rspots,2),length(alltifs));
GrI = zeros(size(Rspots,2),length(alltifs));
% Find the spots in the coordinate system of the other (green) channel:
if ~isempty(tform)
    Gspots = tform.FRETmapInv(Rspots);
    % HACK because I'm working with a cutout right now:
    % Gspots = zeros(size(Rspots));
    % RedSpotsGlobalCoords = zeros(size(Rspots));
    % for hh=1:size(Rspots,2)
    %     RedSpotsGlobalCoords(:,hh) = GlobalToROICoords([],Rspots(:,hh),[136;187],88,88);
    % end
    % GspotsGlobal = tform.FRETmapInv(RedSpotsGlobalCoords);
    % for hh=1:size(Rspots,2)
    %     Gspots(:,hh) = GlobalToROICoords(GspotsGlobal(:,hh),[],[136;187],88,88);
    % end
else
    Gspots = Rspots;
end

% It can happen that a red channel spot, when transformed to the green
% channel, is too close to the edge to be useable.  Remove any such spots:
[Gspots,Rspots,spotVars] = CheckSpotBoundaries(Gspots,...
    Rspots,spotVars,[],params,PathToMovie);

    % Subfunction to do the actual intensity calculation:
    function I = CalcSpotIntensityInternal(kind,img,spotcen,spotvar,params,bkgnd)

        if ~exist('bkgnd','var')
            bkgnd = 0;
        end
        if ~isempty(params.FixSpotVar)
            spotvar = params.FixSpotVar;
        end

        if strcmpi(kind,'Gauss')
            spotsize = size(img(:,:,1));
            mask = PlotGauss2D(spotsize,...
               [spotcen(1), spotcen(2), spotvar(1), spotvar(2), bkgnd, params.GaussWeightAmp]);
        else
            SpotDiam = 5;
            mask = zeros(SpotDiam,SpotDiam,class(img));
            mask(2:4,:) = ones(3,SpotDiam,class(img));
            mask(:,2:4) = ones(SpotDiam,3,class(img));
        end

        % Instead of using a for-loop here, reshape so that we can take
        % advantage of Matlab's fast matrix multiplication:
        % (though this would probably matter more if I could load more than
        % 100 movie frames at once ... )
        newmask = reshape(mask,1,size(mask,1)*size(mask,2));
        data = reshape(img,size(img,1)*size(img,2),size(img,3));
        %clear mask
        %mask = repmat(newmask,size(img,3),1);
        I = newmask*data;
        
    end

for jj = 1:100:length(alltifs)
    [imgR,imgG,bkgndR,bkgndG] = LoadScaledMovie(PathToMovie,[jj jj+99]);
    
    for kk = 1:size(Rspots,2)
        
       if params.IntensityGaussWeight==1
            % Get ROI in red channel
           [spotimgR,localcenR] = ExtractROI(imgR,params.DNASize,Rspots(:,kk));
           [spotRbkgnd,~] = ExtractROI(bkgndR,params.DNASize,Rspots(:,kk));
           % Get ROI in green channel:
           [spotimgG,localcenG] = ExtractROI(imgG,params.DNASize,Gspots(:,kk));
           [spotGbkgnd,~] = ExtractROI(bkgndG,params.DNASize,Rspots(:,kk));
           
           RedI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityInternal('Gauss',...
               spotimgR-spotRbkgnd,localcenR,spotVars(:,kk),params);
           GrI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityInternal('Gauss',...
               spotimgG-spotGbkgnd,localcenG,spotVars(:,kk),params);

       else
           % Get ROI in red channel
           [spotimgR,localcenR] = ExtractROI(imgR,5,Rspots(:,kk));
           % Get ROI in green channel:
           [spotimgG,localcenG] = ExtractROI(imgG,5,Gspots(:,kk));

           RedI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityInternal('NoGauss',...
               spotimgR,localcenR,[],params);
           GrI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensityInternal('NoGauss',...
               spotimgG,localcenG,[],params);
       end
       clear spotimgG spotimgR
    end
    
   clear imgR imgG bkgndR bkgndG
   disp(sprintf('Calculated intensity for frames %d to %d of %d', jj, jj+99,length(alltifs)))
end
end
