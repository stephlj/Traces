% function I = CalcSpotIntensityV3(img,spotcen)
%
% Given a movie and a spot center, returns the total intensity of the spot 
% for each frame in the movie in vector I. Similar to CalcSpotIntensityV2,
% except fits a Gaussian to each spot to refine the total intensity estimation.
%
% My spot finding function (and all related functions that manipulate found
% spots) returns spots in the coordinate system of the matrix, that is, with
% the spot center given as (row;col), though row and col don't have to be integers.  
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [I,varargout] = CalcSpotIntensityV3(img,spotcen,SpotSize,varargin)

% The part of the image that contains the spot is going to be in:
SpotCenRound = round(spotcen);
SpotBox = [SpotCenRound(1)-floor(SpotSize/2), SpotCenRound(1)+floor(SpotSize/2);...
    SpotCenRound(2)-floor(SpotSize/2), SpotCenRound(2)+floor(SpotSize/2)];

% Make sure the spot isn't too close to the edge of the image: 
if SpotBox(1,1)<=0 || SpotBox(1,2)>size(img,1) || ...
        SpotBox(2,1)<=0 || SpotBox(2,2)>size(img,2)
    disp(strcat('Spot must be further than: ',int2str(SpotDiam/2),' pxls from edge of img'))
    I = -1;
    return
end

% If this spot has already been fit to a Gaussian, use those prior
% parameters as initial parameters for the next fit:
if ~isempty(varargin)
    init_params = varargin{1};
else
    init_params = [];
end

spotimg = img(SpotBox(1,1):SpotBox(1,2),SpotBox(2,1):SpotBox(2,2),:);

I = zeros(1,size(img,3),class(img));

for i = 1:size(img,3)
    [Xcen_out, Ycen_out, Xvar_out, Yvar_out, bkgnd_out, A_out] = Fit2DGaussToSpot(spotimg(:,:,i),'StartParams',...
        init_params);
    % Use a Gaussian with the same parameters except with amplitude 1 to 
    % create a weighted sum of the spot's intensities. Also subtract the
    % local background.
    spotimg_subbkgnd = spotimg(:,:,i) - bkgnd_out;
    I(i) = sum(sum(spotimg(:,:,i).*PlotGauss2D(size(spotimg(:,:,i)),...
        [Xcen_out, Ycen_out, Xvar_out, Yvar_out, 0, 1])));
    init_params = [Xcen_out, Ycen_out, Xvar_out, Yvar_out, bkgnd_out, A_out];
    
    if A_out>=0.4
        
        mask = zeros(5,5,class(img));
        mask(2:4,:) = ones(3,5,class(img));
        mask(:,2:4) = ones(5,3,class(img));
        
        figure('Position',[200,300,900,400])
        subplot(1,2,1)
        surf(spotimg(:,:,i))
        colormap gray
        title('Original image','Fontsize',14)
        zlim([0 1])

        subplot(1,2,2)
        surf(spotimg(:,:,i).*PlotGauss2D(size(spotimg(:,:,i)),...
        [Xcen_out, Ycen_out, Xvar_out, Yvar_out, 0, 1]))
        colormap jet
        title('Multiplied by Gaussian','Fontsize',14)
        zlim([0 1])
        
        totIraw = sum(sum(mask.*spotimg(3:end-2,3:end-2,i)));
        
        totIrefine = I(i);

        pause
        close
        clear mask
    end
    
end

varargout{1} = init_params;