% function [RefinedCenters,Vars,bkgnd] = RefineCensByGauss(spots,img,ROIsize,Debug)
%
% Given a (2xnumspots) matrix of spots, and the image they come from
% ('img'), return a 2xnumspots set of spot centers refined by a Gaussian
% fit. Also return the x- and y-variances for each spot (also size
% 2xnumspots), and the background value from the Gaussian fit in bkgnd
% (1xnumspots size).
%
% The other inputs are: 
% ROIsize: dimension (side length, in pixels) of a square over which to fit
%   the Gaussian
% Debug: if this is 1, this will plot some figures of the fit, for
%   debugging
%
% Stephanie 4/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [RefinedCenters,Vars,bkgnd] = RefineCensByGauss(spots,img,ROIsize,Debug)

if ~exist('Debug','var')
    Debug=0;
elseif Debug~=1
    Debug=0;
end

% Setting up some parameters
cen_tolerance = 2; % To judge the goodness of Gaussian fit, don't let the 
    % refined center position differ from the maximum pixel by more than
    % this value. Realistically 1 is probably fine too, the refinement is
    % usually within 1 pixel.
defaultXvar = 0.3; %If the fit fails, set the spot variance to this value
    % For now taking this value from the Ha lab IDL code, but it's probably
    % system (camera, microscope, etc) specific
defaultYvar = defaultXvar;

RefinedCenters = zeros(2,size(spots,2));
Vars = zeros(2,size(spots,2));
bkgnd = zeros(1,size(spots,2));

for ss = 1:size(spots,2)
   % Extract an ROI over which to calculate the fit:
   spotimg = ExtractROI(img,ROIsize,spots(:,ss));
   
   % Fit Gaussian:
   if Debug
        disp(strcat('Composite for spot number',int2str(ss)))
        [Xcen, Ycen, Xvar, Yvar, bkgnd(ss), Amp] = Fit2DGaussToSpot(spotimg,'full','Debug',1);
        pause
        close all
   else
        [Xcen, Ycen, Xvar, Yvar, bkgnd(ss), Amp] = Fit2DGaussToSpot(spotimg);
   end
   % Switch to global coordinates (Xcen and Ycen are currently in the
   % coordinates of the spot)
   Xcen = Xcen-floor(ROIsize/2)-1+spots(1,ss);
   Ycen = Ycen-floor(ROIsize/2)-1+spots(2,ss);
   % Check that the fit was reasonable: Note that Fit2DGaussToSpot returns
   % Amp = 0.0001 if the fit errored out
   if Amp == 0.0001 || (FindSpotDists(spots(:,ss),[Xcen;Ycen])>cen_tolerance)
       % TODO: also check that amplitude goes to some small value of max at
       % a distance of ~5 pixels away from the max ... (like in the Ha code)
       RefinedCenters(:,ss) = spots(:,ss);
       Vars(:,ss) = [defaultXvar;defaultYvar];
   else
       RefinedCenters(:,ss) = [Xcen;Ycen];
       Vars(:,ss) = [Xvar,Yvar];
   end
   % TODO: if the fit fails, reject the spot?
   % TODO: histogram the Vars and bkgnd to see how much they vary
   clear Xcen Ycen Xvar Yvar Amp
end