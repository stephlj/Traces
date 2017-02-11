% function [RefinedCenters, Vars,varargout] = FindRefinedSpotCenters(imgs,spots,bkgnd_tolerance,params,debug)
%
% Given a frame (best to have it be an average of ~10 frames) and the locations of
% bright pixels identified as spots in that movie, fit a 2D Gaussian to 
% each spot to obtain a sub-pixel location of the center and variance of
% each spot. 
%
% The assumption is that this is called separately on spots in each
% channel.
%
% This also implements several further refinements of spots that used to be
% part of FindSpotsV5: (1) check that the new spot center isn't too far
% away from the old one; (2) check that each spot is surrounded by a
% sufficiently low background; (3) check that the maximum isn't just one
% pixel big.
%
% Inputs:
% imgs: x-by-y or x-by-y-by-frame total fields of view in which spots have
%   already been located
% spots: (x,y)-by-number locations of spots in imgs
% bkgnd_tolerance: mean intensity value around the spot has to go to less than the
    % fitted background value plus this tolerance. For beads, a good number is 0.01;
    % for DNAs with background subtracted, 0.01 is also good; for non-background subtracted,
    % 0.2 is better.
% params: file saved by smFRETsetup
% debug: If this is 1, plots the figures in Fit2DGaussToSpot. Default is 0.
%
% Outputs:
% RefinedCenters: (xcen,ycen)-by-spots sub-pixel resolution spot centers
% Vars: (xvar,yvar)-by-spots variances, which are the mean values of these
%   variances across however many frames are in imgs
% optional output: amplitudes from the fits, in case that information is
%   used to choose which set of variances (red channel vs green channel) to
%   keep in the function that calls this one.
%
% The MIT License (MIT)
% 
% Copyright (c) 2014 Stephanie Johnson, University of California, San Francisco
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [RefinedCenters, Vars,varargout] = FindRefinedSpotCenters(imgs,spots,bkgnd_tolerance,params,debug)

if ~exist('debug','var') debug = 0; end
if size(imgs,3)>1 imgs = mat2gray(mean(imgs,3)); end

cen_tolerance = 2; % To judge the goodness of Gaussian fit, don't let the 
    % refined center position differ from the maximum pixel by more than
    % this value. Realistically 1 is probably fine too, the refinement is
    % usually within 1 pixel.
VarTolerance = 1;

ROIsize = params.DNASize;
symGauss = params.UseSymGauss;

defaultXvar = 0.7; %If the fit fails, set the spot variance to this value
    % For now taking this value (0.3) from the Ha lab IDL code, but it's probably
    % system (camera, microscope, etc) specific.
    % Update 5/2014: Changed default to 0.7 which is more consistent with
    % our setup
defaultYvar = defaultXvar;

tempCenters = zeros(2,size(spots,2));
tempVars = zeros(2,size(spots,2)); % Vars will be:
    % First row will contain xvars, second row yvars
    % Columns will be spots
tempbkgnd = zeros(1,size(spots,2));
tempAmps = zeros(1,size(spots,2));

RefinedCenters = [];
Vars = [];
bkgnd = [];
Amps = [];

    for ss = 1:size(spots,2)
       % Extract an ROI over which to calculate the fit:
       [spotimg,localcen] = ExtractROI(imgs,ROIsize,spots(:,ss));
       % Note if the fit fails, Fit2DGaussToSpot will return the
       % StartParams, so make sure the StartParams for xvar and yvar are
       % the defaults hardcoded above
       [RefinedLocalCenX,RefinedLocalCenY,tempVars(1,ss),tempVars(2,ss),tempbkgnd(ss),tempAmps(ss)] = Fit2DGaussToSpot(spotimg,'Full',...
           'StartParams',[localcen(1),localcen(2),defaultXvar,defaultYvar,min(spotimg(:)),max(spotimg(:))],...
           'Debug',0,'symGauss',symGauss);
       tempCenters(:,ss) = GlobalToROICoords([],[RefinedLocalCenX;RefinedLocalCenY],ceil(spots(:,ss)),ROIsize,ROIsize);
       % for debugging
       % figure
%        imshow(imgs,[])
%        hold on
%        plot(spots(2,ss),spots(1,ss),'og')
%        plot(RefinedCenters(2,ss,i),RefinedCenters(1,ss,i),'xr')
%        pause 
%        close

        % Some further refinements:
        % Check that the maximum is surrounded by a sufficiently low
        % background:
        spotboundary = (sum(spotimg(:,1))+sum(spotimg(:,end))+...
            sum(spotimg(1,:))+sum(spotimg(end,:)))/(2*size(spotimg,1)+2*size(spotimg,2));
        % And check that the fit didn't fail and move the spot center somewhere
        % ridiculuous and that this isn't a one-pixel maximum by checking
        % the variance:
        if (FindSpotDists(tempCenters(:,ss),spots(:,ss))<cen_tolerance) && ... 
                spotboundary<=tempbkgnd(ss)+bkgnd_tolerance && ...
                tempVars(1,ss) < VarTolerance
            RefinedCenters(:,end+1) = tempCenters(:,ss);
            Vars(:,end+1) = tempVars(:,ss);
            bkgnd(end+1) = tempbkgnd(ss);
            Amps(end+1) = tempAmps(ss);
        else
            if debug
                if FindSpotDists(tempCenters(:,ss),spots(:,ss))>cen_tolerance
                    disp(sprintf('Spot center refinement failed: Bad fit convergence (new spot center is %d away)',FindSpotDists(tempCenters(:,ss),spots(:,ss))))
                end
                if spotboundary>bkgnd+bkgnd_tolerance
                    disp(sprintf('Spot center refinement failed: Too high background (boundary background is %d)',spotboundary))
                end
                if tempVars(1,ss) < VarTolerance
                    disp(sprintf('Spot center refinement failed: Variance is below tolerance (%d; tolerance is %d)',tempVars(1,ss),VarTolerance))
                end
            end
        end
        clear spotimg localcen tempVars tempCenters spotboundary empbkgnd
    end

varargout{1} = bkgnd;
varargout{2} = Amps;