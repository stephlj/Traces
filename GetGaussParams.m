% function [RefinedCenters,Vars] = GetGaussParams(spotsR,composite,imgG,imgR,A,b,ROIsize,Debug)
%
% Given a set of spot centers (spots) and the image they're found in
% (composite, the single-channel image produced by CalcCombinedImage), fit 2D
% Gaussians to each spot and return a refined center position for the spot
% in the coordinate of composite image (which for me is the acceptor
% channel), plus the X and Y variances of each spot.
%
% Inputs:
% spotsR: (x,y)-by-numspots matrix of spots in the coordinate system of the
%   composite image
% composite: a single-channel, composite image. This will have the
%   coordinate system of one of the channels; for me it's the acceptor
%   channel.
% A, b: matrix that performs the transformation from the other channel into
%   the coordinate system of the composite image's channel
% imgR, imgG: the two channels separately. This function compares the
%   results of fitting a 2D Gaussian to the composite image versus the
%   channels individually to best refine the spot's position and variance.
% ROIsize: how big an area around the spot center to use for Gaussian
%   fitting
% Debug: enter 1 to show some figures to check these fits and the decision
%   about which spots to keep
%
% Outputs: 
% RefinedCenters: new x,y coordinates for all spots where a "good" Gaussian
%   could be fit, returned in the coordinate system of the acceptor (red)
%   channel (i.e. the composite image)
% Vars: x and y variances for each spot
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [RefinedCenters,Vars] = GetGaussParams(spotsR,composite,imgG,imgR,A,b, ROIsize,varargin)

if isempty(varargin)
    Debug=0;
elseif varargin{1}==1
    Debug=1;
else
    Debug=0;
end

RefinedCenters = zeros(2,size(spotsR,2));
Vars = zeros(2,size(spotsR,2));
% Find the spots in the coordinate system of the other (green) channel:
Gspots = CalcSpotTransform([],spotsR,A,b);

for ss = 1:size(spotsR,2)
   % Get ROI in composite image
   spotimg = ExtractROI(composite,ROIsize,spotsR(:,ss));
   % Get ROI in red channel
   spotimgR = ExtractROI(imgR,ROIsize,spotsR(:,ss));
   % Get ROI in green channel:
   % Get coordinates of this spot in the other channel:
   spotimgG = ExtractROI(imgG,ROIsize,Gspots(:,ss));

   % Fit Gaussian in composite channel:
   if Debug
        disp(strcat('Composite for spot number',int2str(ss)))
        [Xcen, Ycen, Xvar, Yvar, bkgnd, Amp] = Fit2DGaussToSpot(spotimg,'Debug',1);
   else
       [Xcen, Ycen, Xvar, Yvar, bkgnd, Amp] = Fit2DGaussToSpot(spotimg);
   end
   Xcen = Xcen-floor(ROIsize/2)-1+spotsR(1,ss);
   Ycen = Ycen-floor(ROIsize/2)-1+spotsR(2,ss);

   % Fit Gaussian in red channel:
   if Debug
        disp(strcat('Red for spot number',int2str(ss)))
        [XcenR, YcenR, XvarR, YvarR, bkgndR, AmpR] = Fit2DGaussToSpot(spotimgR,'Debug',1);
   else 
       [XcenR, YcenR, XvarR, YvarR, bkgndR, AmpR] = Fit2DGaussToSpot(spotimgR);
   end
   XcenR = XcenR-floor(ROIsize/2)-1+spotsR(1,ss);
   YcenR = YcenR-floor(ROIsize/2)-1+spotsR(2,ss);

   % Fit Gaussian in green channel:
   if Debug
        disp(strcat('Green for spot number',int2str(ss)))
        [XcenG, YcenG, XvarG, YvarG, bkgndG, AmpG] = Fit2DGaussToSpot(spotimgG,'Debug',1);
   else 
       [XcenG, YcenG, XvarG, YvarG, bkgndG, AmpG] = Fit2DGaussToSpot(spotimgG);
   end
   XcenG = XcenG-floor(ROIsize/2)-1+Gspots(1,ss);
   YcenG = YcenG-floor(ROIsize/2)-1+Gspots(2,ss);

   % Pick which best represent the true spot center:
   % As far as I've seen, there are four possibilities here:
   % (1) The spot in the composite image corresponds to a spot mostly found
   % in the green channel
   % (2) or mostly in the red channel
   % (3) The spot has significant intensities in both channels (mid-FRET)
   % (4) Or this is a junk spot--e.g., a "spot" was found in the composite
   % image that's actually one spot in the red channel and one in green
   % that have nearby centers but are not the same
   % To distinguish among these possibilities, first compare AmpG and AmpR
   % to Amp. If neither AmpG nor AmpR corresponds to most of Amp, then
   % compare the centers.
   % TODO: Further refinement, like checking background
   
   if AmpG >= 0.6*Amp && AmpR < 0.6*Amp
       % Trust the green channel fit best
       % Check that picking out an ROI based on the composite image didn't
       % mess up the fit (if the center picked out by spot-finding on the
       % composite image was way off from the "true" center):
       if FindSpotDists([XcenG,YcenG],[Xcen,Ycen]) >= 3
           spotimgG = ExtractROI(imgG,ROIsize,[XcenG,YcenG]);
           oldcen = [XcenG,YcenG];
           if Debug
                [XcenG, YcenG, XvarG, YvarG, bkgndG, AmpG] = Fit2DGaussToSpot(spotimgG,'Debug',1);
           else
               [XcenG, YcenG, XvarG, YvarG, bkgndG, AmpG] = Fit2DGaussToSpot(spotimgG);
           end
           XcenG = XcenG-floor(ROIsize/2)-1+oldcen(1);
           YcenG = YcenG-floor(ROIsize/2)-1+oldcen(2);
           clear oldcen
       end
       % Convert back to the coordinate system of the other channel:
       RefinedCenters(:,ss) = CalcSpotTransform([XcenG,YcenG],[],A,b);
       Vars(:,ss) = [XvarG,YvarG];
       if Debug
            disp('Green channel wins.')
       end
   elseif AmpR >= 0.6*Amp && AmpG < 0.6*Amp
       % Trust the red channel fit best
       % Again check that the fit used the best ROI:
       if FindSpotDists([XcenR,YcenR],[Xcen,Ycen]) >= 3
           spotimgR = ExtractROI(imgR,ROIsize,[XcenR,YcenR]);
           oldcen = [XcenR,YcenR];
           if Debug
                [XcenR, YcenR, XvarR, YvarR, bkgndR, AmpR] = Fit2DGaussToSpot(spotimgR,'Debug',1);
           else 
               [XcenR, YcenR, XvarR, YvarR, bkgndR, AmpR] = Fit2DGaussToSpot(spotimgR);
           end
           XcenR = XcenR-floor(ROIsize/2)-1+oldcen(1);
           YcenR = YcenR-floor(ROIsize/2)-1+oldcen(2);
           clear oldcen
       end
       RefinedCenters(:,ss) = [XcenR,YcenR];
       Vars(:,ss) = [XvarR,YvarR];
       if Debug
            disp('Red channel wins.')
       end
   elseif FindSpotDists([XcenG,YcenG],[Xcen,Ycen]) <= 2 && ...
           FindSpotDists([XcenR,YcenR],[Xcen,Ycen]) <= 2
       % This is the same spot in both channels, mid-FRET probably. Go with
       % whichever fit has higher amplitude, as probably more accurate
       if AmpG > AmpR
            % Trust the green channel fit best
            % Convert back to the coordinate system of the other channel:
            RefinedCenters(:,ss) = CalcSpotTransform([XcenG,YcenG],[],A,b);
            Vars(:,ss) = [XvarG,YvarG];
            if Debug
                disp('Green channel wins.')
            end
       else
           % Trust the red channel fit best
           RefinedCenters(:,ss) = [XcenR,YcenR];
           Vars(:,ss) = [XvarR,YvarR];
           if Debug
                disp('Red channel wins.')
           end
       end
   else
       RefinedCenters(:,ss) = [0,0];
       Vars(:,ss) = [0,0];
       if Debug
            disp('Bad spot, nobody wins.')
       end
   end
   
   if Debug
        pause
        close all
   end
   
   clear spotimg spotimgG spotimgR
   clear Xcen Ycen XcenG YcenG XcenR YcenR Xvar Yvar XvarG YvarG
   clear XvarR YvarR bkgnd bkgndR bkgndG Amp AmpR AmpG
   
end

RefinedCenters = RefinedCenters(:,find(RefinedCenters(1,:)));
Vars = Vars(:,find(Vars(1,:)));
disp(sprintf('Kept %d of %d spots', size(RefinedCenters,2), size(spotsR,2)))
