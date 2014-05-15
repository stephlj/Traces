% function ROI = ExtractROI(img,boxdim,spotcen)
%
% Given an image (img), extracts and returns a region of interest (ROI) of
% dimension boxdim on each side, centered at spotcen.  This is called by a
% number of functions to pick out the region that contains (hopefully) one
% spot, so that, for example, that spot's total intensity can be
% calculated.
%
% Also returns the center of the spot in the frame of reference of the ROI
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [ROI,localcen] = ExtractROI(img,boxdim,spotcen)

SpotCenRound = round(spotcen);

SpotBox = [SpotCenRound(1)-floor(boxdim/2), SpotCenRound(1)+floor(boxdim/2);...
    SpotCenRound(2)-floor(boxdim/2), SpotCenRound(2)+floor(boxdim/2)];
% Make sure the spot isn't too close to the edge of the image: 
if SpotBox(1,1)<=0 || SpotBox(1,2)>size(img,1) || ...
    SpotBox(2,1)<=0 || SpotBox(2,2)>size(img,2)
        disp(strcat('ExtractROI: Spot must be further than: ',int2str(boxdim/2),' pxls from edge of img'))
        ROI=-1;
        return
end
ROI = img(SpotBox(1,1):SpotBox(1,2),SpotBox(2,1):SpotBox(2,2),:);
localcen = spotcen-(round(spotcen)-[floor(boxdim)/2; floor(boxdim)/2]);