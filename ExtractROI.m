% function ROI = ExtractROI(img,boxdim,spotcen)
%
% Given an image (img), extracts and returns a region of interest (ROI) of
% dimension boxdim+1 on each side, centered at spotcen.  This is called by a
% number of functions to pick out the region that contains (hopefully) one
% spot, so that, for example, that spot's total intensity can be
% calculated.
%
% Also returns the center of the spot in the frame of reference of the ROI.
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [ROI,localcen] = ExtractROI(img,boxdim,spotcen)

% In keeping with Matlab's convention that a point n = (n_x,n_y) will plot
% onto the pixel ceil(n), I'm using ceil here, not round. NOTE that if you
% use round here, GlobalToRoiCoords will not give you the right conversion
% for global to local coordinates!
SpotCenRound = ceil(spotcen);

SpotBox = [SpotCenRound(1)-floor(boxdim/2), SpotCenRound(1)+floor(boxdim/2);...
    SpotCenRound(2)-floor(boxdim/2), SpotCenRound(2)+floor(boxdim/2)];

% Make sure the edges of the ROI aren't off the edges of the image:
numreductions = 0;
while SpotBox(1,1)<1 || SpotBox(1,2)>size(img,1) || ...
    SpotBox(2,1)<1 || SpotBox(2,2)>size(img,2)
        % disp('ExtractROI: Box edge out of bounds, reducing box size.')
        boxdim = boxdim-1;
        SpotBox = [SpotCenRound(1)-floor(boxdim/2), SpotCenRound(1)+floor(boxdim/2);...
            SpotCenRound(2)-floor(boxdim/2), SpotCenRound(2)+floor(boxdim/2)];
        numreductions = numreductions+1;
end
if numreductions > 0
    disp(strcat('ExtractROI: Box edge out of bounds, reduced boxdim to ',int2str(boxdim),'(reduced by ',int2str(numreductions),')'))
end
ROI = img(SpotBox(1,1):SpotBox(1,2),SpotBox(2,1):SpotBox(2,2),:);

localcen = GlobalToROICoords(spotcen,[],ceil(spotcen),boxdim,boxdim);

% For debugging
% subplot(2,1,1)
% imshow(ROI)
% hold on
% plot(localcen(2),localcen(1),'xr')