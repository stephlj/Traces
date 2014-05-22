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
% if SpotBox(1,1)<=0 || SpotBox(1,2)>size(img,1) || ...
%     SpotBox(2,1)<=0 || SpotBox(2,2)>size(img,2)
%         disp(strcat('ExtractROI: Spot must be further than: ',int2str(boxdim/2),' pxls from edge of img'))
%         ROI=-1;
%         return
% end
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

% Not sure if there's a better way to write a universal formula?
localcen = zeros(size(spotcen));
if round(spotcen(1))<spotcen(1)
	localcen(1) = floor(boxdim/2)+1+abs(round(spotcen(1))-spotcen(1));
else
	localcen(1) = floor(boxdim/2)+1-abs(round(spotcen(1))-spotcen(1));
end
if round(spotcen(2))<spotcen(2)
	localcen(2) = floor(boxdim/2)+1+abs(round(spotcen(2))-spotcen(2));
else
	localcen(2) = floor(boxdim/2)+1-abs(round(spotcen(2))-spotcen(2));
end
% For debugging
% subplot(2,1,1)
% imshow(ROI)
% hold on
% plot(localcen(2),localcen(1),'xr')