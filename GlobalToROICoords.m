% function newcoords = GlobalToROICoords(spot_global_coords,spot_local_coords,ROI_center_global,boxSizeX,boxSizeY)
%
% Since I switch back and forth between the location of a spot in the
% coordinates of the whole (256x512) image, and its location in an ROI,
% this function encapsulates that calculation. If you have the coordinates
% in the ROI (i.e. local coordinates), pass spot_global_coords = [] and 
% spot_local_coords = [xcen;ycen]. If you have the global
% coordinates, pass spot_global_coords = [xcen;ycen] and spot_local_coords = [].
%
% ROI_center_global is the center of the ROI in GLOBAL coordinates.
% BoxSizeX and boxSizeY are the width and height of the ROI respectively,
% where width refers to the X-dimension and height to the y-dimension.
% (The terminology I am using here is that the first element of spot_global_coords,
% spot_local_coords, and ROI_center_global is the "x" dimension, and the second
% element is the "y" dimension; note that this is not consistent with the way 
% Matlab plots coordinates on top of images, but it is consistent with how I've 
% been using the terms "x" and "y" in the rest of my code).
%
% Steph 6/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function newcoords = GlobalToROICoords(spot_global_coords,spot_local_coords,ROI_center_global,boxSizeX,boxSizeY)

if isempty(spot_global_coords)
    % User has local coords, wants global
    newcoords = [spot_local_coords(1)+ceil(ROI_center_global(1))-floor(boxSizeX)/2-1;...
        spot_local_coords(2)+ceil(ROI_center_global(2))-floor(boxSizeY)/2-1];
elseif isempty(spot_local_coords)
    % User has global, wants local
    newcoords = [spot_global_coords(1)-(ceil(ROI_center_global(1))-floor(boxSizeX)/2-1);...
        spot_global_coords(2)-(ceil(ROI_center_global(2))-floor(boxSizeY)/2-1)];
else
    disp('GlobalToROICoords: Either the first or second input must be an empty matrix.')
    newcoords = [];
    return
end