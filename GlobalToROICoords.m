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
% Copyright (C) 2014 Stephanie Johnson, University of California, San Francisco
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% A copy of the GNU General Public License can be found in the LICENSE.txt 
% file that accompanies this software; it can also be found at 
% <http://www.gnu.org/licenses/>.

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