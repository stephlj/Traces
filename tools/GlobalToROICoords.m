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