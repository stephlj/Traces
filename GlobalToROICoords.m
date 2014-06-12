% function newcoords = GlobalToROICoords(globalcoords,ROIcoords,center,boxSizeX,boxSizeY)
%
% Since I switch back and forth between the location of a spot in the
% coordinates of the whole (256x512) image, and its location in an ROI,
% this function encapsulates that calculation. If you have the coordinates
% in the ROI, pass globalcoords = [] and ROIcoords=[xcen;ycen]; if you have the global
% coordinates, pass globalcoords = [xcen;ycen] and ROIcoords = [].
%
% Center is the center of the ROI in GLOBAL coordinates.
% BoxSizeX and boxSizeY are the width and height of the ROI respectively,
% where width refers to the X-dimension and height to the y-dimension.
% It is assumed that globalcoords, ROIcoords, and center are passed as
% (x;y).
%
% Steph 6/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function newcoords = GlobalToROICoords(globalcoords,ROIcoords,center,boxSizeX,boxSizeY)

if isempty(globalcoords)
    % User has local coords, wants global
    newcoords = [round(center(1))-floor(boxSizeX/2)-1+ROIcoords(1);...
        round(spots(2))-floor(boxSizeY/2)-1+ROIcoords(2)];
elseif isempty(ROIcoords)
    % User has global, wants local
    % Not sure if there's a better way to write a universal formula?
    newcoords = zeros(size(center));
    if round(globalcoords(1))<globalcoords(1)
        newcoords(1) = floor(boxSizeX/2)+1+abs(round(center(1))-globalcoords(1));
    else
        newcoords(1) = floor(boxSizeX/2)+1-abs(round(center(1))-globalcoords(1));
    end
    if round(globalcoords(2))<globalcoords(2)
        newcoords(2) = floor(boxSizeY/2)+1+abs(round(center(2))-globalcoords(2));
    else
        newcoords(2) = floor(boxSizeY/2)+1-abs(round(center(2))-globalcoords(2));
    end
else
    disp('GlobalToROICoords: Either the first or second input must be an empty matrix.')
    newcoords = [];
    return
end