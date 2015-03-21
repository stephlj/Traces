% function center = FindSpotCentroid(spot,varargin)
%
% Takes as its input an image of a single spot and returns the [x,y] coordinates
% of the centroid of the spot.
%
% Optional input is a threshold value to apply to the spots before
% calculating the centroid.
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

function center = FindSpotCentroid(spot,varargin)

center = [0,0];

if ~isempty(varargin)
    spot = mat2gray(spot,double([varargin{1}*max(max(spot)) max(max(spot))]));
end

[sy,sx] = size(spot);
tot = sum(spot(:)); %Total intensity

%The centroid or center of mass is a weighted average of the pixels and
%their intensities
center(1) = sum(sum(spot,1).*(1:sx))/tot;%The pixels go from 1:sx; if we were doing
    %the centroid like region props, then we would do ... ?.  The inner sum
    %weights by intensities.
    %By the way, sum(A,1) sums all the columns, which is x in an image
center(2) = sum(sum(spot,2)'.*(1:sy))/tot;

%For debugging
% imgwboxes(:,:,1) = spot;
% imgwboxes(:,:,2) = spot;
% imgwboxes(:,:,3) = spot;
% 
% imgwboxes(round(center(2)),round(center(1)),1) = 1;
% imgwboxes(round(center(2)),round(center(1)),2) = 0;
% imgwboxes(round(center(2)),round(center(1)),3) = 0;
% 
% figure
% imshow(imgwboxes)
% pause
% close
% clear imgwboxes


