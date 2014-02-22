%function center = FindSpotCentroid(spot,varargin)
%
%Takes as its input an image of a single spot and returns the [x,y] coordinates
%of the centroid of the spot.
%
%Optional input is a threshold value to apply to the spots before
%calculating the centroid.
%
%Steph 6/2013

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


