% function I = CalcSpotIntensityV2(img,spotcen)
%
% This version uses a circular area to define a spot, rather than a
% square.  The size of the circle is hard-coded in right now.  You cannot
% pass this function a spot that is right at the edge of the image!--needs
% to be more than SpotDiam/2 = 2 pxls from the edges of the image.
%
% Given a movie and a spot center, returns the sum of the intensities in a 
% circle around the spot center for each frame in the movie.
%
% My spot finding function (and all related functions that manipulate found
% spots) returns spots in the coordinate system of the matrix, that is, with
% the spot center given as (row;col), though row and col don't have to be integers.  
%
% Steph 10/2013
% Copyright 2013 Stephanie Johnson, University of California, San Francisco

function I = CalcSpotIntensityNoGauss(img,spotcen)

% The way I'm going to define a circle, or really a disk, that contains the
% spot is by creating a square mask where some elements are zero, and
% multiplying that by the image.  DNAs in our hands are about 3x3 pixels, so
% doing a 5 pxl diameter circle around each spot. (4 is hard.)

SpotDiam = 5;

mask = zeros(SpotDiam,SpotDiam,class(img));
mask(2:4,:) = ones(3,SpotDiam,class(img));
mask(:,2:4) = ones(SpotDiam,3,class(img));

% The part of the image that contains the spot is going to be in:
SpotCenRound = round(spotcen);
SpotBox = [SpotCenRound(1)-floor(SpotDiam/2), SpotCenRound(1)+floor(SpotDiam/2);...
    SpotCenRound(2)-floor(SpotDiam/2), SpotCenRound(2)+floor(SpotDiam/2)];

% Make sure the spot isn't too close to the edge of the image: 
if SpotBox(1,1)<=0 || SpotBox(1,2)>size(img,1) || ...
        SpotBox(2,1)<=0 || SpotBox(2,2)>size(img,2)
    disp(strcat('Spot must be further than: ',int2str(SpotDiam/2),' pxls from edge of img'))
    I = -1;
    return
end

spotimg = img(SpotBox(1,1):SpotBox(1,2),SpotBox(2,1):SpotBox(2,2),:);

I = zeros(1,size(img,3),class(img));

% %For debugging:
%     figure
%     imshow(mean(img(newbox(1,1):newbox(2,1),newbox(1,2):newbox(2,2),1:20),3),[])
%     pause
%     close

for i = 1:size(img,3) %Can I avoid a for loop?
    I(i) = sum(sum(mask.*spotimg(:,:,i)));
end