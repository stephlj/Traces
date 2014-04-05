% function I = CalcSpotIntensityV4(img,spotcen,spotvar)
%
% Calculates the total intensity of a fluorescent spot for every frame in a
% movie. The intensity at every pixel is weighted by a Gaussian centered at
% spotcen, with x-variance spotvar(1) and y-variance spotvar(2).
%
% Note this assumes background has been subtracted already, that is, the
% Gaussian that is used does not have a background offset. Updated 4/2014
% so it doesn't make this assumption but allows the user to decide whether
% the background is 0 or not.
%
% Inputs:
% img: x-by-y-by-frames for a REGION OF INTEREST of a movie. The Gaussian
%   will be fit over the entire size of img(:,:,1) so pass it only a small
%   box containing the spot!
% spotcen: (x,y) of the center of the spot. 
% spotvar: (xvar,yvar)
% bkgnd: the Gaussian will be shifted in the z direction by this amount.
%   Defaults to 0 if no bkgnd is passed.
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function I = CalcSpotIntensityV4(img,spotcen,spotvar,bkgnd)

if ~exist('bkgnd','var')
    bkgnd = 0;
end

I = zeros(1,size(img,3),class(img));

spotsize = size(img(:,:,1));

for i=1:size(img,3)
    I(i) = sum(sum(img(:,:,i).*PlotGauss2D(spotsize,...
        [spotcen(1), spotcen(2), spotvar(1), spotvar(2), bkgnd, 2])));
    % Note: The Ha lab IDL code weights by a Gaussian of Amp 2 instead of
    % 1...
end