%function composite = CalcCombinedImage(A,b,StartImg,EndImg)
%
%Given two images and the transformation matrix and vectors A, b that 
%define an affine map from StartImg to EndImg (see outputs of CalcChannelMapping), 
%calculate the transformed version of StartImg and return a
%composite of EndImg and the transformed image.
%
%This is called by smFRET so that you can find spots in a composite image
%from both channels, the idea being that medium-FRET spots will have
%decreased intensities in the two channels and may be missed by a
%spotfinding algorithm that only looks for spots in green and spots in red
%separately.
%
%Note: you could also modify this code to find tform using built-in Matlab
%functions:
%   Use fitgeotrans with pairs of matching points:
%   tform = fitgeotrans(matchRall',matchGall','Affine');
%   where matchRall, matchRall are, for example, the pairs of points found
%   in smFRET from bead images. Or, use the Matlab function
%   cpselect to generate points manually.
%
%Update 1/2014: fitgeotrans does a bit better than my handwritten code in
%CalcChannelMapping; it might be iterating to find optimal solutions or
%something like that.  Regardless, I now use fitgeotrans, but, for
%backwards compatibility, turn the output of fitgeotrans into A, b, and
%then re-form tform here.  
%
%Copyright 2014 Stephanie Johnson, University of California, San Francisco

function composite = CalcCombinedImage(A,b,StartImg,EndImg)

tformMatrix = [A; -b'];
tformMatrix = [tformMatrix,[0; 0; 1]];

tform = affine2d(tformMatrix);

Rfixed = imref2d(size(EndImg));

alignedimgG = imwarp(StartImg,tform,'OutputView',Rfixed);

composite = imfuse(EndImg,alignedimgG,'blend');

figure, imshowpair(EndImg,alignedimgG,'blend')