% function composite = CalcCombinedImage(tform,StartImg,EndImg,ShowResults)
%
% Given two images and the transformation encoded in the matlab-specific
% structure tform that maps from StartImg to EndImg, 
% calculate the transformed version of StartImg and return a
% composite of EndImg and the transformed image.
%
% This is called by smFRET so that you can find spots in a composite image
% from both channels, the idea being that medium-FRET spots will have
% decreased intensities in the two channels and may be missed by a
% spotfinding algorithm that only looks for spots in green and spots in red
% separately.
%
% NOTE this currently uses the built-in Matlab function imwarp, which was
% released with R2012a. Older Matlab versions cannot run this function.
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

function composite = CalcCombinedImage(tform,StartImg,EndImg,ShowResults)

if ~exist('ShowResults','var') ShowResults = 0; end

if isnan(tform)
    % Traces passes in a NaN for tform if Kind = Poly in FRETmap.
    disp('CalcCombinedImage: Your version of Matlab does not support this function.')
    composite = -1;
    return
end

try 
    Rfixed = imref2d(size(EndImg));
    alignedimgG = imwarp(StartImg,tform,'OutputView',Rfixed);
catch
    %alignedimgG = imtransform(StartImg,tform); % With R2011b and older, you
        %can use imtransform instead of imwarp, but there's no imfuse
        %equivalent (or imshowpair).
    disp('CalcCombinedImage: Your version of Matlab does not support this function.')
    composite = -1;
    return
end

composite = imfuse(EndImg,alignedimgG,'blend');

if ShowResults
    figure, imshowpair(EndImg,alignedimgG,'blend')
    title('Green overlayed on red')
end