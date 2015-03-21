% function [ImgRMinusBkgnd,ImgGMinusBkgnd] = SubBkgnd(imgR,imgG,params,varargin)
%
% Subtracts background from bead images for spotfinding.
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

function [imgRbkgnd,imgGbkgnd,imgRMinusBkgnd,imgGMinusBkgnd] = SubBkgnd(imgR,imgG,params,varargin)

BkgndSubSigma = 4;

% Medfilt2 works best with doubles:
if ~strcmpi(class(imgR),'double') || ~strcmpi(class(imgG),'double')
    disp('SubBkgnd: Warning: works best with doubles!')
end

% (1) Smooth with a boxcar average of width 2, with edge truncation
% Matlab doesn't have a boxcar filter built-in to the image processing
% toolbox, but it does have a Gaussian filter, with edge truncation as
% default:
h = fspecial('gaussian',params.DNASize,BkgndSubSigma);
imgRsmooth = imfilter(imgR,h);
imgGsmooth = imfilter(imgG,h);

% (2) IDL code uses either a median or min filt next (_brief uses med filt); 
% Matlab has a built-in median filter so using that:
imgRmedians = medfilt2(imgRsmooth,[params.DNASize*2, params.DNASize*2]);
imgGmedians = medfilt2(imgGsmooth,[params.DNASize*2, params.DNASize*2]);

% (3) Boxcar smooth again with window of 30 or 60 pixels (depending on
% whether you use _brief); again here going with a Gaussian filter:
h2 = fspecial('gaussian',params.DNASize.*5,BkgndSubSigma.*5);
imgRbkgnd = imfilter(imgRmedians,h2);
imgGbkgnd = imfilter(imgGmedians,h2);

imgRMinusBkgnd = imgR - imgRbkgnd;
imgGMinusBkgnd = imgG - imgGbkgnd;

if ~isempty(varargin)
    figure
    subplot(1,2,1)
    imshow(imgR,[])
    title('Red channel, raw')
    subplot(1,2,2)
    imshow(imgG,[])
    title('Green channel, raw')
    
    figure
    subplot(1,2,1)
    imshow(imgRsmooth,[])
    title('Red channel, smoothed')
    subplot(1,2,2)
    imshow(imgGsmooth,[])
    title('Green channel, smoothed')
    
    figure
    subplot(1,2,1)
    imshow(imgRmedians,[])
    title('Red channel, medians')
    subplot(1,2,2)
    imshow(imgGmedians,[])
    title('Green channel, medians')
    
    figure
    subplot(1,2,1)
    imshow(imgRbkgnd,[])
    title('Red channel, medians, smoothed')
    subplot(1,2,2)
    imshow(imgGbkgnd,[])
    title('Green channel, medians, smoothed')
    
    figure
    subplot(1,2,1)
    imshow(imgRMinusBkgnd,[])
    title('Red channel minus background')
    subplot(1,2,2)
    imshow(imgGMinusBkgnd,[])
    title('Green channel minus background')
    
    pause
    close all
end
