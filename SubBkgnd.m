% function [ImgRMinusBkgnd,ImgGMinusBkgnd] = SubBkgnd(imgR,imgG,params,varargin)
%
% Steph 3/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [imgRMinusBkgnd,imgGMinusBkgnd] = SubBkgnd(imgR,imgG,params,varargin)

% Medfilt2 works best with doubles:
if ~strcmpi(class(imgR),'double') || ~strcmpi(class(imgG),'double')
    disp('SubBkgnd: Warning: works best with doubles!')
end

% (1) Smooth with a boxcar average of width 2, with edge truncation
% Matlab doesn't have a boxcar filter built-in to the image processing
% toolbox, but it does have a Gaussian filter, with edge truncation as
% default:
h = fspecial('gaussian',params.DNASize,params.BkgndSubSigma);
imgRsmooth = imfilter(imgR,h);
imgGsmooth = imfilter(imgG,h);

% (2) IDL code uses either a median or min filt next (_brief uses med filt); 
% Matlab has a built-in median filter so using that:
imgRmedians = medfilt2(imgRsmooth,[params.DNASize, params.DNASize]);
imgGmedians = medfilt2(imgGsmooth,[params.DNASize, params.DNASize]);

% (3) Boxcar smooth again with window of 30 or 60 pixels (depending on
% whether you use _brief); again here going with a Gaussian filter:
h2 = fspecial('gaussian',params.DNASize.*5,params.BkgndSubSigma.*5);
imgRbkgnd = imfilter(imgRmedians,h2);
imgGbkgnd = imfilter(imgGmedians,h2);

imgRMinusBkgnd = imgR - imgRbkgnd;
imgGMinusBkgnd = imgG - imgGbkgnd;

if ~isempty(varargin{1})
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
