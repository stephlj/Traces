% function [ImgRMinusBkgnd,ImgGMinusBkgnd] = SubBkgndFullImg(img,params,debug)
%
% Subtracting background from smFRET images. The diffrences between this
% and SubBkgnd are: 
% (1) SubBkgnd accepts two images at once, ideally a red image and a green 
% one, and is also slightly out of date (though it works well with the beads, 
% I haven't really tested the parameter values and settings with DNAs.) 
% (2) SubBkgnd cannot deal with 3d inputs, only 2D images.
%
% Inputs:
% img: 2D or 3D matrix of intensity values, assumed to be of the form
%   x-by-y-by frames (if 3D).
% params: output of smFRETsetup, allows efficient passing of the parameters
%   used in this function (filter sizes, etc)
% debug: pass 1 to plot some figures at the end. Default is 0.
% optional input: If you pass a 3D image and want to average the background
%   across some number of frames, pass that number of frames as the last
%   input
% 
% Outputs:
% imgBkgnd: the background that was subtracted
% imgMinusBkgnd: image minus the background
%
% A note about algorithm used here: 
% The Ha lab does background subtraction in 3 steps:
% (1) Smooth with a boxcar average of width 2, with edge truncation
% Matlab doesn't have a boxcar filter built-in to the image processing
% toolbox, but it does have a Gaussian filter, with edge truncation as
% default. Imfilter, however, can be quite slow.
% (2) IDL code uses either a median or min filt next (_brief uses med filt); 
% Matlab has a built-in median filter (medfilt2), but it doesn't accept 3D
% images. You can do a minimum filter using ordfilt2, but again, no 3D.
% (3) Boxcar smooth again with window of 30 or 60 pixels (depending on
% whether you use _brief).
%
% I have found that using Matlab's imopen works as well as or better than 
% the above, and requires fewer parameters to adjust.  Though I've added 
% a smoothing filter before and after calling imopen, for best results.
%
% Steph 6/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [imgBkgnd,imgMinusBkgnd] = SubBkgndFullImg(img,params,debug)

if ~exist('debug','var') debug = 0; end

if params.PxlsToExclude<sqrt(params.DNANeighborhood)/2
    disp('SubBkgndFullImg: Will work better if PxlsToExclude is less than sqrt(params.DNANeighborhood)/2.')
end

% Medfilt2 works best with doubles:
if ~strcmpi(class(img),'double')
    disp('SubBkgndFullImg: Warning: works best with doubles!')
end

% We want to remove background, so average away some of the noise first, if given a 3D image:
if size(img,3)>1 && params.FramesToAvg>1
    imgOrig = img;
    clear img
    for j=params.FramesToAvg/2:1:(size(imgOrig,3)-params.FramesToAvg/2)
        img(:,:,j-params.FramesToAvg/2+1:j+params.FramesToAvg/2) = repmat(median(imgOrig(:,:,j-params.FramesToAvg/2+1:j+params.FramesToAvg/2),3),...
            1,1,length(j-params.FramesToAvg/2+1:j+params.FramesToAvg/2));
        clear meds
    end
end

% (1) Smooth a little before calling imopen:
% Some notes about choosing these parameters in fspecial:
% (1) The bigger you make the second input, which is the size of the
% neighborhood over which the Gaussian is applied, the slower imfilter will
% return the result.
% (2) As the third input, which is the width of the Gaussian to apply,
% approaches the second input, the result stops changing much.
% (3) As the second input gets much larger than the third, it also stops
% mattering much (except for being slower).
%F = fspecial('gaussian',sqrt(params.DNANeighborhood),params.DNASize/2);
% Note that if img is 3-d, imfilter smoothes each frame independently,
% which at this point is what I want.
%imgSmooth = imfilter(img,F);
imgSmooth = img;

% Using imopen here, because it gives roughly the same result, but is much
% faster than medfilt2. I'm probably not using it as it's supposed to, but
% again it does do pretty well ... 
% % Unfortunately, medfilt2 must be 2D only ... can I avoid a for-loop here?
% imgMedians = zeros(size(img));
% for i=1:size(img,3)
%     imgMedians(:,:,i) = medfilt2(imgSmooth(:,:,i),[sqrt(params.DNANeighborhood), sqrt(params.DNANeighborhood)]);
% end
%imgMedians = imopen(imgSmooth,strel('disk',params.DNASize/2));
% Trying to avoid a for-loop with medfilt2:
if params.splitx
    imgSmoothtemp = reshape(imgSmooth,size(imgSmooth,1),size(imgSmooth,2)*size(imgSmooth,3));
else
    imgSmoothtemp = reshape(imgSmooth,size(imgSmooth,1)*size(imgSmooth,2),size(imgSmooth,3));
end
% imgMedianstemp = medfilt2(imgSmoothtemp,[sqrt(params.DNANeighborhood), sqrt(params.DNANeighborhood)]);
imgMedianstemp = ordfilt2(imgSmoothtemp,1,ones(sqrt(params.DNANeighborhood),sqrt(params.DNANeighborhood)));
if params.splitx
    imgMedians = reshape(imgMedianstemp,size(imgSmooth,1),size(imgSmooth,2),size(imgSmooth,3));
else
    imgMedians = reshape(imgMedianstemp,size(imgSmooth,1),size(imgSmooth,2),size(imgSmooth,3));
end

F2 = fspecial('gaussian',sqrt(params.DNANeighborhood)*2,params.DNASize);
imgBkgnd = imfilter(imgMedians,F2);

imgMinusBkgnd = img - imgBkgnd;

if debug
    
    imtool(img(:,:,1))
    title('Raw')
    
    imtool(imgSmooth(:,:,1))
    title('Smoothed')
    
    imtool(imgMedians(:,:,1))
    title('Smoothed and median filtered')
    
    imtool(imgBkgnd(:,:,1))
    title('Final background to be subtracted')
    
    imtool(imgMinusBkgnd(:,:,1))
    title('Img minus background')
    
    pause
    close all
end
