% function [ImgRMinusBkgnd,ImgGMinusBkgnd] = SubBkgndFullImg(img,params,debug,varargin)
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
% Steph 6/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [imgBkgnd,imgMinusBkgnd] = SubBkgndFullImg(img,params,debug,varargin)

if ~exist('debug','var') debug = 0; end

% Medfilt2 works best with doubles:
if ~strcmpi(class(img),'double')
    disp('SubBkgndFullImg: Warning: works best with doubles!')
end

% (1) Smooth with a boxcar average of width 2, with edge truncation
% Matlab doesn't have a boxcar filter built-in to the image processing
% toolbox, but it does have a Gaussian filter, with edge truncation as
% default.
% Some notes about choosing these parameters:
% (1) The bigger you make the second input, which is the size of the
% neighborhood over which the Gaussian is applied, the slower imfilter will
% return the result.
% (2) As the third input, which is the width of the Gaussian to apply,
% approaches the second input, the result stops changing much.
% (3) As the second input gets much larger than the third, it also stops
% mattering much (except for being slower).
F = fspecial('gaussian',sqrt(params.DNANeighborhood),params.DNASize);
% Note that if img is 3-d, imfilter smoothes each frame independently,
% which at this point is what I want.
imgSmooth = imfilter(img,F);

% (2) IDL code uses either a median or min filt next (_brief uses med filt); 
% Matlab has a built-in median filter so using that:
% Unfortunately, medfilt2 must be 2D only ... can I avoid a for-loop here?
imgMedians = zeros(size(img));
for i=1:size(img,3)
    imgMedians(:,:,i) = medfilt2(imgSmooth(:,:,i),[sqrt(params.DNANeighborhood), sqrt(params.DNANeighborhood)]);
end

% (3) Boxcar smooth again with window of 30 or 60 pixels (depending on
% whether you use _brief); again here going with a Gaussian filter:
F2 = fspecial('gaussian',sqrt(params.DNANeighborhood)*2,params.DNASize*2);
imgBkgnd = imfilter(imgMedians,F2);

% Allow the option of smoothing in time, not just in space:
if ~isempty(varargin) && size(img,3)>1
    for j=varargin{1}/2:1:size(img,3)-varargin{1}/2
        imgBkgnd(:,:,j-varargin{1}/2+1:j+varargin{1}/2) = median(imgBkgnd(:,:,j-varargin{1}/2+1:j+varargin{1}/2),3);
    end
end

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
