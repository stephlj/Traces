% function [imgRBkgnd,imgGBkgnd] = CalcBkgnd(imgR,imgG,params)
%
% Calculates an image that CalcIntensitiesV3 uses to subtract a background.
%
% Inputs:
% imgR,imgG: the scaled intensity matrices saved by ScaleMovieV2.
% params: output of smFRETsetup, allows efficient passing of the parameters
%   used in this function (filter sizes, etc)
% 
% Outputs:
% imgRBkgnd,imgGBkgnd: background images. These are the background, they
%    are not images with background subtracted (in contrast to SubBkbnd)
%
% A note about algorithm used here: 
% The Ha lab creates a background image from only an average of the
% first 20 frames, I believe, or at least only the first frame. Here
% I'm calculating one image per frame for the entire movie. Otherwise, I'm
% mostly following their algorithm:
% (1) Smooth with a boxcar average of width 2, with edge truncation
% Matlab doesn't have a boxcar filter built-in to the image processing
% toolbox, but it does have a Gaussian filter, with edge truncation as
% default.
% (2) IDL code uses either a median or min filt next (_brief uses med filt); 
% Matlab has a built-in median filter (medfilt2), but it doesn't accept 3D
% images. You can do a minimum filter using ordfilt2, but again, no 3D.
% See work-around below. 
% (3) Boxcar smooth again with window of 30 or 60 pixels (depending on
% whether you used _brief).
%
% Steph 6/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [imgRBkgnd,imgGBkgnd] = CalcBkgnd(imgR,imgG,params)

if ~strcmpi(class(imgR),'double') || ~strcmpi(class(imgG),'double')
    % Medfilt2 works best with doubles:
    disp('CalcBkgnd: Warning: works best with doubles!')
end

% Setting up filter paramters:
% Some notes about choosing these parameters in fspecial:
% (1) The bigger you make the second input, which is the size of the
% neighborhood over which the Gaussian is applied, the slower imfilter will
% return the result.
% (2) As the third input, which is the width of the Gaussian to apply,
% approaches the second input, the result stops changing much.
% (3) As the second input gets much larger than the third, it also stops
% mattering much (except for being slower).
F = fspecial('gaussian',sqrt(params.DNANeighborhood),sqrt(params.DNANeighborhood)/4);
F2 = fspecial('gaussian',sqrt(params.DNANeighborhood)*2,params.DNASize);

% We want to remove background, so average away some of the noise first, if
% given a 3D image. Or do I want to smooth the resultant background at the
% end?
% if size(imgR,3)>1 && params.FramesToAvg>1
%     imgROrig = imgR;
%     imgGOrig = imgG;
%     clear imgG imgR
%     for j=params.FramesToAvg/2:params.FramesToAvg/2:(size(imgROrig,3)-params.FramesToAvg/2)
%         imgR(:,:,j-params.FramesToAvg/2+1:j+params.FramesToAvg/2) = repmat(median(imgROrig(:,:,j-params.FramesToAvg/2+1:j+params.FramesToAvg/2),3),...
%             1,1,length(j-params.FramesToAvg/2+1:j+params.FramesToAvg/2));
%         imgG(:,:,j-params.FramesToAvg/2+1:j+params.FramesToAvg/2) = repmat(median(imgGOrig(:,:,j-params.FramesToAvg/2+1:j+params.FramesToAvg/2),3),...
%             1,1,length(j-params.FramesToAvg/2+1:j+params.FramesToAvg/2));
%     end
% end
    
% Note that if img is 3-d, imfilter smoothes each frame independently,
% which at this point is what I want.
imgRSmooth = imfilter(imgR,F);
imgGSmooth = imfilter(imgG,F);

% % Unfortunately, medfilt2 must be 2D only ... can I avoid a for-loop here?
% imgRMedians = zeros(size(imgR));
% imgGMedians = zeros(size(imgG));
% for i=1:size(imgR,3)
%     imgRMedians(:,:,i) = medfilt2(imgSmooth(:,:,i),[sqrt(params.DNANeighborhood), sqrt(params.DNANeighborhood)]);
%     imgGMedians(:,:,i) = medfilt2(imgSmooth(:,:,i),[sqrt(params.DNANeighborhood), sqrt(params.DNANeighborhood)]);
% end
% Trying to avoid a for-loop with medfilt2:
% I should pad with zero's on the sides before reshaping if I'm going to do
% this:
% if params.splitx
%     imgRSmoothtemp = reshape(imgRSmooth,size(imgRSmooth,1),size(imgRSmooth,2)*size(imgRSmooth,3));
%     imgGSmoothtemp = reshape(imgGSmooth,size(imgGSmooth,1),size(imgGSmooth,2)*size(imgGSmooth,3));
% else
%     imgRSmoothtemp = reshape(imgRSmooth,size(imgRSmooth,1)*size(imgRSmooth,2),size(imgRSmooth,3));
%     imgGSmoothtemp = reshape(imgGSmooth,size(imgGSmooth,1)*size(imgGSmooth,2),size(imgGSmooth,3));
% end
% % This is a median filter:
% imgRMedianstemp = medfilt2(imgRSmoothtemp,[sqrt(params.DNANeighborhood), sqrt(params.DNANeighborhood)]);
% imgGMedianstemp = medfilt2(imgGSmoothtemp,[sqrt(params.DNANeighborhood), sqrt(params.DNANeighborhood)]);
% % whereas this would be a min filter:
% % imgRMedianstemp = ordfilt2(imgRSmoothtemp,1,ones(sqrt(params.DNANeighborhood),sqrt(params.DNANeighborhood)));
% % imgGMedianstemp = ordfilt2(imgGSmoothtemp,1,ones(sqrt(params.DNANeighborhood),sqrt(params.DNANeighborhood)));
% imgRMedians = reshape(imgRMedianstemp,size(imgRSmooth,1),size(imgRSmooth,2),size(imgRSmooth,3));
% imgGMedians = reshape(imgGMedianstemp,size(imgGSmooth,1),size(imgGSmooth,2),size(imgGSmooth,3));

% Trying imopen here, because it gives roughly the same 
% result, but is much faster than medfilt2. 

imgRMedians = imopen(imgRSmooth,strel('disk',params.DNASize/2));
imgGMedians = imopen(imgGSmooth,strel('disk',params.DNASize/2));

imgRBkgnd = imfilter(imgRMedians,F2);
imgGBkgnd = imfilter(imgGMedians,F2);

% I would kind of like to do some kind of time-averaging here, but this
% takes forever ... there's probably a smarter and faster way to do it?
% if size(imgR,3)>1 && params.FramesToAvg>1
%     for i=1:size(imgR,1)
%         for j=1:size(imgR,2)
%             imgRBkgnd(i,j,:) = smooth(imgRBkgnd(i,j,:),params.FramesToAvg);
%             imgGBkgnd(i,j,:) = smooth(imgGBkgnd(i,j,:),params.FramesToAvg);
%         end
%     end
% end
