% function [imgR, imgG] = SplitImg(img,params_struct)
%
% Getting red and green channels out of a combined image. Called by smFRET,
% UserSpotSelection, and ScaleMovie.
%
% Inputs are the image or set of images to scale, and the output of
% smFRETsetup.m.
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

function [imgR,imgG] = SplitImg(img,params_struct)
    if params_struct.splitx
        if ~params_struct.Acceptor
            imgG = img(:,(size(img,2)/2)+1+params_struct.PxlsToExclude:end-params_struct.PxlsToExclude,:); 
            imgR = img(:,1+params_struct.PxlsToExclude:(size(img,2)/2-params_struct.PxlsToExclude),:);
        else
            imgR = img(:,(size(img,2)/2)+1+params_struct.PxlsToExclude:end-params_struct.PxlsToExclude,:); 
            imgG = img(:,1+params_struct.PxlsToExclude:(size(img,2)/2)-params_struct.PxlsToExclude,:);
        end
    else
        if ~params_struct.Acceptor
            imgG = img((size(img,2)/2)+1+params_struct.PxlsToExclude:end-params_struct.PxlsToExclude,:,:); 
            imgR = img(1+params_struct.PxlsToExclude:(size(img,2)/2)-params_struct.PxlsToExclude,:,:);
        else
            imgR = img((size(img,2)/2)+1+params_struct.PxlsToExclude:end-params_struct.PxlsToExclude,:,:); 
            imgG = img(1+params_struct.PxlsToExclude:(size(img,2)/2)-params_struct.PxlsToExclude,:,:);
        end
    end
end