% function [imgR, imgG] = SplitImg(img,params_struct)
%
% Getting red and green channels out of a combined image. Called by smFRET,
% UserSpotSelection, and ScaleMovie.
%
% Inputs are the image or set of images to scale, and the output of
% smFRETsetup.m.
%
% Copyright (C) 2014 Stephanie Johnson, University of California, San Francisco
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% A copy of the GNU General Public License can be found in the LICENSE.txt 
% file that accompanies this software; it can also be found at 
% <http://www.gnu.org/licenses/>.

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