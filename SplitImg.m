% function [imgR, imgG] = SplitImg(img,params_struct)
%
% Getting red and green channels out of a combined image. Called by smFRET,
% UserSpotSelection, and ScaleMovie.
%
% Inputs are the image or set of images to scale, and the output of
% smFRETsetup.m.
%
% Update 2/2014: Make use of the new PxlsToExclude parameter.
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

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