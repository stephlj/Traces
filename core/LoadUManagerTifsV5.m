% function allimgs = LoadUManagerTifsV5(D,varargin)
%
% Loads images in a uManager-created folder into one 3-d matrix.
% D is the full path to the folder where they're stored; output is the
% image matrix.  In V5, images returned are NOT scaled between 0 and 1.
% 
% Optional inputs: Pass these as '<name>',<val> pairs
% 'FramesToLoad',[start end]: allows the user to specify how many images to
%    load, in the form of [start end] vector.
% 'FrameSize',[xpxls ypxls]: size of each frame to be loaded. 
%
% Updated 1/2014 to return an allimgs array of the same integer type as 
% the original file (in our case, to keep it as uint16).
%
% Updated 7/2014 to add a second optional input containing the xdim and ydim
% of the images to be loaded, so you don't have to call GetInfoFromMetadata
% every time if you're calling this in a loop.
%
% Note that this does not make use of the FrameLoadMax parameter in 
% smFRETsetup!  You can load as many frames as you want with this function,
% including so many it'll crash Matlab if Matlab doesn't have enough
% memory ...
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

function allimgs = LoadUManagerTifsV5(D,varargin)

    alltifs = dir(fullfile(D,'img*.tif'));
    
    % Deal with inputs
    if ~isempty(varargin)
        for p = 1:2:length(varargin)
            if strcmpi(varargin{p},'FramesToLoad')
                StartStop = sort(varargin{p+1});
                if StartStop(1)<=0
                    StartStop(1)=1;
                end
                if StartStop(2)>length(alltifs)
                    StartStop(2) = length(alltifs);
                end
            elseif strcmpi(varargin{p},'FrameSize')
                xpxls = varargin{p+1}(1);
                ypxls = varargin{p+1}(2);
            end
        end
    end
    
    % In folder D uManager will have saved a bunch of image files, and two
    % text files that contain information about the data and its acquisition.
    
    % First figure out the size of the images:
    % Update 7/2014:
    if ~exist('xpxls','var')
        val = GetInfoFromMetaData(D,'imgsize');
        xpxls = val(1);
        ypxls = val(2);
    end
    
    % Update 1/2014: Load one file to get the integer type class:
    temp = imread(fullfile(D,alltifs(1).name));
    classtype = class(temp);
    clear temp

    if ~exist('StartStop','var') || (StartStop(2)-StartStop(1)+1)>=length(alltifs)
        allimgs = zeros(xpxls,ypxls,length(alltifs),classtype);
    else
        allimgs = zeros(xpxls,ypxls,(StartStop(2)-StartStop(1)+1),classtype);
    end

    % Load all the images into a 3d matrix:
    if ~exist('StartStop','var')
        for i = 1:length(alltifs)
            img = imread(fullfile(D,alltifs(i).name));
            % img = mat2gray(img);
            allimgs(:,:,i) = img;
            clear img
        end
    else
        incr = 1;
        for i = StartStop(1):StartStop(2);
            img = imread(fullfile(D,alltifs(i).name));
            % img = mat2gray(img);
            allimgs(:,:,incr) = img;
            clear img
            incr = incr+1;
        end
    end

end