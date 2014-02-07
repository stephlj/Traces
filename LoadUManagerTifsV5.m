% function allimgs = LoadUManagerTifsV5(D)
%
% Loads all the images in a uManager-created folder into one 3-d matrix.
% D is the full path to the folder where they're stored; output is the
% image matrix.  In V5, images returned are NOT scaled between 0 and 1.
% Also, the optional input allows the user to specify how many images to
% load, in the form of [start end] vector.
%
% Updated 1/2014 to return an allimgs array of the same integer type as 
% the original file (in our case, to keep it as uint16).
%
% Steph 3/2013, updated 10/2013
% Copyright 2013 Stephanie Johnson, University of California, San Francisco

function allimgs = LoadUManagerTifsV5(D,varargin)

    %Input error handling
    if ~isempty(varargin)
        StartStop = sort(varargin{1});
        if StartStop(1)<=0
            StartStop(1)=1;
        end
    end

    %In folder D uManager will have saved a bunch of image files, and two
    %text files that contain information about the data and its acquisition.
    
    %First figure out the size of the images:
    val = GetInfoFromMetaData(D,'imgsize');
    xpxls = val(1);
    ypxls = val(2);
    
    %Next figure out how many tif files there are:
    
    alltifs = dir(fullfile(D,'img*.tif'));
    
    %Update 1/2014: Load one file to get the integer type class:
    temp = imread(fullfile(D,alltifs(1).name));
    classtype = class(temp);
    clear temp

    if isempty(varargin) || (StartStop(2)-StartStop(1)+1)>=length(alltifs)
        allimgs = zeros(xpxls,ypxls,length(alltifs),classtype);
    else
        allimgs = zeros(xpxls,ypxls,(StartStop(2)-StartStop(1)+1),classtype);
    end

    %Load all the bead images into a 3d matrix:
    if isempty(varargin)
        for i = 1:length(alltifs)
            img = imread(fullfile(D,alltifs(i).name));
            %img = mat2gray(img);
            allimgs(:,:,i) = img;
            clear img
        end
    else
        incr = 1;
        if StartStop(2)>length(alltifs)
            StartStop(2) = length(alltifs);
        end
        for i = StartStop(1):StartStop(2);
            img = imread(fullfile(D,alltifs(i).name));
            %img = mat2gray(img);
            allimgs(:,:,incr) = img;
            clear img
            incr = incr+1;
        end
    end

end