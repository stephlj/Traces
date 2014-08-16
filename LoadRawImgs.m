% function [allimgs,numframes] = LoadRawImgs(D,varargin)
%
% Calls either LoadPMA or LoadUManagerTifs and returns an intensity matrix
% as allimgs, and the total number of frames as numframes.  
% More details about the inputs can be found by calling "help
% LoadPMA" and "help LoadUManagerTifs" from the command line.
%
% The default is to load the pma, if one exists.
%
% Inputs:
% D: the full path to the directory containg the pma file or the tifs
%
% Optional inputs: Pass these as '<name>',<val> pairs, in any order.
% 'FramesToLoad',[start end]: allows the user to specify how many images to
%    load, in the form of [start end] vector. Optional for both LoadPMA and
%    LoadUManagerTifs.
% 'FrameSize',[xpxls ypxls]: size of each frame to be loaded. Only an option
%    for LoadUManagerTifs.
%
% Output image matrix returned is NOT scaled between 0 and 1 
% (or scaled at all), and is returned as the same integer type as the raw
% images.
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

function [allimgs,numframes] = LoadRawImgs(D,varargin)

    % Input handling:
    if ~isempty(varargin)
        for p = 1:2:length(varargin)
            if strcmpi(varargin{p},'FramesToLoad')
                StartStop = sort(varargin{p+1});
                if StartStop(1)<=0
                    StartStop(1)=1;
                end
            elseif strcmpi(varargin{p},'FrameSize')
                xpxls = varargin{p+1}(1);
                ypxls = varargin{p+1}(2);
                FrameSize = [xpxls, ypxls];
            end
        end
    end

    % Decide whether to call LoadPMA or LoadUManagerTifs
    % irritatingly, can't use regular expressions with exist:
    templist = dir(fullfile(D,'*.pma'));
    if ~isempty(templist)
        numtype = GetInfoFromMetaData(D,'precision');
        dirinfo = dir(fullfile(D,'*.pma'));
        if ~exist('StartStop','var')
            [allimgs,numframes] = LoadPMA(fullfile(D,dirinfo.name),numtype);
        else
            [allimgs,numframes] = LoadPMA(fullfile(D,dirinfo.name),numtype,...
                'FramesToLoad',StartStop);
        end
    else
        % Get the total number of tifs
        alltifs = dir(fullfile(D,'img*.tif'));
        numframes = length(alltifs);
        if ~exist('StartStop','var') && ~exist('FrameSize','var')
            allimgs = LoadUManagerTifsV5(D);
        elseif ~exist('StartStop','var')
            allimgs = LoadUManagerTifsV5(D,'FrameSize',FrameSize);
        elseif ~exist('FrameSize','var')
            allimgs = LoadUManagerTifsV5(D,'FramesToLoad',StartStop);
        else
            allimgs = LoadUManagerTifsV5(D,'FrameSize',FrameSize,...
                'FramesToLoad',StartStop);
        end
    end
end
