% function PlayMovie(PathToMovie,frames,params,varargin)
%
% Plays a movie of smFRET data from frames(1) to frames(2), with the
% optional input allowing the user to display this movie in an existing
% figure.
%
% Inputs:
% PathToMovie: Full path to folder with the image files
% frames: [start end] vector of frames to show
% params: parameter structure saved by smFRETsetup
% Optional inputs: must have either the first three only, or all seven
%   varargin{1}: handle to a figure to play the movie into
%   varargin{2}: a string containing 'subplot(blah)' to play the red
%       channel into
%   varargin{3}: same as {2} but for green channel
%   varargin{4}: If a movie of a zoom-in on a particular spot is also
%       desired, this must contain an [x;y] vector of the spot's location
%       in the red channel
%   varargin{5}: same as {4} but for green channel
%   varargin{6}: If {4} and {5} are passed, then a string of
%       'subplot(blah)' must also be passed for where to plot the zoomed
%       movie in the red channel
%   varargin{7}: same as {6} but for green channel
%   varargin{8}: size of the ROI to extract (one side of a square)
%   varargin{9}: if {4} and {5} are passed, [xdim;ydim] of size of circle 
%       to put on the zoomed image
%
% Outputs:
% LastRedFrame = last frame of the red channel that was played
% LastGreenFrame = last frame of the green channel
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

function [LastRedFrame,LastGreenFrame] = PlayMovie(PathToMovie,frames,params,varargin)

    % Input error handling
    frames = sort(frames);
    if frames(1)<=0
        frames(1)=1;
    end
    if isempty(varargin)
        h2p = figure('Position',[650,800,500,650]);
    else 
        if length(varargin)~=3 && length(varargin)~=9
            disp('PlayMovie: Optional input must contain either 3 or 9 elements.');
            LastRedFrame=-1;
            LastGreenFrame = -1;
            return
        else
            h2p = varargin{1};
            figure(h2p)
        end
    end
    % subfunction for putting circles around a spot:
    function boxfun(currspot,circlesize,markercolor)
        % CalcSpotIntensityNoGauss puts a circle of diameter 5 around each spot:
        t = 0:pi/100:2*pi;
        plot(currspot(2)+circlesize(2)/2.*cos(t),currspot(1)+circlesize(1)/2.*sin(t),strcat('-',markercolor))
        clear t
    end
    
    lastframe = frames(1);
    
    while lastframe<=frames(2)
        disp('Loading movie part ... ')
        [movRed,movGreen,~,~,lastframe] = LoadScaledMovie(PathToMovie,...
            [lastframe min(lastframe+params.FrameLoadMax, frames(2))],params);
        
        for i=1:size(movRed,3)
            if isempty(varargin)
                subplot('Position',[0.08 0.23 0.39 0.39*512/256])
                imshow(movRed(:,:,i))
                subplot('Position',[0.54 0.23 0.39 0.39*512/256])
                imshow(movGreen(:,:,i))
                drawnow
            else
                eval(varargin{2})
                imshow(movRed(:,:,i))
                hold on
                boxfun(varargin{4},varargin{9},'r');
                hold off
                title('Red','Fontsize',12)
                eval(varargin{3})
                imshow(movGreen(:,:,i))
                hold on
                boxfun(varargin{5},varargin{9},'g');
                hold off
                title('Green','Fontsize',12)
                if length(varargin)>3
                    [imgRzoom,zoomcenR] = ExtractROI(movRed(:,:,i),varargin{8},varargin{4});
                    [imgGzoom,zoomcenG] = ExtractROI(movGreen(:,:,i),varargin{8},varargin{5});
                    eval(varargin{6})
                    imshow(imgRzoom)
                    hold on
                    boxfun(zoomcenR,varargin{9},'r');
                    hold off
                    eval(varargin{7})
                    imshow(imgGzoom)
                    hold on
                    boxfun(zoomcenG,varargin{9},'g');
                    hold off
                end
                drawnow
            end
            clear imgRzoom imgGzoom zoomcenR zoomcenG
        end
            
        lastframe = lastframe+1;
        LastRedFrame = movRed(:,:,end);
        LastGreenFrame = movGreen(:,:,end);
        clear movRed movGreen
        
    end
end
