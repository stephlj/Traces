% function [movRed,movGreen,movRedBkgnd,movGrBkgnd,firstframe,lastframe] = LoadScaledMovie(PathToMovie,frames)
%
% Loads the scaled and background-subtracted movie frames saved by
% ScaleMovie.m. 
%
% Inputs:
% PathToMovie: Full path to folder with the "ScaledMovieFrames..." files
% frames: [start end] vector of frames to return. 
% Optional: you must also pass 'bkgnd' as the last input, if you want it
%   to return the third and fourth outputs (see below).
%
% Outputs:
% movRed, movGreen: x-by-y-by-frames matrix of intensities for the red and
%   green channels
% movRedBkgnd, movGrBkgnd: images used to compute background for intensity
%   calculation (currently only called by CalcIntensitiesV3)
% lastframe: The last frame returned.  This is frames(2) if the user asked
%   for fewer than 100 frames to be returned; else this is start+99.
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

function [movRed,movGreen,movRedBkgnd,movGrBkgnd,lastframe] = LoadScaledMovie(PathToMovie,frames,varargin)

    % Input error handling
    frames = sort(frames);
    if frames(1)<=0
        frames(1)=1;
    end
    if frames(1)==frames(2)
        frames(2)=frames(1)+1;
    end
    % Make sure all files and frames exist
    [~,totframes] = LoadRawImgs(PathToMovie,'FramesToLoad',[1 1]);
    if frames(1) > totframes
        disp('LoadScaledMovie: Movie not that long?')
        movRed = -1;
        movGreen = -1;
        lastframe = -1;
        movRedBkgnd = -1;
        movGrBkgnd = -1;
        return
    end
    if frames(2) > totframes
        frames(2) = totframes;
        if frames(1) == frames(2)
            frames(2) = frames(2)+1;
        end
    end
    
    % Make sure the scaling information is available:
    if ~exist(fullfile(PathToMovie,strcat('ScalingInfo.mat')),'file')
        disp('LoadScaledMovie: Scaling info not saved?')
        movRed = -1;
        movGreen = -1;
        lastframe = -1;
        movRedBkgnd = -1;
        movGrBkgnd = -1;
        return
    else
        load(fullfile(PathToMovie,strcat('ScalingInfo.mat')));
    end

    totfiles = length(allfiles);
    % If there's only one file, life is pretty easy:
    if totfiles == 1;
        if frames(1) > FrameSaveIncr
            disp('LoadScaledMovie: Movie not that long?')
            movRed = -1;
            movGreen = -1;
            lastframe = -1;
            movRedBkgnd = -1;
            movGrBkgnd = -1;
            return
        end
        if isempty(varargin)
            temp = load(fullfile(PathToMovie,strcat('ScaledMovieFrames1to',int2str(FrameSaveIncr),'.mat')),'imgR','imgG');
        else
            temp = load(fullfile(PathToMovie,strcat('ScaledMovieFrames1to',int2str(FrameSaveIncr),'.mat')));
        end
        lastframe = min(frames(2),size(temp.imgR,3)); %Can't load frames that don't exist, even if the user asked for them
        movRed = temp.imgR(:,:,frames(1):lastframe);
        movGreen = temp.imgG(:,:,frames(1):lastframe);
        if ~isempty(varargin)
            movRedBkgnd = temp.imgRBkgnd(:,:,frames(1):lastframe);
            movGrBkgnd = temp.imgGBkgnd(:,:,frames(1):lastframe);
        else
            movRedBkgnd = [];
            movGrBkgnd = [];
        end
        return
    end
    
    firstframe = frames(1);
    lastframe = firstframe;
    
    % For some reason, if you initialize movRed to be
    % movRed = [];
    % size(movRed,3) = 1 ... 
    movRed = zeros(0,0,0);
    movGreen = zeros(0,0,0);
    movRedBkgnd = zeros(0,0,0);
    movGrBkgnd = zeros(0,0,0);
    
    % Converts global frame number to per-file frame number
    FileFrame = @(frame) rem(frame - 1, FrameSaveIncr) + 1;

    % Figure out which file(s) to load:
    while lastframe<frames(2)
        FileNum = firstframe-FileFrame(firstframe) + 1;
        
        if isempty(varargin)
            temp = load(fullfile(PathToMovie,sprintf('ScaledMovieFrames%dto%d.mat',...
                FileNum,FileNum+FrameSaveIncr-1)),'imgR','imgG');
        else
            temp = load(fullfile(PathToMovie,sprintf('ScaledMovieFrames%dto%d.mat',...
                FileNum,FileNum+FrameSaveIncr-1)));
        end
        
        lastframe = min(FileNum + FrameSaveIncr-1, frames(2));
        
        %disp(sprintf('movRed(:,:,%d:%d) = temp.imgR(:,:,%d:%d)', size(movRed, 3)+1, size(movRed, 3)+1+lastframe-firstframe, FileFrame(firstframe), FileFrame(lastframe)))
        movRed(  :,:,end+1:end+1+lastframe-firstframe) = temp.imgR(:,:,FileFrame(firstframe):FileFrame(lastframe));
        movGreen(:,:,end+1:end+1+lastframe-firstframe) = temp.imgG(:,:,FileFrame(firstframe):FileFrame(lastframe));
        if ~isempty(varargin)
            movRedBkgnd(:,:,end+1:end+1+lastframe-firstframe) = temp.imgRBkgnd(:,:,FileFrame(firstframe):FileFrame(lastframe));
            movGrBkgnd (:,:,end+1:end+1+lastframe-firstframe) = temp.imgGBkgnd(:,:,FileFrame(firstframe):FileFrame(lastframe));
        else
            movRedBkgnd = [];
            movGrBkgnd = [];
        end
        firstframe = lastframe+1;
        
        clear temp
        
    end
end
    