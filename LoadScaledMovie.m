% function [movRed,movGreen, firstframe,lastframe] = LoadScaledMovie(PathToMovie,frames)
%
% Loads the scaled and background-subtracted movie frames saved by
% CalcIntensities.m. 
%
% Inputs:
% PathToMovie: Full path to folder with the "ScaledMovieFrames..." files
% frames: [start end] vector of frames to return. End-start must be <=100,
%   else it will return only the first 100 frames specified by the interval
%   [start start+99];
%
% Outputs:
% movRed, movGreen: x-by-y-by-frames matrix of intensities for the red and
%   green channels
% lastframe: The last frame returned.  This is frames(2) if the user asked
%   for fewer than 100 frames to be returned; else this is start+99.
%
% Stephanie 4/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [movRed,movGreen,lastframe] = LoadScaledMovie(PathToMovie,frames)

    % Input error handling
    frames = sort(frames);
    if frames(1)<=0
        frames(1)=1;
    end
    % My laptop can't handle more than 100 frames at a time
    if frames(2)-frames(1)>100
        frames(2)=frames(1)+99;
    end
    totfiles = length(dir(fullfile(PathToMovie,'ScaledMovie*.mat')));
    % Make sure the file exists
    if frames(2)>totfiles*100
        disp('LoadScaledMovie: Movie not that long?')
    end
    
    firstframe = frames(1);
    lastframe = firstframe;
    
    % For some reason, if you initialize movRed to be
    % movRed = [];
    % size(movRed,3) = 1 ... 
    movRed = zeros(0,0,0);
    movGreen = zeros(0,0,0);
    
    % Converts global frame number to per-file frame number
    FileFrame = @(frame) rem(frame - 1, 100) + 1;

    % Figure out which file(s) to load:
    while lastframe<frames(2)
        FileNum = firstframe-FileFrame(firstframe) + 1;
        
        temp = load(fullfile(PathToMovie,sprintf('ScaledMovieFrames%dto%d.mat',...
            FileNum,FileNum+99)));
        
        lastframe = min(FileNum + 99, frames(2));
        
        %disp(sprintf('movRed(:,:,%d:%d) = temp.imgR(:,:,%d:%d)', size(movRed, 3)+1, size(movRed, 3)+1+lastframe-firstframe, FileFrame(firstframe), FileFrame(lastframe)))
        movRed(  :,:,end+1:end+1+lastframe-firstframe) = temp.imgR(:,:,FileFrame(firstframe):FileFrame(lastframe));
        movGreen(:,:,end+1:end+1+lastframe-firstframe) = temp.imgG(:,:,FileFrame(firstframe):FileFrame(lastframe));
        
        firstframe = lastframe+1;
        
        clear temp
        
    end
end
    