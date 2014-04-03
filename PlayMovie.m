% function PlayMovie(PathToMovie,frames,varargin)
%
% Plays a movie of smFRET data from frames(1) to frames(2), with the
% optional input allowing the user to display this movie in an existing
% figure.
%
% Inputs:
% PathToMovie: Full path to folder with the "ScaledMovieFrames..." files
% frames: [start end] vector of frames to show
%
% Stephanie 4/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function PlayMovie(PathToMovie,frames,varargin)

    % Input error handling
    frames = sort(frames);
    if frames(1)<=0
        frames(1)=1;
    end
    
    if isempty(varargin)
        h2 = figure('Position',[650,800,500,650]);
    end
    
    lastframe = frames(1);
    
    while lastframe<frames(2)
        disp('Loading movie part ... ')
        [movRed,movGreen,lastframe] = LoadScaledMovie(PathToMovie,...
            [lastframe frames(2)]);
        
        for i=1:size(movRed,3)
            if isempty(varargin)
                subplot('Position',[0.08 0.23 0.39 0.39*512/256])
                imshow(movRed(:,:,i),[])
                subplot('Position',[0.54 0.23 0.39 0.39*512/256])
                imshow(movGreen(:,:,i),[])
                drawnow
            end
        end
            
        lastframe = lastframe+1;
        clear movRed movGreen
        
    end
