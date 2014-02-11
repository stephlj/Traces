% function [MovieMin,MovieMax] = ScaleMovie(PathToMovie,numframes)
%
% Given the path to a directory with images, PathToMovie, and the total number of
% frames in this directory, numframes, calculate the min and max over the
% whole movie.
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [MovieMin,MovieMax] = ScaleMovie(PathToMovie,numframes)
 
    % First, set the initial max to the minimum possible value, then the initial minimum to
    % the maximum possible value:
    MovieMax = 0;
    % To get the maximum possible value, load one frame and figure out the
    % numeric type:
    tempfr = LoadUManagerTifsV5(PathToMovie,[1 1]);
    if strcmpi(class(tempfr),'uint16')
        MovieMin = 2^16-1;
    elseif strcmpi(class(tempfr),'uint8')
        MovieMin = 2^8-1;
    else
        MovieMin = 100000; % Random very large number!
    end
    clear tempfr;

    for jj = 1:100:numframes
        moviebit = double(LoadUManagerTifsV5(PathToMovie,[jj jj+99]));
        % Do I want to scale the entire image to the same max and min, or
        % do I want to scale the two channels separately?  The Ha lab code
        % does the whole image together, I believe.  So only calculating
        % one max and one min.
        
        MovieMax = max(MovieMax,double(max(max(max(moviebit)))));
        MovieMin = min(MovieMin,double(min(min(min(moviebit)))));
        
    end
end