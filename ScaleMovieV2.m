% function ScaleMovie(PathToMovie,numframes,params)
%
% Given the path to a directory with images, PathToMovie, and the total number of
% frames in this directory, numframes, calculate the min and max over the
% whole movie, and then scale the whole movie between this min and max.
% Saves it in 100-frame increments to the same folder as PathToMovie.
%
% Steph 2/2014, updated 4/2014 to actually do the scaling in this function
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function ScaleMovieV2(PathToMovie,numframes,params)
 
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
        clear moviebit
    end
    % Re-load everything and actually do the scaling:
    for jj = 1:100:numframes
        moviebit = LoadUManagerTifsV5(PathToMovie,[jj jj+99]);
        ScaledMovie = mat2gray(moviebit,[MovieMin MovieMax]); % This also converts it to double precision
        [imgR,imgG] = SplitImg(ScaledMovie,params);
        
        % If I were going to do background subtraction, I should do it
        % here, but I think this doesn't work well:
        % [imgRraw,imgGraw] = SplitImg(moviebit,params);
            % Subtract time-local background: not actually sure this is the right
            % thing to do. And if I do this I should probably scale movie after?
            % Commenting out for now.
            % For debugging:
            %[imgRbkgnd,imgGbkgnd,imgR,imgG] = SubBkgnd(mean(imgRraw,3),mean(imgGraw,3),params,1);
        %     [imgRbkgnd,imgGbkgnd,~,~] = SubBkgnd(mean(imgRraw,3),mean(imgGraw,3),params);
        %     imgR = imgRraw-repmat(imgRbkgnd,1,1,size(imgRraw,3));
        %     imgG = imgGraw-repmat(imgGbkgnd,1,1,size(imgRraw,3));

            % For debugging the background subtraction:
        %     temp = reshape(imgRraw,1,size(imgRraw,1)*size(imgRraw,2)*size(imgRraw,3));
        %     hist(temp,1000)
        %     xlim([0 1])
        %     title('Red, Raw')
        %     clear temp
        %     temp = reshape(imgRbkgnd,1,size(imgRbkgnd,1)*size(imgRbkgnd,2));
        %     figure,hist(temp,1000)
        %     xlim([0 1])
        %     title('Red, Bkgnd')
        %     clear temp
        %     temp = reshape(imgR,1,size(imgR,1)*size(imgR,2)*size(imgR,3));
        %     figure,hist(temp,1000)
        %     xlim([0 1])
        %     title('Red, Minus bkgnd')
        %     clear temp
        %     pause 
        %     close all
        
        save(fullfile(PathToMovie,strcat('ScaledMovieFrames',int2str(jj),...
        'to',int2str(jj+99),'.mat')),'imgR','imgG')
        
        clear moviebit imgR imgG ScaledMovie
    end
end