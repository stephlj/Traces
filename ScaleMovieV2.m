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
    
    allMedians = zeros(1,numframes);
    allMaxes = zeros(1,numframes);
    % Just for visualization:
    allMins = zeros(1,numframes);
        
    % For debugging
    %totI = zeros(1,numframes);
    %allMeans = zeros(1,numframes);

    for jj = 1:100:numframes
        moviebit = double(LoadUManagerTifsV5(PathToMovie,[jj jj+99]));
        % Do I want to scale the entire image to the same max and min, or
        % do I want to scale the two channels separately?  The Ha lab code
        % does the whole image together, I believe.  So only calculating
        % one max and one min.
        
        % Update 4/2014: Allowing a normalization option--see note about
        % NormImage in smFRETsetup.m. (The real normalization is done
        % below, this just makes sure MovieMax and MovieMin are
        % appropriately scaled if necessary):
        allMedians(jj:jj+size(moviebit,3)-1) = median(median(moviebit,2),1);
        if params.NormImage
            moviebit = moviebit./repmat(median(median(moviebit,2),1),...
                size(moviebit,1),size(moviebit,2),1);
        end
        
        MovieMax = max(MovieMax,double(max(max(max(moviebit)))));
        MovieMin = min(MovieMin,double(min(min(min(moviebit)))));

        allMaxes(jj:jj+size(moviebit,3)-1) = max(max(moviebit,[],2),[],1);
        allMins(jj:jj+size(moviebit,3)-1) = min(min(moviebit,[],2),[],1);
        
        % For debugging: 
        %totIR(jj:jj+size(moviebit,3)-1) = sum(sum(moviebit,2),1);
        %allMeans(jj:jj+size(moviebit,3)-1) = sum(sum(moviebit,2),1)./(size(moviebit,1)+size(moviebit,2));
        
        clear moviebit
    end
    % Check that the calculated Max isn't way larger than the rest of the
    % maxes (it happens):
    sortedMaxes = sort(allMaxes); % Sort defaults to ascending order
    MaxDiffs = sortedMaxes(2:end)-sortedMaxes(1:end-1);
    meanDiff = mean(MaxDiffs);
    stdDiff = std(MaxDiffs);
    while (sortedMaxes(end)-sortedMaxes(end-1))>(meanDiff+3*stdDiff) ||...
            MovieMax>(mean(allMaxes)+6*std(allMaxes))
        MovieMax = sortedMaxes(end-1);
        sortedMaxes = sortedMaxes(1:end-1);
        if length(sortedMaxes)<(length(allMaxes)-5)
            disp('ScaleMovieV2: A lot of max outliers?')
            figure
            plot(1:numframes,allMaxes,'ob',1:numframes,allMins,'xr')
            keyboard
        end
    end
    clear sortedMaxes MaxDiffs meanDiff stdDiff
    
    % Plot a figure so the user can check whether normalization was
    % necessary or not:
    disp(sprintf('Using max=%d, min=%d',MovieMax,MovieMin))
    figure
    subplot(2,1,1)
    plot(1:numframes,allMaxes,'ob',1:numframes,allMins,'xr')
    legend('Maxes','Mins')
    xlabel('Frame')
    if params.NormImage
        ylabel('Normalized Intensity (a.u.)')
    else
        ylabel('Raw intensity (a.u.)')
    end
    subplot(2,1,2)
    plot(1:numframes,allMedians,'.g')
    legend('Median')
    xlabel('Frame')
    ylabel('Raw Intensity (a.u.)')
    pause
    close
    
    % Re-load everything and actually do the scaling:
    for jj = 1:100:numframes
        moviebit = double(LoadUManagerTifsV5(PathToMovie,[jj jj+99]));
                
        % Update 4/2014: Allowing a normalization option--see note about
        % NormImage in smFRETsetup.m
        if params.NormImage
            tempMeds(1,1,:) = allMedians(jj:jj+size(moviebit,3)-1);
            moviebit = moviebit./repmat(tempMeds,size(moviebit,1),size(moviebit,2),1);
            clear tempMeds
        end
        
        ScaledMovie = mat2gray(moviebit,[MovieMin MovieMax]); % This also converts it to double precision,
            % but need to explicitly do so earlier in case NormImage is 1
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
        % HACK saving this info
        %save(fullfile(PathToMovie,'MovieTotIntensities.mat'),'totIR','totIG')
    end
end