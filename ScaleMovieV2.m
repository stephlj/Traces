% function ScaleMovie(PathToMovie,numframes,params)
%
% Given the path to a directory with images (PathToMovie) and the total number of
% frames in this directory (numframes) calculate the min and max over the
% whole movie, and then scale the whole movie between this min and max.
% Saves it in FrameLoadMax-frame increments to the same folder as PathToMovie.
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

function ScaleMovieV2(PathToMovie,numframes,params)

    % Get some preliminary info:
    framesize = GetInfoFromMetaData(PathToMovie,'imgsize');
 
    % First, set the initial max to the minimum possible value, then the initial minimum to
    % the maximum possible value:
    MovieMax = 0;
    % To get the maximum possible value, load one frame and figure out the
    % numeric type:
    tempfr = LoadRawImgs(PathToMovie,'FramesToLoad',[1 1]);
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
    
    if params.FrameLoadMax > numframes
        FrameLoadMax = numframes;
    else
        FrameLoadMax = params.FrameLoadMax;
    end

    for jj = 1:FrameLoadMax:numframes
        moviebit = double(LoadRawImgs(PathToMovie,'FramesToLoad',[jj jj+FrameLoadMax-1],...
            'FrameSize',framesize));
        % Do I want to scale the entire image to the same max and min, or
        % do I want to scale the two channels separately?  The Ha lab code
        % does the whole image together, I believe.  So only calculating
        % one max and one min.
        
        % Update 4/2014: Allowing a normalization option--see note about
        % NormImage in smFRETsetup.m. (The real normalization is done
        % below, this just makes sure MovieMax and MovieMin are
        % appropriately scaled if necessary):
        rmoviebit = reshape(moviebit,size(moviebit,1)*size(moviebit,2),size(moviebit,3));
        allMedians(jj:jj+size(moviebit,3)-1) = median(rmoviebit,1);
        if params.NormImage
            moviebit = bsxfun(@rdivide,moviebit,...
                reshape(allMedians(jj:jj+size(moviebit,3)-1),1,1,size(moviebit,3)));
        end
        
        MovieMax = max(MovieMax,double(max(moviebit(:))));
        MovieMin = min(MovieMin,double(min(moviebit(:))));

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
            if params.NormImage
                ylabel('Normalized Intensity (a.u.)','FontSize',14)
            else
                ylabel('Raw intensity (a.u.)','FontSize',14)
            end
            legend('Maxes','Mins')
            set(gca,'FontSize',14)
            xlabel('Frame','FontSize',14)
            keyboard
            close
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
    set(gca,'FontSize',14)
    xlabel('Frame','FontSize',14)
    if params.NormImage
        ylabel('Normalized Intensity (a.u.)','FontSize',14)
    else
        ylabel('Raw intensity (a.u.)','FontSize',14)
    end
    subplot(2,1,2)
    plot(1:numframes,allMedians,'.g')
    legend('Median')
    set(gca,'FontSize',14)
    xlabel('Frame','FontSize',14)
    ylabel('Raw Intensity (a.u.)','FontSize',14)
    print('-depsc',fullfile(PathToMovie,'ScalingFig'))
    pause
    close
    
    disp('Continuing with the scaling ...')
    
    % Re-load everything and actually do the scaling:
    profile on
    for jj = 1:FrameLoadMax:numframes
        moviebit = double(LoadRawImgs(PathToMovie,'FramesToLoad',[jj jj+FrameLoadMax-1],...
            'FrameSize',framesize));
                
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

        [imgRBkgnd,imgGBkgnd] = CalcBkgnd(imgR,imgG,params);
        
        save(fullfile(PathToMovie,strcat('ScaledMovieFrames',int2str(jj),...
            'to',int2str(jj+FrameLoadMax-1),'.mat')),'imgR','imgG','imgRBkgnd','imgGBkgnd')
        
        clear moviebit imgR imgG ScaledMovie imgRBkgnd imgGBkgnd
        
    end
    keyboard
end