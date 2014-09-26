% function ScaleMovie(PathToMovie,params)
%
% Given the path to a directory with images (PathToMovie) and a parameter
% structure created by smFRETsetup (params), calculate the min and max over the
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

function ScaleMovieV2(PathToMovie,params)

        % Subfunction to check for and remove intensity spikes:
        function newMax = CheckMaxes(Maxes,Mins,oldMax,xaxislim,NormImage)
            sortedMaxes = sort(Maxes); % Sort defaults to ascending order
            MaxDiffs = sortedMaxes(2:end)-sortedMaxes(1:end-1);
            meanDiff = mean(MaxDiffs);
            stdDiff = std(MaxDiffs);
            newMax = oldMax;
            while (sortedMaxes(end)-sortedMaxes(end-1))>(meanDiff+3*stdDiff) ||...
                    newMax>(mean(Maxes)+5*std(Maxes))
                newMax = sortedMaxes(end-1);
                sortedMaxes = sortedMaxes(1:end-1);
                if length(sortedMaxes)<(length(Maxes)-7)
                    disp('ScaleMovieV2: A lot of max outliers?')
                    figure
                    plot(1:xaxislim,Maxes,'ob',1:xaxislim,Mins,'xr')
                    if NormImage
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
        end

    % Get some preliminary info:
    framesize = GetInfoFromMetaData(PathToMovie,'imgsize');
 
    % First, set the initial max to the minimum possible value, then the initial minimum to
    % the maximum possible value:
    MovieMax = 0;
    % To get the maximum possible value, load one frame and figure out the
    % numeric type:
    [tempfr,numframes] = LoadRawImgs(PathToMovie,'FramesToLoad',[1 1]);
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
        
    % Update 8/2014: Now allowing the option of scaling each channel
    % separately. Note that because of the random spikes in maxes that we
    % currently observe, I have to record the max for each channel in each
    % frame ... 
    allGrMaxes = zeros(1,numframes);
    allRedMaxes = zeros(1,numframes);
    
    if params.FrameLoadMax > numframes
        FrameLoadMax = numframes;
    else
        FrameLoadMax = params.FrameLoadMax;
    end
    
    % Using SplitImg below to get red and green channels--note that this will
    % clip the image by params.PxlsToExclude before calculating the min
    % and the max, unless I mess with the params structure a bit before
    % calling SplitImg:
    realPxlsToExclude = params.PxlsToExclude;
    params.PxlsToExclude = 0;

    for jj = 1:FrameLoadMax:numframes
        [moviebit,~] = LoadRawImgs(PathToMovie,'FramesToLoad',[jj jj+FrameLoadMax-1],...
            'FrameSize',framesize);
        moviebit = double(moviebit);
        
        [redbit,grbit] = SplitImg(moviebit,params);
        
        % Update 4/2014: Allowing a normalization option--see note about
        % NormImage in smFRETsetup.m. (The real normalization is done
        % in LoadScaledMovie, this just makes sure MovieMax and MovieMin are
        % appropriately scaled if necessary):
        rmoviebit = reshape(moviebit,size(moviebit,1)*size(moviebit,2),size(moviebit,3));
        allMedians(jj:jj+size(moviebit,3)-1) = median(rmoviebit,1);
        if params.NormImage
            moviebit = bsxfun(@rdivide,moviebit,...
                reshape(allMedians(jj:jj+size(moviebit,3)-1),1,1,size(moviebit,3)));
            redbit = bsxfun(@rdivide,redbit,...
                reshape(allMedians(jj:jj+size(redbit,3)-1),1,1,size(redbit,3)));
            grbit = bsxfun(@rdivide,grbit,...
                reshape(allMedians(jj:jj+size(grbit,3)-1),1,1,size(grbit,3)));
        end
        
        MovieMax = max(MovieMax,double(max(moviebit(:))));
        MovieMin = min(MovieMin,double(min(moviebit(:))));
        
        MovieMaxGr = max(MovieMax,double(max(grbit(:))));
        MovieMinGr = min(MovieMin,double(min(grbit(:))));
        MovieMaxRed = max(MovieMax,double(max(redbit(:))));
        MovieMinRed = min(MovieMin,double(min(redbit(:))));

        allMaxes(jj:jj+size(moviebit,3)-1) = max(max(moviebit,[],2),[],1);
        allMins(jj:jj+size(moviebit,3)-1) = min(min(moviebit,[],2),[],1);
        
        allGrMaxes(jj:jj+size(moviebit,3)-1) = max(max(grbit,[],2),[],1);
        allRedMaxes(jj:jj+size(moviebit,3)-1) = max(max(redbit,[],2),[],1);
        allGrMins(jj:jj+size(moviebit,3)-1) = min(min(grbit,[],2),[],1);
        allRedMins(jj:jj+size(moviebit,3)-1) = min(min(redbit,[],2),[],1);
        
        clear moviebit redbit grbit
    end
    
    params.PxlsToExclude = realPxlsToExclude;
    clear realPxlsToExclude
    
    % Check that the calculated Max(es) isn't/aren't way larger than the rest of the
    % maxes (it happens):
    MovieMax = CheckMaxes(allMaxes,allMins,MovieMax,numframes,params.NormImage);
    MovieMaxGr = CheckMaxes(allGrMaxes,allGrMins,MovieMaxGr,numframes,params.NormImage);
    MovieMaxRed = CheckMaxes(allRedMaxes,allRedMins,MovieMaxRed,numframes,params.NormImage);
    
    % Plot a figure so the user can check whether normalization was
    % necessary or not:
    if ~params.ScaleChannelsSeparately
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
    else
        disp(sprintf('Using acceptor max=%d, min=%d',MovieMaxRed,MovieMinRed))
        disp(sprintf('Using donor max=%d, min=%d',MovieMaxGr,MovieMinGr))
        
        figure
        subplot(3,1,1)
        plot(1:numframes,allMaxes,'ob',1:numframes,allMins,'xr')
        legend('Maxes','Mins')
        set(gca,'FontSize',14)
        xlabel('Frame','FontSize',14)
        if params.NormImage
            ylabel('Normalized Intensity (a.u.)','FontSize',14)
        else
            ylabel('Raw intensity (a.u.)','FontSize',14)
        end
        subplot(3,1,2)
        plot(1:numframes,allRedMaxes,'or',1:numframes,allGrMaxes,'og',...
            1:numframes,allRedMins,'xr',1:numframes,allGrMins,'xg')
        legend('Acceptor maxes','Donor maxes','Acceptor mins','Donor mins')
        set(gca,'FontSize',14)
        xlabel('Frame','FontSize',14)
        if params.NormImage
            ylabel('Normalized Intensity (a.u.)','FontSize',14)
        else
            ylabel('Raw intensity (a.u.)','FontSize',14)
        end
        subplot(3,1,3)
        plot(1:numframes,allMedians,'.g')
        legend('Median')
        set(gca,'FontSize',14)
        xlabel('Frame','FontSize',14)
        ylabel('Raw Intensity (a.u.)','FontSize',14)
        print('-depsc',fullfile(PathToMovie,'ScalingFig'))
    end
    pause
    close
    
    NormImage = params.NormImage;
    ScaleChannelsSeparately = params.ScaleChannelsSeparately;
    % Save the scaling information to disk for use by LoadScaledMovie:
    save(fullfile(PathToMovie,strcat('ScalingInfo.mat')),'allMedians',...
        'MovieMin','MovieMax','MovieMinRed','MovieMaxRed','MovieMinGr','MovieMaxGr',...
        'NormImage','ScaleChannelsSeparately')
    clear NormImage ScaleChannelsSeparately
end