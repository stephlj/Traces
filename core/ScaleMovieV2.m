% function params = ScaleMovie(PathToMovie,params)
%
% Given the path to a directory with images (PathToMovie) and a parameter
% structure created by smFRETsetup (params), calculate the min and max over the
% whole movie, and then scale the whole movie between this min and max.
% Saves it in FrameLoadMax-frame increments to the same folder as PathToMovie.
%
% The MIT License (MIT)
% 
% Copyright (c) 2014 Stephanie Johnson, University of California, San Francisco
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function params = ScaleMovieV2(PathToMovie,params)

        % Subfunction to check for and remove intensity spikes:
        function newMax = CheckMaxes(Maxes,Mins,oldMax,xaxislim,NormImage)
            sortedMaxes = sort(Maxes); % Sort defaults to ascending order
            MaxDiffs = sortedMaxes(2:end)-sortedMaxes(1:end-1);
            meanDiff = mean(MaxDiffs);
            stdDiff = std(MaxDiffs);
            newMax = oldMax;
            stdDiffmultiplier = 3;
            meanmultiplier = 4;
            while (sortedMaxes(end)-sortedMaxes(end-1))>(meanDiff+stdDiffmultiplier*stdDiff) ||...
                    newMax>(mean(Maxes)+meanmultiplier*std(Maxes))
                newMax = sortedMaxes(end-1);
                sortedMaxes = sortedMaxes(1:end-1);
                if length(sortedMaxes)<(length(Maxes)-7)
                    disp('ScaleMovieV2: A lot of max outliers?')
                    disp('To continue reducing the max value, enter: dbcont')
                    disp('To stop here, enter: stdDiffmultiplier = stdDiffmultiplier+1; meanmultiplier = meanmultiplier+1; dbcont')
                    figure
                    plot(1:xaxislim,Maxes,'ob',1:xaxislim,Mins,'xr',...
                        1:xaxislim,ones(1,xaxislim).*newMax,'--k')
                    if NormImage
                        ylabel('Normalized Intensity (a.u.)','FontSize',14)
                    else
                        ylabel('Raw intensity (a.u.)','FontSize',14)
                    end
                    legend('99.99 percentiles','Mins','Max to scale to')
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
    MovieMaxGr = MovieMax;
    MovieMaxRed = MovieMax;
    MovieMinGr = MovieMin;
    MovieMinRed = MovieMin;
    
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
    % Update 1/2015: Just looking at the maxes doesn't help very much,
    % because if you histogram all the intensity values for a movie, it has
    % a very very long tail of high values. So now looking at intensity
    % histograms for each movie, and calculating the 99.99 percentile of
    % intensity values, instead of the maxima.
    rinthist = zeros(1000,1);
    grinthist = zeros(1000,1);
    allinthist = zeros(1000,1);
    allprctiles = zeros(numframes,1);
    rprctiles = zeros(numframes,1);
    gprctiles = zeros(numframes,1);
    
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
        
        MovieMaxGr = max(MovieMaxGr,double(max(grbit(:))));
        MovieMinGr = min(MovieMinGr,double(min(grbit(:))));
        MovieMaxRed = max(MovieMaxRed,double(max(redbit(:))));
        MovieMinRed = min(MovieMinRed,double(min(redbit(:))));

        allMaxes(jj:jj+size(moviebit,3)-1) = max(max(moviebit,[],2),[],1);
        allMins(jj:jj+size(moviebit,3)-1) = min(min(moviebit,[],2),[],1);
        
        allGrMaxes(jj:jj+size(moviebit,3)-1) = max(max(grbit,[],2),[],1);
        allRedMaxes(jj:jj+size(moviebit,3)-1) = max(max(redbit,[],2),[],1);
        allGrMins(jj:jj+size(moviebit,3)-1) = min(min(grbit,[],2),[],1);
        allRedMins(jj:jj+size(moviebit,3)-1) = min(min(redbit,[],2),[],1);
                
        % Update 1/2015: Reshaping red and green to calculate prctile:
        rredbit = reshape(redbit,size(redbit,1)*size(redbit,2),size(redbit,3));
        rgrbit = reshape(grbit,size(grbit,1)*size(grbit,2),size(grbit,3));
        clear rmoviebit
        rmoviebit = reshape(moviebit,size(moviebit,1)*size(moviebit,2),size(moviebit,3));
        % prctile works on columns if given a matrix:
        allprctiles(jj:jj+size(moviebit,3)-1) = prctile(rmoviebit,99.99);
        rprctiles(jj:jj+size(moviebit,3)-1) = prctile(rredbit,99.99);
        gprctiles(jj:jj+size(moviebit,3)-1) = prctile(rgrbit,99.99);
        clear rredbit rgrbit rmoviebit
        
        % Update 1/2015: Now reshaping the red and green intensities so
        % that I can histogram them, and re-reshaping moviebit so I can get
        % an overall histogram as well:
        rmoviebit = reshape(moviebit,size(moviebit,1)*size(moviebit,2)*size(moviebit,3),1);
        rredbit = reshape(redbit,size(redbit,1)*size(redbit,2)*size(redbit,3),1);
        rgrbit = reshape(grbit,size(grbit,1)*size(grbit,2)*size(grbit,3),1);
        if jj==1
            [tempallinthist,xoutall] = hist(rmoviebit',1000);
            [temprinthist,xoutr] = hist(rredbit',1000);
            [tempgrinthist,xoutg] = hist(rgrbit',1000);
        else
            [tempallinthist,~] = hist(rmoviebit',xoutall);
            [temprinthist,~] = hist(rredbit',xoutr);
            [tempgrinthist,~] = hist(rgrbit',xoutg);
        end
        allinthist = allinthist+tempallinthist';
        rinthist = rinthist+temprinthist';
        grinthist = grinthist+tempgrinthist';
        clear temprinthist tempgrinthist tempallinthist
        
        clear moviebit redbit grbit rredbit rgrbit rmoviebit
    end
    
    params.PxlsToExclude = realPxlsToExclude;
    clear realPxlsToExclude
    
    % Update 4/2015: Taking advantage of all these maxes and medians
    % already calculated here, by detecting any red flashes here:
    if params.DetectRedFlash>0
        params = DetectRedFlash(rprctiles,gprctiles,allMedians,params);
        % These artificial maxes (from flashing the red laser) mess with
        % the calculation of the true max to which to scale the image. So
        % removing them here:
        for ss = 1:length(params.InjectPoints)
            allprctiles(params.InjectPoints(ss)-3:params.InjectPoints(ss)+3) = median(allprctiles).*ones(7,1);
            rprctiles(params.InjectPoints(ss)-3:params.InjectPoints(ss)+3) = median(rprctiles).*ones(7,1);
            allMedians(params.InjectPoints(ss)-3:params.InjectPoints(ss)+3) = median(allMedians).*ones(7,1);
        end
    end
    
    % Check that the calculated Max(es) isn't/aren't way larger than the rest of the
    % maxes (it happens):
    % MovieMax = CheckMaxes(allMaxes,allMins,MovieMax,numframes,params.NormImage);
    % MovieMaxGr = CheckMaxes(allGrMaxes,allGrMins,MovieMaxGr,numframes,params.NormImage);
    % MovieMaxRed = CheckMaxes(allRedMaxes,allRedMins,MovieMaxRed,numframes,params.NormImage);
    % Think I need to do the same with the percentiles. Note it seems to
    % work better to take the maxes instead of a percentile of the
    % percentiles.
    MovieMax = max(allprctiles);
    MovieMaxGr = max(gprctiles);
    MovieMaxRed = max(rprctiles);
    disp('Checking for outliers in overall percentiles ... ')
    MovieMax = CheckMaxes(allprctiles,allMins,MovieMax,numframes,params.NormImage);
    if params.ScaleChannelsSeparately
        disp('Checking for outliers in donor percentiles ... ')
        MovieMaxGr = CheckMaxes(gprctiles,allGrMins,MovieMaxGr,numframes,params.NormImage);
        disp('Checking for outliers in acceptor percentiles ... ')
        MovieMaxRed = CheckMaxes(rprctiles,allRedMins,MovieMaxRed,numframes,params.NormImage);
    end
    
    % Plot a figure so the user can check whether normalization was
    % necessary or not:
    if ~params.ScaleChannelsSeparately
        disp(sprintf('Using max=%d, min=%d',MovieMax,MovieMin))
        figure
        subplot(2,1,1)
        % plot(1:numframes,allMaxes,'ob',1:numframes,allMins,'xr')
        % legend('Maxes','Mins')
        plot(1:numframes,allprctiles,'ob',1:numframes,ones(1,numframes).*MovieMax,'--k')
        legend('99.99 Percentile Intensity')
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
        
        figure
        bar(xoutall,allinthist'./numframes)
        title('Histogram of all intensities across the movie','Fontsize',14)
        xlabel('Intensity values','Fontsize',14)
        ylabel('Frequency','Fontsize',14)
        set(gca,'Fontsize',14)
        
    else
        disp(sprintf('Using acceptor max=%d, min=%d',MovieMaxRed,MovieMinRed))
        disp(sprintf('Using donor max=%d, min=%d',MovieMaxGr,MovieMinGr))
        
        figure
        subplot(3,1,1)
        % plot(1:numframes,allMaxes,'ob',1:numframes,allMins,'xr')
        % legend('Maxes','Mins')
        plot(1:numframes,allprctiles,'ob')
        legend('99.99 Percentile Intensity')
        set(gca,'FontSize',14)
        xlabel('Frame','FontSize',14)
        if params.NormImage
            ylabel('Normalized Intensity (a.u.)','FontSize',14)
        else
            ylabel('Raw intensity (a.u.)','FontSize',14)
        end
        subplot(3,1,2)
        % plot(1:numframes,allRedMaxes,'or',1:numframes,allGrMaxes,'og',...
        %     1:numframes,allRedMins,'xr',1:numframes,allGrMins,'xg')
        % legend('Acceptor maxes','Donor maxes','Acceptor mins','Donor mins')
        plot(1:numframes,rprctiles,'or',1:numframes,gprctiles,'og',...
            1:numframes,ones(1,numframes).*MovieMaxRed,'-k',1:numframes,ones(1,numframes).*MovieMaxGr,'--k')
        legend('99.99 Percentile Acceptor Intensities','99.99 Percentile Donor Intensities',...
            'Acceptor Max','Donor Max')
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
        
        figure
        subplot(1,2,1)
        bar(xoutr,rinthist'./numframes)
        title('Histogram of all red intensities across the movie','Fontsize',14)
        xlabel('Intensity values','Fontsize',14)
        ylabel('Frequency','Fontsize',14)
        set(gca,'Fontsize',14)
        subplot(1,2,2)
        bar(xoutg,grinthist'./numframes)
        title('Histogram of all green intensities across the movie','Fontsize',14)
        xlabel('Intensity values','Fontsize',14)
        ylabel('Frequency','Fontsize',14)
        set(gca,'Fontsize',14)
        
    end
    pause
    close all
    
    NormImage = params.NormImage;
    ScaleChannelsSeparately = params.ScaleChannelsSeparately;
    % Save the scaling information to disk for use by LoadScaledMovie:
    save(fullfile(PathToMovie,strcat('ScalingInfo.mat')),'allMedians',...
        'MovieMin','MovieMax','MovieMinRed','MovieMaxRed','MovieMinGr','MovieMaxGr',...
        'NormImage','ScaleChannelsSeparately')
    clear NormImage ScaleChannelsSeparately
end