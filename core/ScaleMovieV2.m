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
                    
                    PlotMaxes(xaxislim,[],[],Maxes,newMax,Mins,[],NormImage)
                    keyboard
                    close
                end
            end
            clear sortedMaxes MaxDiffs meanDiff stdDiff
        end
    
        % Plotting function
        function PlotMaxes(n_frames,medians,allprctile,R_prctile,R_max_val,G_prctile,G_max_val,norm)
            figure
            if isempty(medians)
                plot(1:n_frames,R_prctile,'ob',1:n_frames,G_prctile,'xr',...
                    1:n_frames,ones(1,n_frames).*R_max_val,'--k')
                legend('99.99 percentiles','Mins','Max to scale to')
            else
                if ~isempty(G_max_val)
                    subplot(3,1,1)
                    plot(1:n_frames,allprctile,'ob')
                else
                    subplot(2,1,1)
                    plot(1:n_frames,R_prctile,'ob',...
                        1:n_frames,ones(1,n_frames).*R_max_val,'--k')
                end
                
                legend('99.99 Percentile Intensity')
            end
            set(gca,'FontSize',14)
            xlabel('Frame','FontSize',14)
            if norm
                ylabel('Normalized Intensity (a.u.)','FontSize',14)
            else
                ylabel('Raw intensity (a.u.)','FontSize',14)
            end
            
            if ~isempty(medians)
                if ~isempty(G_max_val)
                    subplot(3,1,2)
                    plot(1:n_frames,R_prctile,'or',1:n_frames,G_prctile,'og',...
                        1:n_frames,ones(1,n_frames).*R_max_val,'-k',1:n_frames,ones(1,n_frames).*G_max_val,'--k')
                    legend('99.99 Percentile Acceptor Intensities','99.99 Percentile Donor Intensities',...
                        'Acceptor Max','Donor Max')
                end
                set(gca,'FontSize',14)
                xlabel('Frame','FontSize',14)
                if norm
                    ylabel('Normalized Intensity (a.u.)','FontSize',14)
                else
                    ylabel('Raw intensity (a.u.)','FontSize',14)
                end
                
                if ~isempty(G_max_val)
                    subplot(3,1,3)
                else
                    subplot(2,1,2)
                end
                plot(1:n_frames,medians,'.g')
                legend('Median')
                set(gca,'FontSize',14)
                xlabel('Frame','FontSize',14)
                ylabel('Raw Intensity (a.u.)','FontSize',14)
            end
        end
    
        % Helper functions to compute stuff
        function [out_max,out_min] = ComputeMaxes(in_max,in_min,data)
            out_max = max(in_max,double(max(data(:))));
            out_min = min(in_min,double(min(data(:))));
        end
        function [out_maxes,out_mins] = ComputeDoubleMax(data)
            out_maxes = max(max(data,[],2),[],1);
        	out_mins = min(min(data,[],2),[],1);
        end
        function prctiles = ComputePercentiles(data)
            rdata = reshape(data,size(data,1)*size(data,2),size(data,3));
            prctiles = prctile(rdata,99.99);
        end
        function rdata = ReshapeData(data,meds,index)
            rdata = bsxfun(@rdivide,data,...
                reshape(meds(index:index+size(data,3)-1),1,1,size(data,3)));
        end

    % Get some preliminary info:
    framesize = GetInfoFromMetaData(PathToMovie,'imgsize');
    NormImage = params.NormImage;
    ScaleChannelsSeparately = params.ScaleChannelsSeparately;
 
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
        if NormImage
            moviebit = ReshapeData(moviebit,allMedians,jj);
            redbit = ReshapeData(redbit,allMedians,jj);
            grbit = ReshapeData(grbit,allMedians,jj);
        end

        [MovieMax,MovieMin] = ComputeMaxes(MovieMax,MovieMin,moviebit);
        [MovieMaxGr,MovieMinGr] = ComputeMaxes(MovieMaxGr,MovieMinGr,grbit);
        [MovieMaxRed,MovieMinRed] = ComputeMaxes(MovieMaxRed,MovieMinRed,redbit);

        [allMaxes(jj:jj+size(moviebit,3)-1),...
            allMins(jj:jj+size(moviebit,3)-1)] = ComputeDoubleMax(moviebit);
        [allGrMaxes(jj:jj+size(moviebit,3)-1),...
            allGrMins(jj:jj+size(moviebit,3)-1)] = ComputeDoubleMax(grbit);
        [allRedMaxes(jj:jj+size(moviebit,3)-1),...
            allRedMins(jj:jj+size(moviebit,3)-1)] = ComputeDoubleMax(redbit);
                
        % Update 1/2015: Reshaping red and green to calculate prctile:
        allprctiles(jj:jj+size(moviebit,3)-1) = ComputePercentiles(moviebit);
        rprctiles(jj:jj+size(moviebit,3)-1) = ComputePercentiles(redbit);
        gprctiles(jj:jj+size(moviebit,3)-1) = ComputePercentiles(grbit);
        
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
        InjectPoints = DetectRedFlash(rprctiles,gprctiles,allMedians,params);
        % These artificial maxes (from flashing the red laser) mess with
        % the calculation of the true max to which to scale the image. So
        % removing them here:
        for ss = 1:length(InjectPoints)
            allprctiles(InjectPoints(ss)*(params.fps)-3:InjectPoints(ss)*(params.fps)+3) = median(allprctiles).*ones(7,1);
            rprctiles(InjectPoints(ss)*(params.fps)-3:InjectPoints(ss)*(params.fps)+3) = median(rprctiles).*ones(7,1);
            allMedians(InjectPoints(ss)*(params.fps)-3:InjectPoints(ss)*(params.fps)+3) = median(allMedians).*ones(7,1);
        end
    end
    
    % Check that the calculated Max(es) isn't/aren't way larger than the rest of the
    % maxes (it happens): Note it seems to work better to take the maxes instead of 
    % a percentile of thepercentiles.
    MovieMax = max(allprctiles);
    MovieMaxGr = max(gprctiles);
    MovieMaxRed = max(rprctiles);
    disp('Checking for outliers in overall percentiles ... ')
    MovieMax = CheckMaxes(allprctiles,allMins,MovieMax,numframes,NormImage);
    if ScaleChannelsSeparately
        disp('Checking for outliers in donor percentiles ... ')
        MovieMaxGr = CheckMaxes(gprctiles,allGrMins,MovieMaxGr,numframes,NormImage);
        disp('Checking for outliers in acceptor percentiles ... ')
        MovieMaxRed = CheckMaxes(rprctiles,allRedMins,MovieMaxRed,numframes,NormImage);
    end
    
    % Plot a figure so the user can check whether normalization was
    % necessary or not:
    if ~ScaleChannelsSeparately
        disp(sprintf('Using max=%d, min=%d',MovieMax,MovieMin))
        
        PlotMaxes(numframes,allMedians,[],allprctiles,MovieMax,[],[],NormImage)
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
        
        PlotMaxes(numframes,allMedians,allprctiles,rprctiles,MovieMaxRed,gprctiles,MovieMaxGr,NormImage)
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
   
    % Save the scaling information to disk for use by LoadScaledMovie:
    if exist('InjectPoints','var')
        save(fullfile(PathToMovie,strcat('ScalingInfo.mat')),'allMedians',...
            'MovieMin','MovieMax','MovieMinRed','MovieMaxRed','MovieMinGr','MovieMaxGr',...
            'NormImage','ScaleChannelsSeparately','InjectPoints')
    else
        save(fullfile(PathToMovie,strcat('ScalingInfo.mat')),'allMedians',...
            'MovieMin','MovieMax','MovieMinRed','MovieMaxRed','MovieMinGr','MovieMaxGr',...
            'NormImage','ScaleChannelsSeparately')
    end
end