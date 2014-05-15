% function smFRET
%
% Wrapper function for analyzing smFRET data--loads data, calls the
% functions that do the analysis, etc.
%
% This will analyze all the movies in a set at once, where a set is defined
% as having the same root filename with _1, _2 at the end (the way
% MicroManager saves multiple movies of the same root filename).  They all
% have to be in the same folder.
%
% Inputs are:
% Root filename (the movies it analyzes will be rootname_1, rootname_2 etc)
% Optionally: if you want it to run all the various debugging things, pass
% "1" as the second input.
%
% Steph 9/2013, updated 2/2014 to use a polynomial transformation rather
% than affine to do mapping
% Copyright 2013 Stephanie Johnson, University of California, San Francisco

function smFRET(rootname,debug)

%%%%%%Preliminaries:
    % Since debug is an optional input, set a default
    if ~exist('debug','var') 
        debug=0; 
    else 
        disp('Running in debug mode.')
    end

    % Nested functions for use later:
    
    % Function 1: Getting absolute (in the whole original image) coordinates
    % for spots, rather than local (in a single channel) coordinates.
    function [spotsGglobal,spotsRglobal] = SpotsIntoAbsCoords(spotsGlocal,...
            spotsRlocal,params_struct,imgwidth)
        if params_struct.splitx
                if ~params_struct.Acceptor
                    spotsGglobal(1,:) = spotsGlocal(1,:);
                    spotsGglobal(2,:) = spotsGlocal(2,:)+imgwidth;
                    spotsRglobal = spotsRlocal;
                else
                    spotsRglobal(1,:) = spotsRlocal(1,:);
                    spotsRglobal(2,:) = spotsRlocal(2,:)+imgwidth;
                    spotsGglobal = spotsGlocal;
                end
            else
                if ~params_struct.Acceptor
                    spotsGglobal(1,:) = spotsGlocal(1,:)+imgwidth;
                    spotsGglobal(2,:) = spotsGlocal(2,:);
                    spotsRglobal = spotsRlocal;
                else
                    spotsRglobal(1,:) = spotsRlocal(1,:)+imgwidth;
                    spotsRglobal(2,:) = spotsRlocal(2,:);
                    spotsGglobal = spotsGlocal;
                end
        end
    end

    % Function 2: Let user keep changing background threshhold for
    % spotfinding till satisfied:
    function newspots = SptFindUserThresh(oldspots,SpotImg,thisn,thisxout,...
            ChName,OrdfiltSize,MinDist,Method)
        newspots = oldspots;
        happy = 0;
        while happy==0
            answer = input('Press enter if satisfied, anything else if not:','s');
            if ~isempty(answer)
                close all
                figure,bar(thisxout,thisn)
                title('Choose a threshold between background and true spots:','Fontsize',12)
                xlabel('<-Background intensities ... Real spot intensities->')
                ylabel('Counts')
                newthresh = input('Enter new threshold to use:');
                close
                % Error handling:
                while newthresh <=0 || newthresh >= 1
                    newthresh = input('Enter new threshold to use:');
                end
                [newspots,thisn,thisxout] = FindSpotsV5(SpotImg,'ShowResults',1,...
                    'UserThresh',newthresh,'ImgTitle',ChName,'NeighborhoodSize',OrdfiltSize,...
                    'maxsize',MinDist,'Method',Method);
            else
                happy = 1;
            end
        end
        close
    end

    smFRETsetup;
    params = load('AnalysisParameters.mat');
    % Error-handling: Check that the user-defined parameters are reasonable:
    params.PxlsToExclude = round(params.PxlsToExclude);
    params.EndInjectFrame = round(params.EndInjectFrame);
    MatlabVer = ver;
    MatlabDate = MatlabVer(1).Date;
    if params.UseCombinedImage == 1 && str2double(MatlabDate(end-1:end))<=11 % Testing for Matlab versions older than 2012
        disp('Warning: This version of Matlab does not support creation of a combined image.')
        params.UseCombinedImage = 0;
    end
    clear MatlabVer MatlabDate

%%%%%%FIRST PART: Channel mapping:
    % Load an old channel mapping, or perform a new one:
    DoMap = input('Press enter to perform channel mapping, anything else to load an old one:','s');

    if ~isempty(DoMap)
        % Default to most recent map:
        prevmapdir = load('PathToRecentMap.mat');
        D_Beads = uigetdir(prevmapdir.MostRecentMapDir,'Select directory with old map');
        if exist(fullfile(D_Beads,'ChannelMapping.mat'),'file')
            Map = load(fullfile(D_Beads,'ChannelMapping.mat'));
            tformPoly = Map.tformPoly;
            tformAffine = Map.tformAffine;
            MappingTolerance = Map.MappingTolerance;
            clear Map prevmapdir
        else
            disp(strcat('Bead map not found in',D_Beads))
            return
        end
    else
        D_Beads = uigetdir(params.defaultdatadir,'Select directory with beads');
        % To easily load the most recent bead map:
        MostRecentMapDir = D_Beads;
        % Figure out how many bead files to analyze there are:
        AllBeads = dir(fullfile(D_Beads,'Bead*'));
        num_BeadDir = input(strcat('How many bead files to use for transformation? Max:',...
            int2str(length(AllBeads)),' (Enter to use max)'));
        if isempty(num_BeadDir) || num_BeadDir >= length(AllBeads)
            BdDir = length(AllBeads);
            num_BeadDir = length(AllBeads);
        elseif num_BeadDir < length(AllBeads)
            checkTform = input(strcat('Use remaining ',int2str(length(AllBeads)-num_BeadDir),...
                ' movies to check transform? (y/n)'),'s');
            if strcmpi(checkTform,'y')
                BdDir = length(AllBeads); % But num_BeadDir will still be the 
                    % number of movie files the user wants to include in
                    % the transform
            else
                BdDir = num_BeadDir;
            end
        end

        matchGall = [];
        matchRall = [];
        allBdImgs = [];
        BeadFilesInMap = cell(num_BeadDir,1);

        for i = 1:BdDir
            if i<=num_BeadDir
                BeadFilesInMap{i} = fullfile(D_Beads,AllBeads(i).name); % Keeps a record of which bead files went into the map
            end
            
            TotImg = LoadUManagerTifsV5(fullfile(D_Beads,AllBeads(i).name),[1 params.FramesToAvg]);
            
            % Updated 2/2014 to account for LoadUManagerTifs returning an
            % image of the same numeric type as the initial files.  For
            % spot-finding, we don't care much about numeric type and image
            % contrast scaling, but in order to get the averaged image scaled
            % between 0 and 1, use mean first:
            
            if size(TotImg,3) == 1
                TotImg2 = TotImg;
            else
                % allBdImgs(:,:,i) = mean(TotImg(:,:,1:10),3);
                % With the new version of LoadUManagerTifs:
                TotImg2 = mean(TotImg,3);
            end
            
            allBdImgs(:,:,i) = mat2gray(TotImg2);
            
            disp(sprintf('Spotfinding in movies %d of %d',i,BdDir))

            % Step 1: Find spots in red and green channels separately, so split the
            % image up into the two channels:
            [imgRedRaw,imgGreenRaw] = SplitImg(allBdImgs(:,:,i),params);
            
            % subtract background:
            [~,~,imgRed,imgGreen] = SubBkgnd(imgRedRaw,imgGreenRaw,params);
            % If you don't want to subtract background, uncomment these
            % lines:
            %imgRedsubB = imgRedRaw;
            %imgGreensubB = imgGreenRaw;

            % Find spots in green channel
            if params.Refine_Bd_Cen
                [spotsG{i},n,xout] = FindSpotsV5(imgGreen,'ShowResults',1,'ImgTitle','Green Channel',...
                    'NeighborhoodSize',params.BeadNeighborhood,'maxsize',params.BeadSize,...
                    'Method','GaussFit');
                spotsG{i} = SptFindUserThresh(spotsG{i},imgGreen,n,xout,'Green Channel',...
                    params.BeadNeighborhood,params.BeadSize,'GaussFit');
            else
                [spotsG{i},n,xout] = FindSpotsV5(imgGreen,'ShowResults',1,'ImgTitle','Green Channel',...
                    'NeighborhoodSize',params.BeadNeighborhood,'maxsize',params.BeadSize);
                spotsG{i} = SptFindUserThresh(spotsG{i},imgGreen,n,xout,'Green Channel',...
                    params.BeadNeighborhood,params.BeadSize,'default');
            end
            clear n xout
            
            % Find spots in red channel
            if params.Refine_Bd_Cen
                [spotsR{i},n,xout] = FindSpotsV5(imgRed,'ShowResults',1,'ImgTitle','Red Channel',...
                    'NeighborhoodSize',params.BeadNeighborhood,'maxsize',params.BeadSize,...
                    'Method','GaussFit');
                spotsR{i} = SptFindUserThresh(spotsR{i},imgRed,n,xout,'Red Channel',...
                    params.BeadNeighborhood,params.BeadSize,'GaussFit');
            else
                [spotsR{i},n,xout] = FindSpotsV5(imgRed,'ShowResults',1,'ImgTitle','Red Channel',...
                    'NeighborhoodSize',params.BeadNeighborhood,'maxsize',params.BeadSize);
                spotsR{i} = SptFindUserThresh(spotsR{i},imgRed,n,xout,'Red Channel',...
                    params.BeadNeighborhood,params.BeadSize,'default');
            end
            clear n xout
            
            % spotsG and spotsR are actually params.PxlsToExclude off from
            % the full image's coordinates along the axis along which the
            % channels are split.  So add those pixels back:
            if params.splitx
                spotsG{i}(1,:) = spotsG{i}(1,:);
                spotsG{i}(2,:) = spotsG{i}(2,:)+params.PxlsToExclude;
                spotsR{i}(1,:) = spotsR{i}(1,:);
                spotsR{i}(2,:) = spotsR{i}(2,:)+params.PxlsToExclude;
            else
                spotsG{i}(1,:) = spotsG{i}(1,:)+params.PxlsToExclude;
                spotsG{i}(2,:) = spotsG{i}(2,:);
                spotsR{i}(1,:) = spotsR{i}(1,:)+params.PxlsToExclude;
                spotsR{i}(2,:) = spotsR{i}(2,:);
            end

            if debug %Figure with all the spots found:
                % plot all the boxes for both channels on a big image:
                [spotsG_abs,spotsR_abs] = SpotsIntoAbsCoords(spotsG{i},...
                        spotsR{i},params,size(allBdImgs(:,:,i),2)/2);
                PutBoxesOnImageV4(allBdImgs(:,:,i),[spotsR_abs,spotsG_abs],params.BeadSize,'0','w');
                pause
                close
                clear spotsG_abs spotsR_abs
            end

            % Step 2: figure out which spots in one channel go with the spots in
            % the other channel.  For our beads, which are brighter in the donor
            % than acceptor channel, the matching works best if you match spots
            % in donor to spots in acceptor:
            [matchG{i},matchR{i}] = FindSpotMatches(spotsG{i},spotsR{i});
            if i<=num_BeadDir
                matchGall = [matchGall, matchG{i}];
                matchRall = [matchRall, matchR{i}];
            end

            if debug
                % Box in green the ones that were matched:
                [matchG_abs,matchR_abs] = SpotsIntoAbsCoords(matchG{i},...
                        matchR{i},params,size(allBdImgs(:,:,i),2)/2);
                PutBoxesOnImageV4(allBdImgs(:,:,i),[matchR_abs';matchG_abs'],params.BeadSize);
                title('Only spots that were matched')
                % Another way of plotting the matching: blue line between points
                % in the two channels, with a green dot for where the point is
                % in the green channel, and the end of the blue line where the
                % point it got matched to in the red channel is.  This is also a
                % great way of looking at the distortion between the two
                % channels
                figure('Position',[200 200 325 625])
                plot([matchR{i}(2,:);matchG{i}(2,:)],0-[matchR{i}(1,:);matchG{i}(1,:)],'-b')
                hold on
                plot(matchR{i}(2,:),0-matchR{i}(1,:),'xr')
                if params.splitx == 1
                    ylim([-size(allBdImgs(:,:,i),1) 0])
                    xlim([0 size(allBdImgs(:,:,i),2)/2])
                else
                    ylim([-size(allBdImgs(:,:,i),1)/2 0])
                    xlim([0 size(allBdImgs(:,:,i),2)])
                end
                title('Red x is center of point in red channel')
                pause
                close all
                clear matchG_abs matchR_abs
            end

            clear TotImg TotImg2 imgGreen imgRed spotsGabs
        end

        % Step three: calculate the transformation using all pairs of beads,
        % from all bead movies or snapshots that were loaded.
        % Update 4/2014 to use my FRETmap class
        tformPoly = FRETmap(matchGall,matchRall,'Green',params.TransformToCalc,...
            params.TformMaxDeg,params.TformTotDeg);
        % Affine tends to do better for overlay images using imwarp, so
        % also calculating:
        tformAffine = FRETmap(matchGall,matchRall,'Green','Affine');
        disp(sprintf('Residuals for %d spots:',size(matchGall,2)))
        ResidualsGtoR = tformPoly.ResidualsFwd
        ResidualsRtoG = tformPoly.ResidualsInv
        
        % Update 4/2014: Any mis-pairings of spots has a big effect on the
        % quality of the fitted transform, even for ~750 spots. With
        % perfect matching and our current alignment, the residuals should
        % be <0.008 per spot (note that the residuals will increase with more spots!):
        InitBdNum = size(matchGall,2);
        PrevResid = ResidualsGtoR;
        steps = 0;
        while ResidualsGtoR/size(matchGall,2)>=params.ResidTolerance ||...
                ResidualsRtoG/size(matchGall,2)>=params.ResidTolerance || ...
                abs(ResidualsGtoR-PrevResid)/PrevResid > 0.05 % Stop if residuals stop changing much
            PrevResid = ResidualsGtoR;
            tempRs = tformPoly.FRETmapFwd(matchGall);
            CurrErrors = sqrt((matchRall(1,:)-tempRs(1,:)).^2+(matchRall(2,:)-tempRs(2,:)).^2);
            % I'm not sure why the above two lines aren't identical to:
            %CurrDists = FindSpotDists(tformPoly.FRETmapFwd(matchGall),matchRall);
            %CurrErrors = min(CurrDists,[],2);
            % Using a threshold based only on the G to R transformation,
            % because while the residuals for the two directions do usually
            % differ, it's not usaully by much
            MismatchThresh = mean(CurrErrors)+5*std(CurrErrors);
            clear tempRs
            
            matchGall = [];
            matchRall = [];
            for p = 1:num_BeadDir
                tempG = matchG{p};
                tempR = matchR{p};
                newRs = tformPoly.FRETmapFwd(tempG);
                tempG = tempG(:,sqrt((tempR(1,:)-newRs(1,:)).^2+(tempR(2,:)-newRs(2,:)).^2)<=MismatchThresh);
                tempR = tempR(:,sqrt((tempR(1,:)-newRs(1,:)).^2+(tempR(2,:)-newRs(2,:)).^2)<=MismatchThresh);
                matchG{p} = tempG;
                matchR{p} = tempR;
                matchGall = [matchGall, matchG{p}];
                matchRall = [matchRall, matchR{p}];
                clear newRs tempR tempG
            end
            
            tformPoly = FRETmap(matchGall,matchRall,'Green',params.TransformToCalc,...
                params.TformMaxDeg,params.TformTotDeg);
            disp(sprintf('Residuals for %d spots:',size(matchGall,2)))
            ResidualsGtoR = tformPoly.ResidualsFwd
            ResidualsRtoG = tformPoly.ResidualsInv
            tformAffine = FRETmap(matchGall,matchRall,'Green','Affine');
            if size(matchGall,2)<0.75*InitBdNum || steps>5
                tformPoly.HistResiduals('fwd');
                disp('Channel mapping: having to exclude lots of beads to get residuals down.')
                keyboard
                % If the histogram looks ok, best thing to do is just to
                % manually increase the tolerance:
                % params.ResidTolerance = ResidualsGtoR/size(matchGall,2)+0.001;
                % or something like that
            end
            clear CurrDists CurrErrors
            steps = steps+1;
        end
        clear steps
        
        % Plot the results for each movie:
        for i = 1:num_BeadDir
            disp(strcat('Iterating through bead images for user to check quality (',...
                int2str(i),' of ',int2str(num_BeadDir),')'))
            % Get green points in absolute coordinates
            [matchG_abs,~] = SpotsIntoAbsCoords(matchG{i},...
                matchR{i},params,size(allBdImgs(:,:,i),2)/2);
            % Use tform to map to where they should be in the red channel:
            newR = tformPoly.FRETmapFwd(matchG{i});
            PutBoxesOnImageV4(allBdImgs(:,:,i),[newR';matchG_abs'],params.BeadSize);
            title('Spots found in green, mapped to red','Fontsize',12)
            tformPoly.PlotTform(newR,matchR{i})
            legend('Green spots mapped to red','Red spots')
            ylim([-512 0])
            xlim([0 256])
            tformPoly.TformResiduals(matchG{i},matchR{i},'fwd')
            xlabel('Distance between green spots mapped to red channel, and real red spots','Fontsize',12)
            
            % Show an overlay of one channel on the other:
            [imgRed,imgGreen] = SplitImg(allBdImgs(:,:,i),params);
            tform = tformAffine.ReturnMatlabTform('fwd');
            CombineErr = CalcCombinedImage(tform,imgGreen,imgRed,1);
            if CombineErr ~= -1
                title('Overlay using affine')
            end
            clear tform CombineErr
            tform = tformPoly.ReturnMatlabTform('fwd');
            CombineErr = CalcCombinedImage(tform,imgGreen,imgRed,1);
            if CombineErr ~= -1
                title('Overlay using polynomial')
            end
            clear tform CombineErr
            
            pause
            close all
            clear newR imgRed imgGreen
            
            % Because I explicitly calculated the transformation both ways,
            % check also that the inverse transformation looks ok:
            % Use tform to map to where they should be in the red channel:
            newG = tformPoly.FRETmapInv(matchR{i});
            % Get green points in absolute coordinates
            [newG_abs,~] = SpotsIntoAbsCoords(newG,...
                matchR{i},params,size(allBdImgs(:,:,i),2)/2);
            PutBoxesOnImageV4(allBdImgs(:,:,i),[matchR{i}';newG_abs'],params.BeadSize);
            title('Spots found in red, mapped to green','Fontsize',12)
            tformPoly.PlotTform(newG,matchG{i})
            ylim([-512 0])
            xlim([0 256])
            legend('Green spots mapped to red','Red spots')
            tformPoly.TformResiduals(matchR{i},matchG{i},'inv')
            xlabel('Distance between red spots mapped to green channel, and real green spots','Fontsize',12)
            
            pause
            close all
            clear newG matchG_abs imgRed imgGreen
        end
        
        % The mapping tolerance is the maximal distance away a spot center
        % could be in the other channel, and still possibly be the same as
        % the spot you're looking at:
        tempRs = tformPoly.FRETmapFwd(matchGall);
        tempGs = tformPoly.FRETmapInv(matchRall);
        DistsG = sqrt((matchRall(1,:)-tempRs(1,:)).^2+(matchRall(2,:)-tempRs(2,:)).^2);
        DistsR = sqrt((matchGall(1,:)-tempGs(1,:)).^2+(matchGall(2,:)-tempGs(2,:)).^2);
        MappingTolerance = ceil(max(mean(DistsG)+5*std(DistsG),mean(DistsR)+5*std(DistsR)))
        
        % Lastly, check the transformation, if the user wishes:
        if checkTform
            GforChecking = [];
            RforChecking = [];
            for j = num_BeadDir:length(AllBeads)
                GforChecking = [GforChecking,matchG{j}];
                RforChecking = [RforChecking,matchR{j}];
            end
            tformPoly.TformResiduals(GforChecking,RforChecking,'fwd');
            title('Residuals after applying the calculated transform to remaining bead movies')
            pause
            close
        end
            
        clear allBdImgs matchG matchR matchGall matchRall Dists tempGs tempRs
        clear num_BeadDir BdDir checkTform

        save(fullfile(D_Beads,'ChannelMapping.mat'),'tformPoly','tformAffine',...
            'BeadFilesInMap','MappingTolerance');
        save('PathToRecentMap','MostRecentMapDir');
    end

%%%%%%SECOND PART: Analyze data:
close all
    D_Data = uigetdir(D_Beads,'Select directory with data to analyze');
    % Figure out how many DNA files to analyze there are:
    ToAnalyze = dir(fullfile(D_Data,strcat(rootname,'*')));
    % Get framerate for plotting:
    if isempty(ToAnalyze)
        disp('Did not find data to analyze; remember not to include _<number> at the end of rootname!') %error handling
        return
    end
    fps = GetInfoFromMetaData(fullfile(D_Data,ToAnalyze(1).name),'fps');
    fps = 1/fps; % This is actually frames per ms
    
    % Make sure not saving over old data:
    if ~exist(fullfile(params.defaultsavedir,rootname),'dir')
        savedir = fullfile(params.defaultsavedir,rootname);
        mkdir(savedir)
    else
        saveover = input('Save directory exists; save over? (y/n)','s');
        if ~strcmpi(saveover,'y')
            savepath = uigetdir(params.defaultsavedir,'Choose directory to save data in:');
            newdirname = input('New directory name:','s');
            savedir = fullfile(savepath,newdirname);
            mkdir(savedir)
            clear savepath newdirname
        else
            savedir = fullfile(params.defaultsavedir,rootname);
        end
        clear saveover
    end

    for i = 1:length(ToAnalyze)
        disp(strcat('Analyzing:',ToAnalyze(i).name))
        
        % Reset parameters in case user had loaded an old one, then changed
        % smFRET and wants to use new parameters for this set
        smFRETsetup;

        % Update 12/2013: If this movie has already been analyzed, provide
        % the option to use the previously found spots, instead of
        % re-finding them
        
        useoldspots = 'n';
        
        if exist(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),'file')
            useoldspots = input('Load previously found spots and their intensities? (y/n)','s');
        end
            
        if strcmpi(useoldspots,'n')
           % Finding spots and local background values
           
           % Update 4/2014: Since the movie scaling now takes a while, provide
           % the option to load the previously scaled movie, but a separate
           % option to re-find spots
        
           UseScaledMov = 'n';
        
           if exist(fullfile(D_Data,ToAnalyze(i).name,'ScaledMovieFrames1to100.mat'),'file')
               UseScaledMov = input('Load scaled movie? (y/n)','s');
           end
           
           if strcmpi(UseScaledMov,'n')
               % Update 4/2014: Scaling the movie first, before finding spots,
               % so that the background value calculated when the spot centers
               % are refined by a GaussFit are meaningful:
               alltifs = dir(fullfile(D_Data,ToAnalyze(i).name,'img*.tif'));
               disp('Scaling movie ...')
               ScaleMovieV2(fullfile(D_Data,ToAnalyze(i).name),length(alltifs),params);
           end
            
           % Load the first FramesToAvg frames for finding spots
    %            % Using the old version of ScaleMovie:
    %            TotImg = LoadUManagerTifsV5(fullfile(D_Data,ToAnalyze(i).name),[1 params.FramesToAvg]);
    %            %TotImg = LoadUManagerTifsV5(fullfile(D_Data,ToAnalyze(i).name),[2000-params.FramesToAvg 2000]);
    %            
    %            if size(TotImg,3) > 1
    %                TotImg = mean(TotImg,3);
    %            end
    %             
    %            TotImg = mat2gray(TotImg);
    %            
    %            [imgRed,imgGreen] = SplitImg(TotImg,params);
    
           % Update 5/2014: Added a parameter in smFRETsetup that allows
           % the user to choose where to start spotfinding (in case, for
           % example, manual flowthroughs cause distortions to early parts
           % of the movie)
           %[imgRed,imgGreen] = LoadScaledMovie(fullfile(D_Data,ToAnalyze(i).name),[1 1+params.FramesToAvg]);
           [imgRed,imgGreen] = LoadScaledMovie(fullfile(D_Data,ToAnalyze(i).name),...
               [params.EndInjectFrame params.EndInjectFrame+params.FramesToAvg]);
           imgRedavg = mat2gray(mean(imgRed,3)); %Do I want to do mat2gray here? Update 4/2014:
                % since I'm going to treat spotfinding as totally separate
                % from Gauss fitting for intensity smoothing, it is best
                % that I do mat2gray here
           imgGreenavg = mat2gray(mean(imgGreen,3));
           
           % Step 0: subtract background:
           [imgRbkgnd,imgGbkgnd,imgRMinusBkgnd,imgGMinusBkgnd] = SubBkgnd(imgRedavg,imgGreenavg,params);
           % If you don't want to subtract background, uncomment these
            % lines:
           %imgRMinusBkgnd = imgRedavg;
           %imgGMinusBkgnd = imgGreenavg;
           
           % Step 1: find spots
           % Find spots in both channels, but don't double-count. Allow
           % user to decide whether to find spots separately in each
           % channel, or using a combined image (see notes on the
           % UseCombinedImage parameter in smFRETsetup.m).
           if params.UseCombinedImage
               disp('smFRET: CombinedImage option is not currently functional!')
               return
               
               % Finding spots in a composite image, so that
               % mid-FRET spots don't get lost.  NOTE that the combined image
               % will have a frame of reference of the acceptor image, which
               % is fine because that's what I pass into UserSpotSelectionV4.
               
               % Step 1.1 Create a combined image
               composite = CalcCombinedImage(tformGtoRAffine,imgGMinusBkgnd,imgRMinusBkgnd);
               % The built-in Matlab function imfuse used to create the output
               % of CalcCombinedImage only returns unit8 images, but fminsearch
               % (called in Fit2DGaussToSpot in GetGaussParams below) needs a 
               % double, so convert back to doubles:
               composite = mat2gray(composite);

               % Step 1.2 Identify spots in this combined image:
               [spotsR,n,xout] = FindSpotsV5(composite,'ShowResults',1,'ImgTitle','Composite Image',...
                     'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize);
               spots = SptFindUserThresh(spotsR,composite,n,xout,'Composite Image',...
                   params.DNANeighborhood,params.DNASize,'default');
               clear n xout
               
               %if params.IntensityGaussWeight
                    disp('Refining spot centers by 2D Gauss fit')
                    [RefinedCenters,Vars,bkgnd] = RefineCensByGauss(spots,composite,params.DNASize,0);
%                else
%                    RefinedCenters = spots;
%                    Vars = -1;
%                    bkgnd = -1;
%                end

               close all
           else
               % Step 1.1 Identify spots in acceptor channel, and refine centers by
               % fitting to a Gaussian, regardless of whether or not user
               % wants to weight intensities by a Gaussian:
               [RefinedCentersR,n,xout] = FindSpotsV5(imgRMinusBkgnd,'ShowResults',1,'ImgTitle','Red Channel',...
                     'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize,...
                     'Method','GaussFit');
               RefinedCentersR = SptFindUserThresh(RefinedCentersR,imgRMinusBkgnd,n,xout,'Red Channel',...
                   params.DNANeighborhood,params.DNASize,'GaussFit');
               clear n xout

               close
               
               % And in donor channel
               [RefinedCentersG,n,xout] = FindSpotsV5(imgGMinusBkgnd,'ShowResults',1,'ImgTitle','Green Channel',...
                     'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize,...
                     'Method','GaussFit');
               RefinedCentersG = SptFindUserThresh(RefinedCentersG,imgGMinusBkgnd,n,xout,'Green Channel',...
                   params.DNANeighborhood,params.DNASize,'GaussFit');
               clear n xout

               close
               
               % Step 1.2: If the user wants to refine spot intensities by
               % a Gaussian, fit variances for those Gaussian weights. Note
               % that this is done on the scaled but not background
               % subtracted, and not averaged, movie
               if params.IntensityGaussWeight
                   disp('Fitting Gaussians to find spot variances')
                   VarsR = FindSpotVars(imgRed,RefinedCentersR,params);
                   VarsG = FindSpotVars(imgGreen,RefinedCentersG,params);
               else
                   VarsR = -1;
                   VarsG = -1;
               end
               
               clear imgRed imgGreen imgRedavg imgGreenavg
               clear imgRbkgnd imgGbkgnd imgRMinusBkgnd imgGMinusBkgnd
               
               % Step 1.3: Check that no spots are double counted, and get
               % rid of spots that are too close together.
               % First convert green channel spot locations to red
               % channel locations, and check that no spots are double
               % counted. Based on the MappingTolerance parameter saved 
               % during the mapping procedure, we know spots that are truly 
               % the same will have their centers at most MappingTolerance
               % apart.
               spotsGinR = tformPoly.FRETmapFwd(RefinedCentersG);
               Dists = FindSpotDists(RefinedCentersR,spotsGinR);
               spotnottooclose = Dists>MappingTolerance;
               % As in FindSpotsV5, each column of spotnottooclose will be all 1's if the 
               % spotGinR represented by the column is more than MappingTolerance away from 
               % the spotR represented by each row. If there's already a spotR too
               % close, one or more elements of the column will be zero.  So now
               % ask if the sum of each column is equal to the length of the
               % column. If so, add it to spots. Note that FindSpots has
               % already checked that no spot in either channel is too
               % close to another spot in that same channel (where "too close" is
               % defined as params.DNASize); but it's possible that there's a high-FRET
               % spot only in the red channel, and a low-FRET spot only in
               % the green channel, which are truly different DNAs but are
               % too close to get FRET out of.  So far I don't deal with
               % that possibility.
               spots = spotsGinR(:,sum(spotnottooclose,1)==size(spotnottooclose,1));
               % Before adding the unique spots found in the red channel:
               if params.IntensityGaussWeight
                   % Problems I haven't really resolved yet:
                    % (1) Assume variance is same in green and red channels?
                    % This is what the Ha lab IDL code does (actually they
                    % hard-code a value for both channels, for all spots)
                    % (2) If spot is found in both channels (e.g. mid-FRET), 
                    % which channel to use for vars? Right now I'm
                    % automatically using red channel. FindSpotVars can
                    % also output the amplitudes of the fits, so I could
                    % also choose based on which channel had a larger
                    % amplitude or something.
                   Vars = VarsG(:,sum(spotnottooclose,1)==size(spotnottooclose,1));
                   Vars(:,end+1:end+size(VarsR,2)) = VarsR;
               else
                   Vars = -1;
               end
               spots(:,end+1:end+size(RefinedCentersR,2)) = RefinedCentersR;
           end
           
           disp(sprintf('Found %d total spots',size(spots,2)))
           
           % Step 2: Load the whole movie in increments and calculate the
           % intensity of each spot in each frame.
           
           disp('Calculating frame-by-frame intensities')
           
           [RedI, GrI] = CalcIntensitiesV2(fullfile(D_Data,ToAnalyze(i).name),...
               spots, Vars, tformPoly,params);
           
           % Save spot positions, intensities and associated GaussFit
           % parameters in case the user wants to re-analyze.
           % Note because spots are never displayed on the full image,
           % only on the image split into two channels, I don't need to
           % add params.PxlsToExclude to get spots into the right
           % coordinates (unlike with the beads)
           
           SpotsInR = spots;
           SpotVars = Vars;
           save(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),'SpotsInR',...
               'SpotVars','RedI','GrI')
           clear SpotsInR SpotVars
           
           if i==1
               % Also save the params structure in the data analysis folder, so
               % you know what analysis parameters were used to analyze the
               % data:
               params.fps = fps;
               save(fullfile(savedir,strcat('AnalysisParameters.mat')),'params');
           end
           
           % Step 3: Display a trace of intensity-vs-time for each spot,
           % with an interactive section for the user to select spots that
           % are true FRET, etc

           disp(strcat('Movie ',int2str(i),'of',int2str(length(ToAnalyze))))
           UserSpotSelectionV4(RedI,GrI,spots,...
               fullfile(D_Data,ToAnalyze(i).name),params,tformPoly,savedir,fps,i);
        else %If the user wants to instead use previously saved data
           oldspots = load(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')));
           params = load(fullfile(savedir,strcat('AnalysisParameters.mat')));
           params = params.params;
           disp(strcat('Movie ',int2str(i),'of',int2str(length(ToAnalyze))))
           UserSpotSelectionV4(oldspots.RedI,oldspots.GrI,oldspots.SpotsInR,...
               fullfile(D_Data,ToAnalyze(i).name),params,tformPoly,savedir,fps,i);
        end
        clear TotImg spots imgRed imgGreen spotsG spotsR spotsG_abs spotsRguess spotstemp
    end
end
