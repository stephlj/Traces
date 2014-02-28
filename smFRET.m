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
% Steph 9/2013
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
    % Also makes use of the new PxlsToExclude parameter.
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
    if round(params.PxlsToExclude)~=params.PxlsToExclude
        params.PxlsToExclude = round(params.PxlsToExclude);
    end

%%%%%%FIRST PART: Channel mapping:
    % Load an old channel mapping, or perform a new one:
    DoMap = input('Press enter to perform channel mapping, anything else to load an old one:','s');

    if ~isempty(DoMap)
        % Default to most recent map:
        prevmapdir = load('PathToRecentMap.mat');
        D_Beads = uigetdir(prevmapdir.MostRecentMapDir,'Select directory with old map');
        if exist(fullfile(D_Beads,'ChannelMapping.mat'),'file')
            Map = load(fullfile(D_Beads,'ChannelMapping.mat'));
            A = Map.A;
            b = Map.b;
            Amatlab = Map.Amatlab;
            bmatlab = Map.bmatlab;
            clear Map prevmapdir
        else
            disp(strcat('Bead map not found in',Beads))
            return
        end
    else
        D_Beads = uigetdir(params.defaultdatadir,'Select directory with beads');
        % To easily load the most recent bead map:
        MostRecentMapDir = D_Beads;
        % Figure out how many bead files to analyze there are:
        AllBeads = dir(fullfile(D_Beads,'Bead*'));
        num_BeadDir = input(strcat('How many bead files to analyze? Max:',...
            int2str(length(AllBeads)),' (Enter to use max; movies load first)'));
        if isempty(num_BeadDir)
            BdDir = length(AllBeads);
        else
            BdDir = num_BeadDir;
        end

        matchGall = [];
        matchRall = [];
        allBdImgs = [];
        BeadFilesInMap = cell(BdDir,1);

        for i = 1:BdDir
            BeadFilesInMap{i} = fullfile(D_Beads,AllBeads(i).name); % Keeps a record of which bead files went into the map
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
            
            disp(strcat('Analyzing beads: ',int2str(i)',' of ',int2str(BdDir)))

            % Step 1: Find spots in red and green channels separately, so split the
            % image up into the two channels:
            [imgRed,imgGreen] = SplitImg(allBdImgs(:,:,i),params);

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
            matchGall = [matchGall, matchG{i}];
            matchRall = [matchRall, matchR{i}];

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
        % from all bead movies or snapshots that were loaded:
        % Calculate transformation:
        % Update 1/2014: the built-in Matlab function fitgeotrans does a bit
        % better than my hand-written code in CalcChannelMapping for overlaying 
        % images using CalcCombinedImage, but mine does better in terms of 
        % calculating where a spot in the donor channel should be in the acceptor
        % channel.  So for now, calculating and saving both:
        [A,b] = CalcChannelMapping(matchGall,matchRall)
        tform = fitgeotrans(matchRall',matchGall','Affine'); % Note different input order for fitgeotrans
        Amatlab = tform.T(1:2,1:2)
        bmatlab = transpose(-tform.T(3,1:2))
        
        % Plot the results for each movie:
        for i = 1:BdDir
            disp(strcat('Iterating through bead images for user to check quality (',...
                int2str(i),' of ',int2str(BdDir),')'))
            [matchG_abs,~] = SpotsIntoAbsCoords(matchG{i},...
                matchR{i},params,size(allBdImgs(:,:,i),2)/2);
            newR = CalcSpotTransform(matchG{i},[],A,b);
            PutBoxesOnImageV4(allBdImgs(:,:,i),[newR';matchG_abs'],params.BeadSize);
            title('Spots found in green, mapped to red','Fontsize',12)
            figure
            errs = FindSpotDists(matchR{i},newR);
            hist(min(errs,[],2),0:0.1:10)
            hold on
            plot([mean(min(errs,[],2)) mean(min(errs,[],2))], [0 size(errs,1)/4],'--k');
            hold off
            % TODO: If the mean error is bigger than, say, 1 pxl, redo the
            % mapping excluding points with too-large errors
            ylabel('Counts','Fontsize',12)
            xlabel('Distance between mapped red bead and real red bead','Fontsize',12)
            if params.splitx
                CalcCombinedImage(Amatlab,bmatlab,...
                    allBdImgs(:,(size(allBdImgs(:,:,i))/2)+1:end,i),allBdImgs(:,1:size(allBdImgs(:,:,i))/2,i),1);
            else
                CalcCombinedImage(Amatlab,bmatlab,...
                    allBdImgs((size(allBdImgs(:,:,i))/2)+1:end,:,i),allBdImgs(1:size(allBdImgs(:,:,i))/2,:,i),1);
            end 
            pause
            close
            close
            close
            clear newR matchG_abs
        end
        clear allBdImgs matchG matchR matchGall matchRall

        save(fullfile(D_Beads,'ChannelMapping.mat'),'A','b','BeadFilesInMap',...
            'Amatlab','bmatlab');
        save('PathToRecentMap','MostRecentMapDir');
    end

%%%%%%SECOND PART: Analyze data:
close all
    D_Data = uigetdir(D_Beads,'Select directory with data to analyze');
    % Figure out how many bead files to analyze there are:
    ToAnalyze = dir(fullfile(D_Data,strcat(rootname,'_*')));
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
        % Update 12/2013: If this movie has already been analyzed, provide
        % the option to use the previously found spots, instead of
        % re-finding them
        
        useoldspots = 'n';
        
        if exist(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),'file')
            useoldspots = input('Load previous analysis? (y/n)','s');
        end
            
        if strcmpi(useoldspots,'n')
           % Finding spots and local background values
            
           % Load the first FramesToAvg frames for finding spots
           TotImg = LoadUManagerTifsV5(fullfile(D_Data,ToAnalyze(i).name),[1 params.FramesToAvg]);
           
           if size(TotImg,3) > 1
               TotImg = mean(TotImg,3);
           end
            
           TotImg = mat2gray(TotImg);
           
           [imgRed,imgGreen] = SplitImg(TotImg,params);
           
           % Step 0: TODO: calculate local background values 
           
           % Find spots in both channels, but don't double-count:
           % Update 1/2014: finding spots in a composite image, so that
           % mid-FRET spots don't get lost.  NOTE that the combined image
           % will have a frame of reference of the acceptor image, which
           % is fine because that's what I pass into UserSpotSelectionV4.
           
           composite = CalcCombinedImage(Amatlab,bmatlab,imgGreen,imgRed);
           % The built-in Matlab function imfuse used to create the output
           % of CalcCombinedImage only returns unit8 images, but fminsearch
           % (called in Fit2DGaussToSpot below) needs a double, so convert
           % back to doubles:
           composite = mat2gray(composite);
           
           % Step 1: identify spots in this combined image:
           [spotsR,n,xout] = FindSpotsV5(composite,'ShowResults',1,'ImgTitle','Composite Image',...
                 'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize);
           spots = SptFindUserThresh(spotsR,composite,n,xout,'Composite Image',...
               params.DNANeighborhood,params.DNASize,'default');
           clear n xout
           
           close all
           
           % Step 2: fit a Gaussian to the spots found in the combined
           % image to get values that will be used to refine the
           % intensity-versus-time calculation later:
           
           [RefinedCenters,Vars] = GetGaussParams(spotsR,composite,imgGreen,...
               imgRed,A,b, params.DNASize);
           
           % How different are these fit values from the fit parameters for
           % each spot in its own channel?
           
           % Step 3: Load the whole movie in increments and calculate the
           % intensity of each spot in each frame.
           
           % Save spot positions, intensities and associated GaussFit
           % parameters in case the user wants to re-analyze.
           % Note because spots are never displayed on the full image,
           % only on the image split into two channels, I don't need to
           % add params.PxlsToExclude to get spots into the right
           % coordinates (unlike with the beads)
            
           save(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),'spots')
           
           % Step 4: Display a trace of intensity-vs-time for each spot,
           % with an interactive section for the user to select spots that
           % are true FRET, etc

           disp(strcat('Movie ',int2str(i),'of',int2str(length(ToAnalyze))))
           UserSpotSelectionV4(spots,fullfile(D_Data,ToAnalyze(i).name),params,A,b,savedir,fps,i);
        else %If the user wants to instead use previously saved data
           oldspots = load(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),'spots');
           disp(strcat('Movie ',int2str(i),'of',int2str(length(ToAnalyze))))
           UserSpotSelectionV4(oldspots.spots,fullfile(D_Data,ToAnalyze(i).name),params,A,b,savedir,fps,i);
        end
        clear TotImg spots imgRed imgGreen spotsG spotsR spotsG_abs spotsRguess spotstemp
    end
end
