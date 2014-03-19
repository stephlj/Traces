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
            tformGtoR = Map.tformGtoR;
            tformRtoG = Map.tformRtoG;
            tformGtoRAffine = Map.tformGtoRAffine;
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
            int2str(length(AllBeads)),' (Enter to use max)'));
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
        % from all bead movies or snapshots that were loaded.
        % Update 2/2014: An affine transformation (using my CalcChannelMapping),
        % does ok at channel mapping--the vast majority of spots will be
        % mapped to within 1.5 pixels of their true center. Matlab's
        % fitgeotrans with the affine option doesn't do well at mapping
        % points, but does better at getting overlaid channels.
        % However, in order to fit Gaussians to each spot in real data
        % (below), the mapping really needs to be good to within 0.5 pixel
        % or less.  For that we need a 3rd-degree polynomial, not an
        % affine.  So using Matlab's cp2tform:
        % tform = cp2tform(matchRall',matchGall','polynomial',3); % Note different input order than my CalcChannelMapping
        % But that does really horribly, and is really annoying, for making
        % overlaid/composite images, which I need for the second half of
        % this code. So also doing an affine transformation using
        % fitgeotrans.
        % Hey wait a sec, even though polynomial isn't listed as an option
        % for fitgeotrans, it can do it!  sheesh
        tformGtoR = fitgeotrans(matchRall',matchGall','polynomial',3);
        % Unlike an affine transformation, I believe this is not
        % invertible, so need to also calculate the other direction:
        tformRtoG = fitgeotrans(matchGall',matchRall','polynomial',3);
        % Ok actually the polynomial transformation does a terrible job of
        % overlaying the two channels to make a combined image. So, lastly,
        % also using fitgeotrans to get an affine transformation:
        tformGtoRAffine = fitgeotrans(matchRall',matchGall','Affine');
        
        % Plot the results for each movie:
        for i = 1:BdDir
            disp(strcat('Iterating through bead images for user to check quality (',...
                int2str(i),' of ',int2str(BdDir),')'))
            % Get green points in absolute coordinates
            [matchG_abs,~] = SpotsIntoAbsCoords(matchG{i},...
                matchR{i},params,size(allBdImgs(:,:,i),2)/2);
            % Use tform to map to where they should be in the red channel:
            newR = transpose(transformPointsInverse(tformGtoR,matchG{i}'));
            PutBoxesOnImageV4(allBdImgs(:,:,i),[newR';matchG_abs'],params.BeadSize);
            title('Spots found in green, mapped to red','Fontsize',12)
            figure
            errs = FindSpotDists(matchR{i},newR);
            hist(min(errs,[],2),[0:0.1:10])
            hold on
            plot([mean(min(errs,[],2)) mean(min(errs,[],2))], [0 size(errs,1)/4],'--k');
            hold off
            % TODO: If the mean error is bigger than, say, 1 pxl, redo the
            % mapping excluding points with too-large errors. 
            % Update 2/2014: it's possible fitgeotrans does this already?
            % The older version, cp2tform, does, I believe
            ylabel('Counts','Fontsize',12)
            xlabel('Distance between mapped red bead and real red bead','Fontsize',12)
            
            % Show an overlay of one channel on the other:
            [imgRed,imgGreen] = SplitImg(allBdImgs(:,:,i),params);
            CalcCombinedImage(tformGtoRAffine,imgGreen,imgRed,1);
            title('Overlay using affine')
            CalcCombinedImage(tformGtoR,imgGreen,imgRed,1);
            title('Overlay using polynomial')
            
            pause
            close
            close
            close
            close
            clear newR errs
            
            % Because I explicitly calculated the transformation both ways,
            % check also that the inverse transformation looks ok:
            % Use tform to map to where they should be in the red channel:
            newG = transpose(transformPointsInverse(tformRtoG,matchR{i}'));
            % Get green points in absolute coordinates
            [newG_abs,~] = SpotsIntoAbsCoords(newG,...
                matchR{i},params,size(allBdImgs(:,:,i),2)/2);
            PutBoxesOnImageV4(allBdImgs(:,:,i),[matchR{i}';newG_abs'],params.BeadSize);
            title('Spots found in red, mapped to green','Fontsize',12)
            figure
            errs = FindSpotDists(matchG{i},newG);
            hist(min(errs,[],2),[0:0.1:10])
            hold on
            plot([mean(min(errs,[],2)) mean(min(errs,[],2))], [0 size(errs,1)/4],'--k');
            hold off
            ylabel('Counts','Fontsize',12)
            xlabel('Distance between mapped green bead and real green bead','Fontsize',12)
            
            pause
            close
            close
            clear newG matchG_abs imgRed imgGreen errs
        end
        clear allBdImgs matchG matchR matchGall matchRall

        save(fullfile(D_Beads,'ChannelMapping.mat'),'tformGtoR','tformRtoG',...
            'tformGtoRAffine','BeadFilesInMap');
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
           
           % Step 0: subtract background:
           [imgRMinusBkgnd,imgGMinusBkgnd] = SubBkgnd(imgRed,imgGreen,params);
           
           % Find spots in both channels, but don't double-count:
           % Update 1/2014: finding spots in a composite image, so that
           % mid-FRET spots don't get lost.  NOTE that the combined image
           % will have a frame of reference of the acceptor image, which
           % is fine because that's what I pass into UserSpotSelectionV4.
           
           composite = CalcCombinedImage(tformGtoRAffine,imgGMinusBkgnd,imgRMinusBkgnd);
           % The built-in Matlab function imfuse used to create the output
           % of CalcCombinedImage only returns unit8 images, but fminsearch
           % (called in Fit2DGaussToSpot in GetGaussParams below) needs a 
           % double, so convert back to doubles:
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
           
           disp('Refining spot centers by 2D Gauss fit')
           
           [RefinedCenters,Vars] = GetGaussParams(spots,composite,imgGMinusBkgnd,...
               imgRMinusBkgnd,tformRtoG,tformGtoR,params.DNASize);
           
           % Some notes about this fitting process:
           % (1) If you do this with beads instead of DNA, in which case
           % the "right answers" are obvious, the green channel is always
           % the best fit (which makes sense since the beads are brightest
           % in green), and it always keeps all the beads as "good"
           % (2) As long as you're consistent, it doesn't seem to matter
           % whether you pass Amatlab, bmatlab or A, b into GetGaussParams.
           % However, whichever parameter set you pass is the one you need
           % to use to convert between channels from now on! It also
           % doesn't seem to make a huge difference if you use A, b even if
           % you use Amatlab,bmatlab to make composite. Again just be
           % consistent from this point onwards.
           
           % Step 3: Load the whole movie in increments and calculate the
           % intensity of each spot in each frame.
           
           disp('Calculating frame-by-frame intensities')
           
           [RedI, GrI] = CalcIntensitiesV2(fullfile(D_Data,ToAnalyze(i).name),...
               RefinedCenters, Vars, tformRtoG,params);
           
           % Save spot positions, intensities and associated GaussFit
           % parameters in case the user wants to re-analyze.
           % Note because spots are never displayed on the full image,
           % only on the image split into two channels, I don't need to
           % add params.PxlsToExclude to get spots into the right
           % coordinates (unlike with the beads)
           
           SpotsInR = RefinedCenters;
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
           
           % Step 4: Display a trace of intensity-vs-time for each spot,
           % with an interactive section for the user to select spots that
           % are true FRET, etc

           disp(strcat('Movie ',int2str(i),'of',int2str(length(ToAnalyze))))
           UserSpotSelectionV4(RedI,GrI,RefinedCenters,...
               fullfile(D_Data,ToAnalyze(i).name),params,tformRtoG,savedir,fps,i);
        else %If the user wants to instead use previously saved data
           oldspots = load(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')));
           params = load(fullfile(savedir,strcat('AnalysisParameters.mat')));
           params = params.params;
           disp(strcat('Movie ',int2str(i),'of',int2str(length(ToAnalyze))))
           UserSpotSelectionV4(oldspots.RedI,oldspots.GrI,oldspots.SpotsInR,...
               fullfile(D_Data,ToAnalyze(i).name),params,tformRtoG,savedir,fps,i);
        end
        clear TotImg spots imgRed imgGreen spotsG spotsR spotsG_abs spotsRguess spotstemp
    end
end
