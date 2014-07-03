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

    % Subfunctions for use later:
    
    % Function 1: Getting absolute (in the whole original image) coordinates
    % for spots, rather than local (in a single channel) coordinates.
    function [spotsGglobal,spotsRglobal] = SpotsIntoAbsCoords(spotsGlocal,...
            spotsRlocal,params_struct,imgwidth)
        % spotsG and spotsR are actually params.PxlsToExclude off from
        % the full image's coordinates along the axis along which the
        % channels are split.  So add those pixels back:
        if params_struct.splitx
            spotsGlocal(1,:) = spotsGlocal(1,:);
            spotsGlocal(2,:) = spotsGlocal(2,:)+params_struct.PxlsToExclude;
            spotsRlocal(1,:) = spotsRlocal(1,:);
            spotsRlocal(2,:) = spotsRlocal(2,:)+params_struct.PxlsToExclude;
        else
            spotsGlocal(1,:) = spotsGlocal(1,:)+params_struct.PxlsToExclude;
            spotsGlocal(2,:) = spotsGlocal(2,:);
            spotsRlocal(1,:) = spotsRlocal(1,:)+params_struct.PxlsToExclude;
            spotsRlocal(2,:) = spotsRlocal(2,:);
        end
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
    function [newspots,newthresh] = SptFindUserThresh(oldspots,SpotImg,thisn,thisxout,...
            ChName,OrdfiltSize,MinDist,Method,oldthresh)
        newspots = oldspots;
        newthresh = oldthresh;
        happy = 0;
        while happy==0
            answer = input('Press enter if satisfied, anything else+enter if not:','s');
            if ~isempty(answer)
                close all
                figure,bar(thisxout,thisn)
                hold on
                plot([newthresh newthresh],[0 max(thisn)],'--k')
                hold off
                title('Choose a threshold between background and true spots:','Fontsize',12)
                xlabel('<-Background intensities ... Real spot intensities->')
                ylabel('Counts')
                newthresh = input('Enter new threshold to use:');
                % Error handling:
                while isempty(newthresh) ||  newthresh <=0 || newthresh >= 1
                    newthresh = input('Enter new threshold to use:');
                end
                close
                [newspots,thisn,thisxout,newthresh] = FindSpotsV5(SpotImg,'ShowResults',1,...
                    'UserThresh',newthresh,'ImgTitle',ChName,'NeighborhoodSize',OrdfiltSize,...
                    'maxsize',MinDist,'Method',Method);
            else
                happy = 1;
            end
        end
        close
    end

    % Function 3: After calculating transform, eliminate any mispaired or
    % badly fit spots and refit, until the fit residuals are sufficiently
    % low. Note matchedG, matchedR need to be cell arrays of matrices. 
    function [newtform,matchedG,matchedR] = RefineTform(matchedG,matchedR,oldtform,params)
        % Update 4/2014: Any mis-pairings of spots has a big effect on the
        % quality of the fitted transform, even for ~750 spots. With
        % perfect matching and our current alignment, the residuals should
        % be <0.008 per spot (note that the residuals will increase with more spots!):
        PrevResid = oldtform.ResidualsFwd;
        ResidGtoR = oldtform.ResidualsFwd;
        ResidRtoG = oldtform.ResidualsInv;
        steps = 0;
        newtform = oldtform;
        allmatchesG = [];
        allmatchesR = [];
        for yy = 1:length(matchedG)
            allmatchesG = [allmatchesG,matchedG{yy}];
            allmatchesR = [allmatchesR,matchedR{yy}];
        end
        InitBdNum = size(allmatchesG,2);
        while ResidGtoR/size(allmatchesG,2)>=params.ResidTolerance ||...
            ResidRtoG/size(allmatchesG,2)>=params.ResidTolerance || ...
            abs(ResidGtoR-PrevResid)/PrevResid > 0.05 % Stop if residuals stop changing much
            PrevResid = ResidGtoR;
            tempRs = newtform.FRETmapFwd(allmatchesG);
            CurrErrors = sqrt((allmatchesR(1,:)-tempRs(1,:)).^2+(allmatchesR(2,:)-tempRs(2,:)).^2);
            % I'm not sure why the above two lines aren't identical to:
            %CurrDists = FindSpotDists(tformPoly.FRETmapFwd(matchGall),matchRall);
            %CurrErrors = min(CurrDists,[],2);
            % Using a threshold based only on the G to R transformation,
            % because while the residuals for the two directions do usually
            % differ, it's not usaully by much
            MismatchThresh = mean(CurrErrors)+5*std(CurrErrors);
            clear tempRs

            allmatchesG = [];
            allmatchesR = [];
            for p = 1:length(matchedG)
                tempG = matchedG{p};
                tempR = matchedR{p};
                newRs = newtform.FRETmapFwd(tempG);
                tempG = tempG(:,sqrt((tempR(1,:)-newRs(1,:)).^2+(tempR(2,:)-newRs(2,:)).^2)<=MismatchThresh);
                tempR = tempR(:,sqrt((tempR(1,:)-newRs(1,:)).^2+(tempR(2,:)-newRs(2,:)).^2)<=MismatchThresh);
                matchedG{p} = tempG;
                matchedR{p} = tempR;
                allmatchesG = [allmatchesG, matchedG{p}];
                allmatchesR = [allmatchesR, matchedR{p}];
                clear newRs tempR tempG
            end

            newtform = FRETmap(allmatchesG,allmatchesR,'Green',params.TransformToCalc,...
                params.TformMaxDeg,params.TformTotDeg);
            disp(sprintf('Residuals for %d spots:',size(allmatchesG,2)))
            ResidGtoR = newtform.ResidualsFwd
            ResidRtoG = newtform.ResidualsInv
            
            if size(allmatchesG,2)<0.75*InitBdNum || steps>5
                newtform.HistResiduals('fwd');
                disp('Channel mapping: having to exclude lots of beads to get residuals down.')
                disp('If histogram looks ok, you can increase tolerance, for example, enter:')
                disp('params.ResidTolerance = ResidualsGtoR/size(matchGall,2)+0.001; dbcont')
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
    end

    % Load and check parameter settings
    smFRETsetup;
    params = load('AnalysisParameters.mat');

%%%%%%FIRST PART: Channel mapping:
    % Load an old channel mapping, or perform a new one:
    DoMap = input('Press enter to perform channel mapping with beads, L+enter to load an old one, D+enter to use data to make a map:','s');

    if strcmpi(DoMap,'L')
        % Default to most recent map:
        if exist('PathToRecentMap.mat','file')
            prevmapdir = load('PathToRecentMap.mat');
        else
            prevmapdir.MostRecentMapDir = params.defaultdatadir;
        end
        D_Beads = uigetdir(prevmapdir.MostRecentMapDir,'Select directory with old map');
        if exist(fullfile(D_Beads,'ChannelMapping.mat'),'file')
            Map = load(fullfile(D_Beads,'ChannelMapping.mat'));
            tformPoly = Map.tformPoly;
            tformAffine = Map.tformAffine;
            MappingTolerance = Map.MappingTolerance;
            if isfield(Map,'PxlsToExclude')
                params.PxlsToExclude = Map.PxlsToExclude;
            end
            clear Map prevmapdir
        else
            disp(strcat('Map not found in',D_Beads))
            return
        end
    else
        if isempty(DoMap)
            D_Beads = uigetdir(params.defaultdatadir,'Select directory with beads');
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
                    [spotsG{i},n,xout,thresh] = FindSpotsV5(imgGreen,'ShowResults',1,'ImgTitle','Green Channel',...
                        'NeighborhoodSize',params.BeadNeighborhood,'maxsize',params.BeadSize,...
                        'Method','GaussFit');
                    spotsG{i} = SptFindUserThresh(spotsG{i},imgGreen,n,xout,'Green Channel',...
                        params.BeadNeighborhood,params.BeadSize,'GaussFit',thresh);
                else
                    [spotsG{i},n,xout,thresh] = FindSpotsV5(imgGreen,'ShowResults',1,'ImgTitle','Green Channel',...
                        'NeighborhoodSize',params.BeadNeighborhood,'maxsize',params.BeadSize);
                    spotsG{i} = SptFindUserThresh(spotsG{i},imgGreen,n,xout,'Green Channel',...
                        params.BeadNeighborhood,params.BeadSize,'default',thresh);
                end
                clear n xout

                % Find spots in red channel
                if params.Refine_Bd_Cen
                    [spotsR{i},n,xout,thresh] = FindSpotsV5(imgRed,'ShowResults',1,'ImgTitle','Red Channel',...
                        'NeighborhoodSize',params.BeadNeighborhood,'maxsize',params.BeadSize,...
                        'Method','GaussFit');
                    spotsR{i} = SptFindUserThresh(spotsR{i},imgRed,n,xout,'Red Channel',...
                        params.BeadNeighborhood,params.BeadSize,'GaussFit',thresh);
                else
                    [spotsR{i},n,xout,thresh] = FindSpotsV5(imgRed,'ShowResults',1,'ImgTitle','Red Channel',...
                        'NeighborhoodSize',params.BeadNeighborhood,'maxsize',params.BeadSize);
                    spotsR{i} = SptFindUserThresh(spotsR{i},imgRed,n,xout,'Red Channel',...
                        params.BeadNeighborhood,params.BeadSize,'default',thresh);
                end
                clear n xout

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
                if matchG{i} == -1 % Greedy algorithm failed
                    [matchR{i},matchG{i}] = UserPickSptsForAffine(imgRed,...
                        imgGreen,spotsR{i},spotsG{i},params.BeadSize);
                end
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
                        ylim([-size(imgRed(:,:,i),1) 0])
                        xlim([0 size(imgRed(:,:,i),2)/2])
                    else
                        ylim([-size(imgRed(:,:,i),1)/2 0])
                        xlim([0 size(imgRed(:,:,i),2)])
                    end
                    title('Red x is center of point in red channel')
                    pause
                    close all
                    clear matchG_abs matchR_abs
                end

                clear TotImg TotImg2 imgGreen imgRed spotsGabs
            end

        else
            % Load paired DNA spots:
            D_PairedDNAs = uigetdir(params.defaultsavedir,'Select directory with paired DNAs');
            PairedSpots = dir(fullfile(D_PairedDNAs,'Spot*.mat'));
            matchGall = [];
            matchRall = [];
            for d = 1:length(PairedSpots)
                if ~isempty(regexpi(PairedSpots(d).name,'Spot\d+.*mat'))
                    tempspot = load(fullfile(D_PairedDNAs,PairedSpots(d).name));
                    matchGall(:,end+1) = tempspot.Gspot;
                    matchRall(:,end+1) = tempspot.Rspot;
                end
            end
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

        % Refine map by removing outliers (caused by mismatched pairs)
        if isempty(DoMap)
            for gg = 1:num_BeadDir
                matchGforTform{gg} = matchG{gg};
                matchRforTform{gg} = matchR{gg};
            end
            [tformPoly,matchGforTform,matchRforTform] = RefineTform(matchGforTform,matchRforTform,tformPoly,params);
            matchGall = [];
            matchRall = [];
            for gg = 1:num_BeadDir
                matchG{gg} = matchGforTform{gg};
                matchR{gg} = matchRforTform{gg};
                matchGall = [matchGall matchGforTform{gg}];
                matchRall = [matchRall matchRforTform{gg}];
            end
        else
            matchG{1} = matchGall;
            matchR{1} = matchRall;
            [tformPoly,matchG,matchR] = RefineTform(matchG,matchR,tformPoly,params);
            matchGall = matchG{1};
            matchRall = matchR{1};
        end
        tformAffine = FRETmap(matchGall,matchRall,'Green','Affine');

        % Finally, let the user check the resulting transformation:
        if isempty(DoMap)
            % Plot the results for each movie:
            for i = 1:num_BeadDir
                disp(strcat('Iterating through bead images for user to check quality (',...
                    int2str(i),' of ',int2str(num_BeadDir),')'))
                % Use tform to map to where they should be in the red channel:
                newR = tformPoly.FRETmapFwd(matchG{i});
                % Get points in absolute coordinates for overlaying on the
                % whole image
                [matchG_abs,newR_abs] = SpotsIntoAbsCoords(matchG{i},...
                    newR,params,size(allBdImgs(:,:,i),2)/2);
                PutBoxesOnImageV4(allBdImgs(:,:,i),[newR_abs';matchG_abs'],params.BeadSize);
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
                % Get points in absolute coordinates
                [newG_abs,matchR_abs] = SpotsIntoAbsCoords(newG,...
                    matchR{i},params,size(allBdImgs(:,:,i),2)/2);
                PutBoxesOnImageV4(allBdImgs(:,:,i),[matchR_abs';newG_abs'],params.BeadSize);
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
        else
            tformPoly.PlotTform(tformPoly.FRETmapFwd(tformPoly.StartData),tformPoly.EndData)
            legend('Green spots mapped to red','Red spots')
            ylim([-512 0])
            xlim([0 256])
            tformPoly.HistResiduals('fwd')
            xlabel('Distance between green spots mapped to red channel, and real red spots','Fontsize',12)
            pause
            close all
        end

        % The mapping tolerance is the maximal distance away a spot center
        % could be in the other channel, and still possibly be the same as
        % the spot you're looking at:
        tempGinRs = tformPoly.FRETmapFwd(matchGall);
        tempRinGs = tformPoly.FRETmapInv(matchRall);
        DistsG = sqrt((matchRall(1,:)-tempGinRs(1,:)).^2+(matchRall(2,:)-tempGinRs(2,:)).^2);
        DistsR = sqrt((matchGall(1,:)-tempRinGs(1,:)).^2+(matchGall(2,:)-tempRinGs(2,:)).^2);
        MappingTolerance = ceil(max(mean(DistsG)+5*std(DistsG),mean(DistsR)+5*std(DistsR)))

        % Lastly, check the transformation with held-out data, if we have some:
        if isempty(DoMap) && strcmpi(checkTform,'y')
            disp('Residuals with held-out data:')
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

        clear allBdImgs matchG matchR matchGall matchRall DistsG DistsR tempRinGs tempGinRs
        clear num_BeadDir BdDir checkTform

        PxlsToExclude = params.PxlsToExclude;
        if isempty(DoMap)
            save(fullfile(D_Beads,'ChannelMapping.mat'),'tformPoly','tformAffine',...
                'BeadFilesInMap','MappingTolerance','PxlsToExclude');
        else
            D_Beads = uigetdir(params.defaultdatadir,...
                'Save this map with the data used to generate it');
            save(fullfile(D_Beads,'ChannelMapping.mat'),'tformPoly','tformAffine',...
                'MappingTolerance','PxlsToExclude');
        end
        % To easily load the most recent bead map:
        MostRecentMapDir = D_Beads;
        save('PathToRecentMap','MostRecentMapDir');
        clear MostRecentMapDir
    end
clear DoMap

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

        % Update 12/2013: If this movie has already been analyzed, provide
        % the option to use the previously found spots, instead of
        % re-finding them
        
        useoldI = 'n';
        
        if exist(fullfile(savedir,strcat('SpotsAndIntensities',int2str(i),'.mat')),'file')
            useoldI = input('Load previously found spots and their intensities? (y/n)','s');
        end
            
        if strcmpi(useoldI,'n')
            
            useoldspots = 'n';
            if exist(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),'file')
                useoldspots = input('Load previously found spots? (y/n)','s');
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

               % Update 5/2014: Added the option to find spots throughout the
               % movie, not just from EndInjectFrame:EndInjectFrame+FramesToAvg
               totframes = params.FrameLoadMax*length(dir(fullfile(D_Data,ToAnalyze(i).name,'ScaledMovie*.mat')));
               if params.FindSpotsEveryXFrames==0
                   SptFindIncrement = totframes;
               else
                   SptFindIncrement = params.FindSpotsEveryXFrames;
               end
%                if params.EndInjectFrame == 1
%                    EndInjectFrame = SptFindIncrement+1;
%                else
                   EndInjectFrame = params.EndInjectFrame;
               %end
               RefinedCentersR = [];
               RefinedCentersG = [];
               VarsR = [];
               VarsG = [];
               newspotsR = [];
               newspotsG = [];
               newVarsR = [];
               newVarsG = [];
               %for ff = [1,EndInjectFrame:SptFindIncrement:totframes]
               for ff = EndInjectFrame:SptFindIncrement:totframes
                   [imgRed,imgGreen] = LoadScaledMovie(fullfile(D_Data,ToAnalyze(i).name),...
                       [ff ff+params.FramesToAvg]);
                   imgRedavg = mat2gray(mean(imgRed,3)); %Do I want to do mat2gray here? Update 4/2014:
                        % since I'm going to treat spotfinding as totally separate
                        % from Gauss fitting for intensity smoothing, it is best
                        % that I do mat2gray here
                   imgGreenavg = mat2gray(mean(imgGreen,3));

                   % Step 0: subtract background:
                   %[imgRbkgnd,imgGbkgnd,imgRMinusBkgnd,imgGMinusBkgnd] = SubBkgnd(imgRedavg,imgGreenavg,params);
                   % If you don't want to subtract background, uncomment these
                    % lines:
                   imgRMinusBkgnd = imgRedavg;
                   imgGMinusBkgnd = imgGreenavg;

                   % Step 1: find spots
                   % Find spots in both channels, but don't double-count. Allow
                   % user to decide whether to find spots separately in each
                   % channel, or using a combined image (see notes on the
                   % UseCombinedImage parameter in smFRETsetup.m).
                   if params.UseCombinedImage
                       disp('Reminder: Using a combined image to find spots does not currently work well.')
                       % Finding spots in a composite image, so that
                       % mid-FRET spots don't get lost.  NOTE that the combined image
                       % will have a frame of reference of the acceptor image, which
                       % is fine because that's what I pass into UserSpotSelectionV4.

                       % Step 1.1 Create a combined image
                       imgRMinusBkgnd = CalcCombinedImage(tformAffine,imgGMinusBkgnd,imgRMinusBkgnd);
                       % The built-in Matlab function imfuse used to create the output
                       % of CalcCombinedImage only returns unit8 images, but fminsearch
                       % (called in Fit2DGaussToSpot in GetGaussParams below) needs a 
                       % double, so convert back to doubles:
                       if imgRMinusBkgnd~=-1
                           imgRMinusBkgnd = mat2gray(imgRMinusBkgnd);
                       else
                           disp('smFRET: Calculation of combined image failed.')
                           return
                       end
                   end

               % Step 1.1 Identify spots in acceptor channel, and refine centers by
                   % fitting to a Gaussian, regardless of whether or not user
                   % wants to weight intensities by a Gaussian:
                   %if ff == 1
                   if ff == EndInjectFrame
                       [newspotsR,nR,xoutR,thresh] = FindSpotsV5(imgRMinusBkgnd,'ShowResults',1,'ImgTitle','Red Channel',...
                             'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize,...
                             'Method','GaussFit');
                       [newspotsR,threshholdR] = SptFindUserThresh(newspotsR,imgRMinusBkgnd,nR,xoutR,'Red Channel',...
                           params.DNANeighborhood,params.DNASize,'GaussFit',thresh);
                       clear thresh
                       close
                       % And in donor channel

                       if ~params.UseCombinedImage
                           [newspotsG,nG,xoutG,thresh] = FindSpotsV5(imgGMinusBkgnd,'ShowResults',1,'ImgTitle','Green Channel',...
                                 'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize,...
                                 'Method','GaussFit');
                           [newspotsG,threshholdG] = SptFindUserThresh(newspotsG,imgGMinusBkgnd,nG,xoutG,'Green Channel',...
                               params.DNANeighborhood,params.DNASize,'GaussFit',thresh);
                           clear thresh
                           close
                       end

                       if debug
                           PutBoxesOnImageV4(mat2gray(imgRMinusBkgnd),newspotsR,params.DNASize);
                           title('Red Channel, start')
                           figure,bar(xoutR,nR)
                           hold on
                           plot([threshholdR threshholdR],[0 max(nR)],'--k')
                           hold off
                           title('Red Channel, start')
                           if ~params.UseCombinedImage && ~isempty(newspotsG)
                               PutBoxesOnImageV4(mat2gray(imgGMinusBkgnd),newspotsG,params.DNASize);
                               title('Green Channel, start')
                               figure,bar(xoutG,nG)
                               hold on
                               plot([threshholdG threshholdG],[0 max(nG)],'--k')
                               hold off
                               title('Green Channel, start')
                           end

                       end
                       clear nR xoutR nG xoutG
                   elseif mod(ff-params.EndInjectFrame,params.CheckSpotFindingEveryXFrames)==0
                       % For each subsequent set of frames, don't fit a
                       % Gaussian to find the center during spotfinding,
                       % since it just makes it take longer. Will do that
                       % below for any new spots:
                        [tempnewspotsR,nR,xoutR,threshholdR] = FindSpotsV5(imgRMinusBkgnd,'ShowResults',1,'ImgTitle','Red Channel',...
                            'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize,...
                            'UserThresh',threshholdR,'Method','GaussFit');
                       [tempnewspotsR,threshholdR] = SptFindUserThresh(tempnewspotsR,imgRMinusBkgnd,nR,xoutR,'Red Channel',...
                           params.DNANeighborhood,params.DNASize,'GaussFit',threshholdR);
                       close
                       % And in donor channel
                       if ~params.UseCombinedImage
                           [tempnewspotsG,nG,xoutG,threshholdG] = FindSpotsV5(imgGMinusBkgnd,'ShowResults',1,'ImgTitle','Green Channel',...
                                 'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize,...
                                 'UserThresh',threshholdG,'Method','GaussFit');
                           [tempnewspotsG,threshholdG] = SptFindUserThresh(tempnewspotsG,imgGMinusBkgnd,nG,xoutG,'Green Channel',...
                               params.DNANeighborhood,params.DNASize,'GaussFit',threshholdG);
                           clear thresh
                           close
                       end
                   else
                       % Don't ask the user to check the spotfinding for
                       % this round
                        [tempnewspotsR,nR,xoutR,~] = FindSpotsV5(imgRMinusBkgnd,...
                            'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize,...
                            'UserThresh',threshholdR,'Method','GaussFit');
                        if ~params.UseCombinedImage
                            [tempnewspotsG,nG,xoutG,~] = FindSpotsV5(imgGMinusBkgnd,...
                                'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize,...
                                'UserThresh',threshholdG,'Method','GaussFit');
                        end

                   end
                   %if ff > 1
                   if ff > EndInjectFrame
                        % Are there any new spots that we didn't find last
                        % time?
                        if ~isempty(tempnewspotsR) && ~isempty(RefinedCentersR)
                            Dists = FindSpotDists(RefinedCentersR,tempnewspotsR);
                            spotnottooclose = Dists>params.DNASize;
                            newspotsR = tempnewspotsR(:,sum(spotnottooclose,1)==size(spotnottooclose,1));
                            clear tempnewspotsR
                            clear Dists spotnottooclose
                        end
                        if ~params.UseCombinedImage && ...
                                ~isempty(RefinedCentersG) && ~isempty(tempnewspotsG)
                            Dists = FindSpotDists(RefinedCentersG,tempnewspotsG);
                            spotnottooclose = Dists>params.DNASize;
                            newspotsG = tempnewspotsG(:,sum(spotnottooclose,1)==size(spotnottooclose,1));
                            clear tempnewspotsG
                            clear Dists spotnottooclose
                        end
                   end

                        if debug && ff>1 && ...
                                mod(ff-params.EndInjectFrame,params.CheckSpotFindingEveryXFrames)==0
                            if ~isempty(newspotsR)
                               PutBoxesOnImageV4(mat2gray(imgRMinusBkgnd),newspotsR,params.DNASize);
                            end
                           title('Red Channel, current')
                           figure('Position',[200 200 325 625])
                           hold on
                           if ~isempty(RefinedCentersR)
                            plot(RefinedCentersR(2,:),0-RefinedCentersR(1,:),'xg')
                           end
                           if ~isempty(newspotsR)
                            plot(newspotsR(2,:),0-newspotsR(1,:),'xr')
                           end
                           hold off
                            ylim([-512 0])
                            xlim([0 256])
                            title('Red Channel,current')
                            if ~params.UseCombinedImage 
                                if ~isempty(newspotsG)
                                   PutBoxesOnImageV4(mat2gray(imgGMinusBkgnd),newspotsG,params.DNASize);
                                end
                               title('Green Channel, current')
                               figure('Position',[200 200 325 625])
                               hold on
                               if ~isempty(RefinedCentersG)
                                    plot(RefinedCentersG(2,:),0-RefinedCentersG(1,:),'xg')
                               end
                               if ~isempty(newspotsG)
                                    plot(newspotsG(2,:),0-newspotsG(1,:),'xr')
                               end
                                hold off
                                ylim([-512 0])
                                xlim([0 256])
                                title('Green Channel,current')

                                figure,bar(xoutR,nR)
                                hold on
                                plot([threshholdR threshholdR],[0 max(nR)],'--k')
                                hold off
                                title('Red Channel,current')
                                figure,bar(xoutG,nG)
                                hold on
                                plot([threshholdG threshholdG],[0 max(nG)],'--k')
                                hold off
                                title('Green Channel,current')
                            end
                            pause
                            close
                            close
                            close
                            if ~params.UseCombinedImage
                                close
                                close
                                close
                            end
                        end
                    clear nR xoutR nG xoutG

                   % Update 5/2014: As noted above, stopped fitting a gaussian
                   % to subsequent sets of frames, since you'll be fitting
                   % largely all the spots you found last round.  Doing that
                   % here instead.  
                   % Some notes on variances: It looks like the fits are much more
                       % likely to be poor if I use single frames, rather than
                       % an average of ~10 frames, to find the variances. It's
                       % also a lot slower to do it frame by frame (probably
                       % because of all the poor fits!).
                       % Doesn't seem to matter a ton to use the background
                       % subtracted image or not.
                       % Lastly, red and green channels have on average roughly
                       % the same variances, so I think it's probably fine to
                       % assume the same variance for a spot in the red channel
                       % as in the green channel.  But I should check by
                       % finding all spots throughout the movie, pairing, and
                       % comparing the variances of true pairs ...
                   if ~isempty(newspotsR)
                       [newspotsRref, newVarsR] = FindRefinedSpotCenters(imgRedavg,newspotsR,0.2,params);
                       clear newspotsR
                       newspotsR = newspotsRref;
                       clear newspotsRref
                   end
                   if ~params.UseCombinedImage && ~isempty(newspotsG)
                       [newspotsGref, newVarsG] = FindRefinedSpotCenters(imgGreenavg,newspotsG,0.2,params);
                       clear newspotsG
                       newspotsG = newspotsGref;
                       clear newspotsGref
                   end

                   RefinedCentersR = [RefinedCentersR,newspotsR];
                   RefinedCentersG = [RefinedCentersG,newspotsG];
                   VarsR = [VarsR,newVarsR];
                   VarsG = [VarsG,newVarsG];
                   newspotsR = []; 
                   newspotsG = []; 
                   newVarsR = []; 
                   newVarsG = [];
                   clear imgRbkgnd imgGbkgnd imgRMinusBkgnd imgGMinusBkgnd imgGreen imgRedavg imgGreenavg

                   if mod((ff-params.EndInjectFrame),100)==0
                       if ~params.UseCombinedImage
                           disp(sprintf('%d red spots and %d green spots found after %d frames',...
                               size(RefinedCentersR,2),size(RefinedCentersG,2),ff+params.FramesToAvg))
                       else
                           disp(sprintf('%d spots found after %d frames',...
                               size(RefinedCentersR,2),ff+params.FramesToAvg))
                       end
                   end

               end
               close all
               clear totframes SptFindIncrement threshhold threshholdR threshholdG
               clear newspotsR newspotsG newVarsR newVarsG
               
               save(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),...
                   'RefinedCentersR','RefinedCentersG','VarsR','VarsG');

               if ~params.UseCombinedImage
                   % Step 1.3: Add any spots in green channel that weren't
                   % paired to a red channel spot to our list of spots. Or, put
                   % equivalently, add everything in RefinedCentersG that's not
                   % in matchG.
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
                   if ~isempty(RefinedCentersG)
                       spotsGinR = tformPoly.FRETmapFwd(RefinedCentersG);
                       Dists = FindSpotDists(RefinedCentersR,spotsGinR);
                       spotnottooclose = Dists>MappingTolerance*2;
                       newspots = spotsGinR(:,sum(spotnottooclose,1)==size(spotnottooclose,1));
                       clear spotsGinR
                       spotsGinR = newspots;
                       Vars = VarsG(:,sum(spotnottooclose,1)==size(spotnottooclose,1));
                   else
                       spotsGinR = [];
                   end

                   % Check that the transformed G spots are reasonable
                   % edges from the red channel boundaries: note that the same
                   % is done for red spots in CalcIntensitiesV2
                   if ~isempty(spotsGinR)
                       [spotsGinR,~,Vars,~] = CheckSpotBoundaries(spotsGinR,...
                            spotsGinR,Vars,Vars,params,fullfile(D_Data,ToAnalyze(i).name));
                   else
                       Vars = [];
                   end
                   Vars(:,end+1:end+size(VarsR,2)) = VarsR;
                   spots = [spotsGinR, RefinedCentersR];
                   clear spotsGinR
                   
               else
                  % From now on I have to use the affine transformation to
                   % map from one channel to the other:
                   tformPoly = tformAffine;
                   Vars = VarsR;
                   spots = RefinedCentersR;
               end

               disp(sprintf('Found %d total spots',size(spots,2)))
               save(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),...
                   'spots','Vars','-append');
            else
               prevspots = load(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')));
               spots = prevspots.spots;
               Vars = prevspots.Vars;
               clear prevspots
            end
           
           % Step 2: Load the whole movie in increments and calculate the
           % intensity of each spot in each frame.
           
           disp('Calculating frame-by-frame intensities')
           
           [RedI, GrI] = CalcIntensitiesV3(fullfile(D_Data,ToAnalyze(i).name),...
               spots, Vars, tformPoly,params);
           
           % Save spot positions, intensities and associated GaussFit
           % parameters in case the user wants to re-analyze.
           % Note because spots are never displayed on the full image,
           % only on the image split into two channels, I don't need to
           % add params.PxlsToExclude to get spots into the right
           % coordinates (unlike with the beads)
           
           SpotsInR = spots;
           SpotVars = Vars;
           save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(i),'.mat')),'SpotsInR',...
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
           if ~params.IntensityGaussWeight 
               % UserSpotSelectionV4 plots 5*the std dev of the weighting Gaussian as a circle
               % around the spot center; if the user didn't weight
               % intensities by a Gaussian, show instead a circle of radius
               % 5 over which the intensity was summed:
               % (It should be noted that what I call "Vars", "spotVars",
               % etc is actually 1/(2*the real variance))
               Vars = ones(size(spots));
           elseif ~isempty(params.FixSpotVar)
               Vars = repmat(params.FixSpotVar,1,size(spots,2));
           end
               UserSpotSelectionV4(RedI,GrI,spots,Vars,...
                   fullfile(D_Data,ToAnalyze(i).name),params,tformPoly,savedir,fps,i);
        else %If the user wants to instead use previously saved data
           oldspots = load(fullfile(savedir,strcat('SpotsAndIntensities',int2str(i),'.mat')));
           useoldparams = input('Use old parameters? (y/n)','s');
           if strcmpi(useoldparams,'y')
                params = load(fullfile(savedir,strcat('AnalysisParameters.mat')));
                params = params.params;
           else
               params.fps = fps;
               save(fullfile(savedir,strcat('AnalysisParameters.mat')),'params');
           end
           disp(strcat('Movie ',int2str(i),'of',int2str(length(ToAnalyze))))
           if params.IntensityGaussWeight
               if isempty(params.FixSpotVar)
                   Vars = oldspots.SpotVars;
               else 
                   Vars = repmat(params.FixSpotVar,1,size(oldspots.SpotsInR,2));
               end
           else
               Vars = ones(size(spots));
           end
           UserSpotSelectionV4(oldspots.RedI,oldspots.GrI,oldspots.SpotsInR,...
               Vars,fullfile(D_Data,ToAnalyze(i).name),params,tformPoly,savedir,fps,i);
        end
        clear TotImg spots imgRed imgGreen spotsG spotsR spotsG_abs spotsRguess spotstemp
    end
end
