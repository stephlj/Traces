% function smFRET(rootname,debug)
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
% Rootname: root filename (the movies it analyzes will be rootname_1, rootname_2 etc)
% Optionally: if you want it to run all the various debugging things, pass
% "1" as the second input.
%
% Copyright (C) 2014 Stephanie Johnson, University of California, San Francisco
% Contact: Stephanie.Johnson@ucsf.edu
% 
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
            answer = input('Press enter if satisfied, anything else+enter if not: ','s');
            if ~isempty(answer)
                close all
                figure,bar(thisxout,thisn)
                hold on
                plot([newthresh newthresh],[0 max(thisn)],'--k')
                hold off
                title('Choose a threshold between background and true spots:','Fontsize',12)
                xlabel('<-Background intensities ... Real spot intensities->')
                ylabel('Counts')
                % Update 10/2014: at Kai-Chun's suggestion, printing out
                % the old threshold
                disp(sprintf('Current threshold: %d',newthresh))
                newthresh = input('Enter new threshold to use: ');
                % Error handling:
                while isempty(newthresh) ||  newthresh <=0 || newthresh >= 1
                    newthresh = input('Enter new threshold to use: ');
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
        if oldtform.A == -1
            newtform = -1;
            matchedG = -1;
            matchedR = -1;
            return
        end
        StdDevMultiplier = 5;
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
            MismatchThresh = mean(CurrErrors)+StdDevMultiplier*std(CurrErrors);
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
            ResidGtoRperSpot = newtform.ResidualsFwd/size(allmatchesG,2)
            
            if size(allmatchesG,2)<0.75*InitBdNum || steps>5
                newtform.HistResiduals('fwd');
                disp('Channel mapping: having to exclude lots of beads to get residuals down.')
                disp('If histogram looks ok, you can increase tolerance, for example, enter:')
                disp('params.ResidTolerance = ResidGtoR/size(allmatchesG,2)+0.001; dbcont')
                disp('If you want it to remove more spots with higher residuals, reduce StdDevMultiplier')
                disp('(e.g., enter StdDevMultiplier = 3; dbcont into the command line)')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% First PART: Channel Mapping: %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            disp('Map not found.')
            return
        end
    else
        if isempty(DoMap)
            D_Beads = uigetdir(params.defaultdatadir,'Select directory with beads');
            % Figure out how many bead files to analyze there are:
            AllBeads = dir(fullfile(D_Beads,'Bead*'));
            num_BeadDir = input(sprintf('How many bead files to use for transformation? Max: %d (Enter to use max) ',...
                length(AllBeads)));
            if isempty(num_BeadDir) || num_BeadDir >= length(AllBeads)
                BdDir = length(AllBeads);
                num_BeadDir = length(AllBeads);
                checkTform = 'n';
            elseif num_BeadDir < length(AllBeads)
                checkTform = input(sprintf('Use remaining %d movies to check transform? (y/n) ',...
                    length(AllBeads)-num_BeadDir),'s');
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

                TotImg = LoadRawImgs(fullfile(D_Beads,AllBeads(i).name),'FramesToLoad',[1 params.FramesToAvg]);

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
                    [spotsG_abs,spotsR_abs] = SpotsIntoAbsCoords(spotsG{i},...
                            spotsR{i},params,size(allBdImgs(:,:,i),2)/2);
                    PlotDebugFigures(1,allBdImgs(:,:,i),spotsR_abs,spotsG_abs,params)
                    clear spotsG_abs spotsR_abs
                end

                % Step 2: figure out which spots in one channel go with the spots in
                % the other channel.  For our beads, which are brighter in the donor
                % than acceptor channel, the matching works best if you match spots
                % in donor to spots in acceptor:
                [matchG{i},matchR{i}] = FindSpotMatches(spotsG{i},spotsR{i});
                if matchG{i} == -1 % Greedy algorithm failed
                    if i~=1 % If we've already gotten through one data set, so we have a temporary map:
                        disp('Using map generated from previous data instead:')
                        [newspotsGtemp,newspotsR] = FindSpotMatches(tformtemp.FRETmapFwd(spotsG{i}),spotsR{i});
                        if size(newspotsGtemp,1) == 1 && newspotsGtemp == -1
                            [matchR{i},matchG{i}] = UserPickSptsForAffine(imgRed,...
                                imgGreen,spotsR{i},spotsG{i},params.BeadSize);
                        else
                            newspotsG = tformtemp.FRETmapInv(newspotsGtemp);
                            matchR{i} = newspotsR;
                            matchG{i} = newspotsG;
                            clear newspotsR newspotsGtemp
                        end
                    else
                        [matchR{i},matchG{i}] = UserPickSptsForAffine(imgRed,...
                            imgGreen,spotsR{i},spotsG{i},params.BeadSize);
                    end
                end
                if i<=num_BeadDir
                    matchGall = [matchGall, matchG{i}];
                    matchRall = [matchRall, matchR{i}];
                end

                if debug % Box in green the ones that were matched:
                    [matchG_abs,matchR_abs] = SpotsIntoAbsCoords(matchG{i},...
                            matchR{i},params,size(allBdImgs(:,:,i),2)/2);
                    PlotDebugFigures(2,allBdImgs(:,:,i),matchG_abs,matchR_abs,params,...
                        matchR{i},matchG{i},size(imgRed))
                    clear matchG_abs matchR_abs
                end
                
                % Update 11/5/2014: creating a temporary map that can be
                % used if the greedy algorithm fails on the next data
                % set(s)
                tformtemp = FRETmap(matchGall,matchRall,'Green',params.TransformToCalc,...
                    params.TformMaxDeg,params.TformTotDeg);

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
        clear tformtemp
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
            if tformPoly == -1
                return
            end
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
                legend('Red spots mapped to green','Green spots')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% SECOND PART: Analyze data: %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
    D_Data = uigetdir(D_Beads,'Select directory with data to analyze');
    % Figure out how many DNA files to analyze there are:
    ToAnalyze = dir(fullfile(D_Data,strcat(rootname,'*')));
    % Get framerate for plotting:
    if isempty(ToAnalyze)
        disp('Did not find data to analyze.') %error handling
        return
    end
    fps = GetInfoFromMetaData(fullfile(D_Data,ToAnalyze(1).name),'fps');
    fps = 1/fps; % This is actually frames per ms
    
    % Make sure not saving over old data:
    if ~exist(fullfile(params.defaultsavedir,rootname),'dir')
        if exist(params.defaultsavedir,'dir')
            savedir = fullfile(params.defaultsavedir,rootname);
        else
            savedirroot = uigetdir('','Save data where?:');
            savedir = fullfile(savedirroot,rootname);
        end
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
                useoldTform = 'n';
                if strcmpi(useoldspots,'y')
                    useoldTform = input('Use same channel mapping as last time but recompute: scaling, background, intensities? (y/n)','s');
                end
            end
            
            % Option 1: keep all old spots (so that data will have all the
            % same labels as last time), just re-scale and recompute
            % background and intensities
            if strcmpi(useoldspots,'y') && strcmpi(useoldTform,'y')
               prevRspots = load(fullfile(savedir,strcat('SpotsAndIntensities',int2str(i),'.mat')));
               spots = prevRspots.SpotsInR;
               % If the user re-located any green spots, want to use those
               % locations:
               if isfield(prevRspots,'SpotsInG')
                    SpotsInG = prevRspots.SpotsInG;
               end
               Vars = prevRspots.SpotVars;
               clear prevRspots
               disp('Scaling movie ...')
               ScaleMovieV2(fullfile(D_Data,ToAnalyze(i).name),params);
               UseScaledMov = 'n'; % This is so the background will get re-computed, below
            
            % Option 2: either re-map, or find all new spots and then re-map
            else
               % Find all new spots
               if strcmpi(useoldspots,'n')
                   % Finding spots and local background values

                   % Update 4/2014: Since the movie scaling now takes a while, provide
                   % the option to load the previously scaled movie, but a separate
                   % option to re-find spots

                   UseScaledMov = 'n';

                   if exist(fullfile(D_Data,ToAnalyze(i).name,strcat('ScalingInfo.mat')),'file')
                       UseScaledMov = input('Load scaled movie? (y/n)','s');
                   end

                   if strcmpi(UseScaledMov,'n')
                       disp('Scaling movie ...')
                       ScaleMovieV2(fullfile(D_Data,ToAnalyze(i).name),params);
                   end

                   % Update 5/2014: Added the option to find spots throughout the movie
                   [~,totframes] = LoadRawImgs(fullfile(D_Data,ToAnalyze(i).name),...
                       'FramesToLoad',[1 1]);
                   if params.FindSpotsEveryXFrames==0
                       SptFindIncrement = totframes;
                   else
                       SptFindIncrement = params.FindSpotsEveryXFrames;
                   end
                   % This option was to always find spots in the first
                   % FramesToAvg frames, even if EndInjectFrame is bigger than
                   % 1. Decided not to allow this.
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
                   % This goes with the option of insisting spotfinding happen
                   % in the first FramesToAvg frames, regardless of the value
                   % of EndInjectFrame.
                   %for ff = [1,EndInjectFrame:SptFindIncrement:totframes]
                   for ff = EndInjectFrame:SptFindIncrement:totframes
                       [imgRed,imgGreen] = LoadScaledMovie(fullfile(D_Data,ToAnalyze(i).name),...
                           [ff ff+params.FramesToAvg],params);
                       imgRedavg = mat2gray(mean(imgRed,3)); 
                       imgGreenavg = mat2gray(mean(imgGreen,3));

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
                           % of CalcCombinedImage only returns uint8 images, but fminsearch
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

                           clear nR xoutR nG xoutG

                       elseif mod(ff-params.EndInjectFrame,params.CheckSpotFindingEveryXFrames)==0
                           % This is one of the frames at which the user wants
                           % to check the threshold value
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
                               PlotDebugFigures(3,RefinedCentersR,newspotsR)
                               PlotDebugFigures(3,RefinedCentersG,newspotsG)
                            end
                        clear nR xoutR nG xoutG

                       % For any new spot found, calculate its variance: 
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

                   % Save all spots found in both channels:
                   save(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),...
                       'RefinedCentersR','RefinedCentersG','VarsR','VarsG');
               else
                   prevspots = load(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')));
                   %spots = prevspots.spots;
                   %Vars = prevspots.Vars;
                   RefinedCentersR = prevspots.RefinedCentersR;
                   RefinedCentersG = prevspots.RefinedCentersG;
                   VarsR = prevspots.VarsR;
                   VarsG = prevspots.VarsG;
                   clear prevspots
               end
             % Now that all DNA spots are found, make a single list of them (in the
             % coordinates of the acceptor channel) to pass into the function that
             % calculates intensities:
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
                   % edges from the red channel boundaries: 
                   if ~isempty(spotsGinR)
                       [spotsGinR,~,Vars,~] = CheckSpotBoundaries(spotsGinR,...
                            spotsGinR,Vars,Vars,params,fullfile(D_Data,ToAnalyze(i).name));
                   else
                       Vars = [];
                   end
                   % Do the same check for transformed R spots and the
                   % green channel boundaries:
                   if ~isempty(RefinedCentersR)
                       spotsRinG = tformPoly.FRETmapInv(RefinedCentersR);
                       [spotsRinG,~,VarsR,~] = CheckSpotBoundaries(spotsRinG,...
                            spotsRinG,VarsR,VarsR,params,fullfile(D_Data,ToAnalyze(i).name));
                       clear RefinedCentersR
                       RefinedCentersR = tformPoly.FRETmapFwd(spotsRinG);
                   end
                   Vars(:,end+1:end+size(VarsR,2)) = VarsR;
                   spots = [spotsGinR, RefinedCentersR];
                   clear spotsGinR spotsRinG
                   
               else
                  % From now on I have to use the affine transformation to
                   % map from one channel to the other:
                   tformPoly = tformAffine;
                   Vars = VarsR;
                   spots = RefinedCentersR;
               end

               disp(sprintf('Keeping %d total spots',size(spots,2)))
               save(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),...
                   'spots','Vars','-append');
               
            end
           
           % Step 2: Load the whole movie in increments and calculate the
           % intensity of each spot in each frame.
           
           % Update 8/2014: Calculating background here--but only if you
           % have to (i.e. if these files don't exist or you re-scaled the
           % movie; whether you re-mapped spots or whatever doesn't matter
           % for the background)
           tempfiles = dir(fullfile(D_Data,ToAnalyze(i).name,'BackgroundImgs*.mat'));
           tempinfo = load(fullfile(D_Data,ToAnalyze(i).name,strcat('ScalingInfo.mat')));
           if isempty(tempfiles) || (exist('UseScaledMov','var') && strcmpi(UseScaledMov,'n'))
               ComputeBackgroundImgs(fullfile(D_Data,ToAnalyze(i).name),params)
           elseif params.ScaleChannelsSeparately ~= tempinfo.ScaleChannelsSeparately || ...
                   params.NormImage ~= tempinfo.NormImage
               disp('NormImage and/or ScaleChannelsSeparately do not match what was used to calculate background.')
               redobkgnd = input('Recalculate background with new parameters? (y/n): ','s');
               if strcmpi(redobkgnd,'n')
                   params.NormImage = tempinfo.NormImage;
                   params.ScaleChannelsSeparately = tempinfo.ScaleChannelsSeparately;
               else
                   ComputeBackgroundImgs(fullfile(D_Data,ToAnalyze(i).name),params)
               end
               clear redobkgnd
           end
           clear tempinfo tempfiles
           
           disp('Calculating frame-by-frame intensities ... ')
           
           if ~exist('SpotsInG','var') 
               [RedI, GrI] = CalcIntensitiesV3(fullfile(D_Data,ToAnalyze(i).name),...
                   spots, Vars, tformPoly,params);
           else % Means the user wanted just to re-scale and recalculate background,
               % etc--so keep any re-located green spots!
               disp('Calculating donor and acceptor channel intensities separately,') 
               disp('to use re-located green spot locations from previous analysis')
               [RedI, ~] = CalcIntensitiesV3(fullfile(D_Data,ToAnalyze(i).name),...
                   spots, Vars, -1,params);
               [~, GrI] = CalcIntensitiesV3(fullfile(D_Data,ToAnalyze(i).name),...
                 SpotsInG, Vars, 1,params);
           end
           
           % Save spot positions, intensities and associated GaussFit
           % parameters in case the user wants to re-analyze.
           % Note because spots are never displayed on the full image,
           % only on the image split into two channels, I don't need to
           % add params.PxlsToExclude to get spots into the right
           % coordinates (unlike with the beads)
           
           SpotsInR = spots;
           SpotVars = Vars;
           % Update 1/2015: Since I added a feature in UserSpotSelection
           % that saves all the changes the user made to a trace (like
           % xlims, etc), I don't want to completely save over the
           % SpotsAndIntensities file if the user didn't re-find spots or
           % re-match them with a new transformation. I want to keep the
           % same spot indices in that case:
           if strcmpi(useoldspots,'y') && strcmpi(useoldTform,'y')
               save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(i),'.mat')),...
                   'RedI','GrI','-append') % Append overwrites the indicated variables, if they already exist
           else
               save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(i),'.mat')),'SpotsInR',...
                   'SpotVars','RedI','GrI');
           end
           clear SpotsInR SpotVars SpotsInG
           
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
               Vars = repmat(params.FixSpotVar',1,size(spots,2));
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
                   Vars = repmat(params.FixSpotVar',1,size(oldspots.SpotsInR,2));
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
