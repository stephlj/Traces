%function smFRET
%
%Wrapper function for analyzing smFRET data--loads data, calls the
%functions that do the analysis, etc.
%
%This will analyze all the movies in a set at once, where a set is defined
%as having the same root filename with _1, _2 at the end (the way
%MicroManager saves multiple movies of the same root filename).  They all
%have to be in the same folder.
%
%Inputs are:
%Root filename (the movies it analyzes will be rootname_1, rootname_2 etc)
%Optionally: if you want it to run all the various debugging things, pass
%"1" as the second input.
%
%Steph 9/2013
%Copyright 2013 Stephanie Johnson, University of California, San Francisco

function smFRET(rootname,debug)

%%%%%%Preliminaries:
    %Since debug is an optional input, set a default
    if ~exist('debug','var') 
        debug=0; 
    elseif debug == 1
        disp('Right now debug assume channels are split L-R and acceptor is on the L!')
    end

    %Nested functions for use later:
    %Function 1: Getting red and green channels out of a combined image:
    function [imgR,imgG] = SplitImg(img,params_struct)
        if params_struct.splitx
                if ~params_struct.Acceptor
                    imgG = img(:,(size(img,2)/2)+1:end,:); 
                    imgR = img(:,1:(size(img,2)/2),:);
                else
                    imgR = img(:,(size(img,2)/2)+1:end,:); 
                    imgG = img(:,1:(size(img,2)/2),:);
                end
            else
                if ~params_struct.Acceptor
                    imgG = img((size(img,2)/2)+1:end,:,:); 
                    imgR = img(1:(size(img,2)/2),:,:);
                else
                    imgR = img((size(img,2)/2)+1:end,:,:); 
                    imgG = img(1:(size(img,2)/2),:,:);
                end
        end
    end

    %Function 2: Let user keep changing background threshhold for
    %spotfinding till satisfied:
    function newspots = SptFindUserThresh(oldspots,SpotImg,thisn,thisxout,...
            ChName,OrdfiltSize,MinDist)
        newspots = oldspots;
        happy = 0;
        while happy==0
            answer = input('Press enter if satisfied, anything else if not:','s');
            if ~isempty(answer)
                close
                bar(thisxout,thisn)
                title('Choose a threshold between background and true spots:','Fontsize',12)
                xlabel('<-Background intensities ... Real spot intensities->')
                ylabel('Counts')
                newthresh = input('Enter new threshold to use:');
                close
                %Error handling:
                while newthresh <=0 || newthresh >= 1
                    newthresh = input('Enter new threshold to use:');
                end
                [newspots,thisn,thisxout] = FindSpotsV5(SpotImg,'ShowResults',1,...
                    'UserThresh',newthresh,'ImgTitle',ChName,'NeighborhoodSize',OrdfiltSize,...
                    'maxsize',MinDist);
            else
                happy = 1;
            end
        end
        close
    end

    smFRETsetup;
    params = load('AnalysisParameters.mat');

%%%%%%FIRST PART: Channel mapping:
    %Load an old channel mapping, or perform a new one:
    DoMap = input('Press enter to perform channel mapping, anything else to load an old one:','s');

    if ~isempty(DoMap)
        %Default to most recent map:
        prevmapdir = load('PathToRecentMap.mat');
        D_Beads = uigetdir(prevmapdir.MostRecentMapDir,'Select directory with old map');
        if exist(fullfile(D_Beads,'ChannelMapping.mat'),'file')
            Map = load(fullfile(D_Beads,'ChannelMapping.mat'));
            A = Map.A;
            b = Map.b;
            clear Map prevmapdir
        else
            disp(strcat('Bead map not found in',Beads))
            return
        end
    else
        D_Beads = uigetdir(params.defaultdatadir,'Select directory with beads');
        %To easily load the most recent bead map:
        MostRecentMapDir = D_Beads;
        %Figure out how many bead files to analyze there are:
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
            BeadFilesInMap{i} = fullfile(D_Beads,AllBeads(i).name); %Keeps a record of which bead files went into the map
            TotImg = LoadUManagerTifsV5(fullfile(D_Beads,AllBeads(i).name),[1 params.FramesToAvg]);
            if size(TotImg,3) == 1
                allBdImgs(:,:,i) = TotImg;
            else
                %allBdImgs(:,:,i) = mean(TotImg(:,:,1:10),3);
                %With the new version of LoadUManagerTifs:
                allBdImgs(:,:,i) = mean(TotImg,3);
            end

            disp(strcat('Analyzing beads: ',int2str(i)',' of ',int2str(BdDir)))

            %Step 1: Find spots in red and green channels separately, so split the
            %image up into the two channels:
            [imgRed,imgGreen] = SplitImg(TotImg,params);

            %Find spots in green channel
            [spotsG{i},n,xout] = FindSpotsV5(imgGreen,'ShowResults',1,'ImgTitle','Green Channel',...
                'NeighborhoodSize',params.BeadNeighborhood,'maxsize',params.BeadSize);
            spotsG{i} = SptFindUserThresh(spotsG{i},imgGreen,n,xout,'Green Channel',...
                params.BeadNeighborhood,params.BeadSize);
            clear n xout
            
            %Find spots in red channel
            [spotsR{i},n,xout] = FindSpotsV5(imgRed,'ShowResults',1,'ImgTitle','Red Channel',...
                'NeighborhoodSize',params.BeadNeighborhood,'maxsize',params.BeadSize);
            spotsR{i} = SptFindUserThresh(spotsR{i},imgRed,n,xout,'Red Channel',...
                params.BeadNeighborhood,params.BeadSize);
            clear n xout

            if debug %Figure with all the spots found:
                %plot all the boxes for both channels on a big image:
                spotsG_abs(:,1) = spotsG{i}(:,1);
                spotsG_abs(:,2) = spotsG{i}(:,2)+size(TotImg,2)/2;
                if size(TotImg,3)==1
                    PutBoxesOnImageV4(mat2gray(TotImg),[spotsR{i};spotsG_abs],params.BeadSize,'0','w');
                else
                    PutBoxesOnImageV4(mat2gray(mean(TotImg(:,:,1:params.FramesToAvg),3)),[spotsR{i};spotsG_abs],params.BeadSize,'0','w');
                end
                clear spotsG_abs
            end

            %Step 2: figure out which spots in one channel go with the spots in
            %the other channel.  For our beads, which are brighter in the donor
            %than acceptor channel, the matching works best if you match spots
            %in donor to spots in acceptor:
            [matchG{i},matchR{i}] = FindSpotMatches(spotsG{i},spotsR{i});
            matchGall = [matchGall, matchG{i}];
            matchRall = [matchRall, matchR{i}];

            if debug
                %Box in green the ones that were matched:
                matchG_abs(1,:) = matchG{i}(1,:);
                matchG_abs(2,:) = matchG{i}(2,:)+size(TotImg,2)/2;
                figure
                if size(TotImg,3) == 1
                    PutBoxesOnImageV4(mat2gray(TotImg),[matchR{i}';matchG_abs'],params.BeadSize);
                else
                    PutBoxesOnImageV4(mat2gray(mean(TotImg(:,:,1:params.FramesToAvg),3)),[matchR{i}';matchG_abs'],params.BeadSize);
                end
                %Another way of plotting the matching: blue line between points
                %in the two channels, with a green dot for where the point is
                %in the green channel, and the end of the blue line where the
                %point it got matched to in the red channel is.  This is also a
                %great way of looking at the distortion between the two
                %channels
                figure('Position',[200 200 325 625])
                plot([matchR{i}(2,:);matchG{i}(2,:)],0-[matchR{i}(1,:);matchG{i}(1,:)],'-b')
                hold on
                plot(matchR{i}(2,:),0-matchR{i}(1,:),'xr')
                ylim([-512 0])
                xlim([0 256])
                title('Red x is center of point in red channel')
                pause
                close all
                clear matchG_abs
            end

            clear TotImg imgGreen imgRed spotsGabs
        end

        %Step three: calculate the transformation using all pairs of beads,
        %from all bead movies or snapshots that were loaded:
        %Calculate transformation:
        %Update 1/2014: the built-in Matlab function fitgeotrans does a bit
        %better than my hand-written code in CalcChannelMapping, so I
        %switched to using that:
        %[A,b] = CalcChannelMapping(matchGall,matchRall)
        tform = fitgeotrans(matchRall',matchGall','Affine'); %Note different input order for fitgeotrans
        A = tform.T(1:2,1:2)
        b = transpose(-tform.T(3,1:2))
        
        %Plot the results for each movie:
        for i = 1:BdDir
            disp(strcat('Iterating through bead images for user to check quality (',...
                int2str(i),' of ',int2str(BdDir),')'))
            matchG_abs(1,:) = matchG{i}(1,:);
            matchG_abs(2,:) = matchG{i}(2,:)+size(allBdImgs(:,:,i),2)/2;
            newR = A*matchG{i}+repmat(b,1,size(matchG{i},2));
            figure
            PutBoxesOnImageV4(mat2gray(allBdImgs(:,:,i)),[newR';matchG_abs'],params.BeadSize);
            title('Spots found in green, matched to red','Fontsize',12)
            figure
            errs = FindSpotDists(matchR{i},newR);
            hist(min(errs,[],2),0:0.1:10)
            ylabel('Counts','Fontsize',12)
            xlabel('Distance between mapped red bead and real red bead','Fontsize',12)
            CalcCombinedImage(A,b,allBdImgs(:,257:end,i),allBdImgs(:,1:256,i),1);
            pause
            close
            close
            close
            clear newR matchG_abs
        end
        clear allBdImgs matchG matchR matchGall matchRall

        save(fullfile(D_Beads,'ChannelMapping.mat'),'A','b','BeadFilesInMap');
        save('PathToRecentMap','MostRecentMapDir');
    end

%%%%%%SECOND PART: Analyze data:
    D_Data = uigetdir(D_Beads,'Select directory with data to analyze');
    %Figure out how many bead files to analyze there are:
    ToAnalyze = dir(fullfile(D_Data,strcat(rootname,'_*')));
    %Get framerate for plotting:
    if isempty(ToAnalyze)
        disp('Did not find data to analyze; remember not to include _<number> at the end of rootname!') %error handling
        return
    end
    fps = GetInfoFromMetaData(fullfile(D_Data,ToAnalyze(1).name),'fps');
    fps = 1/fps; %This is actually frames per ms
    
    %Make sure not saving over old data:
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
        %Update 12/2013: If this movie has already been analyzed, provide
        %the option to use the previously found spots, instead of
        %re-finding them
        
        useoldspots = 'n';
        
        if exist(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),'file')
            useoldspots = input('Load previous analysis? (y/n)','s');
        end
            
        if strcmpi(useoldspots,'n')
           %Load this movie--just the first 10 frames for finding spots
           TotImg = LoadUManagerTifsV5(fullfile(D_Data,ToAnalyze(i).name),[1 params.FramesToAvg]);
           [imgRed,imgGreen] = SplitImg(TotImg,params);
           
           %Find spots in both channels, but don't double-count:
           %Update 1/2014: finding spots in a composite image, so that
           %mid-FRET spots don't get lost.  NOTE that the combined image
           %will have a frame of reference of the acceptor image, which
           %is fine because that's what I pass into UserSpotSelectionV4.
           
           composite = CalcCombinedImage(A,b,mean(imgGreen,3),...
               mean(imgRed,3));
           
           %Find spots in this new image:
           [spotsR,n,xout] = FindSpotsV5(composite,'ShowResults',1,'ImgTitle','Composite Image',...
                 'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize);
           figure, imshow(mean(imgGreen,3),[]);
           title('Green Channel')
           figure, imshow(mean(imgRed,3),[]);
           title('Red Channel')
            spotsR = SptFindUserThresh(spotsR,composite,n,xout,'Composite Image',...
                params.DNANeighborhood,params.DNASize);
            clear n xout
            
            close all

%             %Find spots in green channel
%             [spotsG,n,xout] = FindSpotsV5(imgGreen,'ShowResults',1,'ImgTitle','Green Channel',...
%                  'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize);
%             spotsG = SptFindUserThresh(spotsG,imgGreen,n,xout,'Green Channel',params.DNANeighborhood,params.DNASize);
%             clear n xout
%             spotsG = spotsG';
%             %Figure out where these spots would be in the red channel:
%             spotsRguess = A*spotsG+repmat(b,1,size(spotsG,2));
%             spotsG_abs(1,:) = spotsG(1,:);
%             spotsG_abs(2,:) = spotsG(2,:)+size(TotImg,2)/2;
%             if debug
%                 figure
%                 PutBoxesOnImageV4(mat2gray(mean(TotImg,3)),[spotsG_abs';spotsRguess'],params.DNASize);
%                 title('Spots found in green, with red pairs','Fontsize',12)
%                 pause
%                 close
%             end
% 
%             %Find spots in red channel
%             [spotsR,n,xout] = FindSpotsV5(imgRed,'ShowResults',1,'ImgTitle','Red Channel',...
%                 'NeighborhoodSize',params.DNANeighborhood,'maxsize',params.DNASize);
%             spotsR = SptFindUserThresh(spotsR,imgRed,n,xout,'Red Channel',params.DNANeighborhood,params.DNASize);
%             clear n xout
%             spotsR = spotsR';
%             %Make sure that none of spotsRguess duplicate spotsR: if any do, keep
%             %only the spotsR version because directly finding the spot is probably
%             %more accurate than using the transformation from the green channel to
%             %guess where it is.  This will also elimninate any spots that are too
%             %close together:
%             SelfDists = FindSpotDists([spotsRguess,spotsR]);
%             spottooclose = SelfDists>params.DNASize;
%             %Each row will be all 1's if the spot represented by the row is more than
%             %sqrt(2)*sqrt(maxsize) away from another spot.  If there's another spot too
%             %close, one or more elements will be zero.  So below, ask if the sum of a
%             %peak's row is equal to the length of the row (minus one, because of the
%             %zero element where it's too close to itself)
%             spotstemp = [];
%             for j = 1:size(spotsRguess,2)
%                 if sum(spottooclose(j,:))==length(spottooclose(j,:))-1
%                     spotstemp(:,end+1) = spotsRguess(:,j);
%                 end
%             end
%             spots = [spotstemp,spotsR];
%             clear spotstemp
% 
%             if debug
%                 figure
%                 PutBoxesOnImageV4(mat2gray(mean(TotImg,3)),[spotsG_abs';spotsRguess';spotsR'],params.DNASize);
%                 title('All spots found in red and green channels, and green matched to red','Fontsize',12)
%                 figure 
%                 PutBoxesOnImageV4(mat2gray(mean(TotImg,3)),spots',params.DNASize);
%                 title('All spots to be analyzed','Fontsize',12)
%                 pause
%                 close all
%             end
            
           save(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),'spots')

           %Iterate through spots; display traces and allow user to select which ones to keep
           disp(strcat('Movie ',int2str(i),'of',int2str(length(ToAnalyze))))
           %Mainly for debugging for now:
           %PutBoxesOnImageV3(mat2gray(mean(TotImg(:,:,1:20),3)),spots',params.SpotSize);
           %title('All spots to be analyzed','Fontsize',12)
           %UserSpotSelectionV3(spots,imgRed,imgGreen,params,A,b,savedir,fps,i);
           UserSpotSelectionV4(spots,fullfile(D_Data,ToAnalyze(i).name),params,A,b,savedir,fps,i);
        else
           oldspots = load(fullfile(savedir,strcat('SpotsFound',int2str(i),'.mat')),'spots');
           disp(strcat('Movie ',int2str(i),'of',int2str(length(ToAnalyze))))
           %Mainly for debugging for now:
           %PutBoxesOnImageV3(mat2gray(mean(TotImg(:,:,1:20),3)),spots',params.SpotSize);
           %title('All spots to be analyzed','Fontsize',12)
           %UserSpotSelectionV3(spots,imgRed,imgGreen,params,A,b,savedir,fps,i);
           UserSpotSelectionV4(oldspots.spots,fullfile(D_Data,ToAnalyze(i).name),params,A,b,savedir,fps,i);
        end
        clear TotImg spots imgRed imgGreen spotsG spotsR spotsG_abs spotsRguess spotstemp
    end
end
