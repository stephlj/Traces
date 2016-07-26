% function UserSpotSelectionV4(allRedI,allGrI,spots,spotVars,PathToMovie,...
%   params,tform,savedir,fps,setnum)
%
% Iterates through all the spots in a movie and allows user to adjust
% background, keep the ones they like, etc. Note that spots are passed in in
% the frame of reference of the ACCEPTOR channel.
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

function UserSpotSelectionV4(allRedI,allGrI,spots,spotVars,PathToMovie,params,tform,savedir,fps,setnum)

zoomsize = 10; % How many pixels across (and high) should the zoomed-in image around
    % a spot be
if params.DNASize>=zoomsize
    zoomsize = params.DNASize+1;
end

%%%Setting up some stuff
% subfunction for putting circles around a spot:
    function boxfun(currspot,circlesize,markercolor)
        % CalcSpotIntensityNoGauss puts a circle of diameter 5 around each spot:
        t = 0:pi/100:2*pi;
        plot(currspot(2)+circlesize(2)/2.*cos(t),currspot(1)+circlesize(1)/2.*sin(t),strcat('-',markercolor))
        clear t
    end

% Make spots 2-by-numspots matrices, if it's not already
if size(spots,1)~=2
    spots = transpose(spots);
end

% Finding an injection point to plot, if any:
ScalingParams = load(fullfile(PathToMovie,'ScalingInfo.mat'));
if isfield(ScalingParams,'InjectPoints')
   InjMarkToPlot = ScalingParams.InjectPoints+params.InjectDelay;
% Backwards compatibility
elseif isfield(params,'InjectPoints')
    InjMarkToPlot = params.InjectPoints+params.InjectDelay;
    InjectPoints = params.InjectPoints;
    save(fullfile(PathToMovie,'ScalingInfo.mat','InjectPoints','-append'));
    clear InjectPoints
elseif isfield(params,'ManualInjectMark')
   InjMarkToPlot = params.ManualInjectMark;
else
   InjMarkToPlot = [];
end
t_Inj = InjMarkToPlot;

% Find spots in green channel:
SpotsAndIstruct = load(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')));
if isfield(SpotsAndIstruct,'SpotsInG')
    GrSpots = SpotsAndIstruct.SpotsInG;
else
    GrSpots = tform.FRETmapInv(spots);
end
% Update 1/2015: saving these following variables (and loading them in if
% they already exist) so that you don't have to redo if you re-load an old 
% data set:
if isfield(SpotsAndIstruct,'xlims')
    Rbkgnd = SpotsAndIstruct.Rbkgnd;
    Gbkgnd = SpotsAndIstruct.Gbkgnd;
    xlims = SpotsAndIstruct.xlims;
    ends = SpotsAndIstruct.ends;
    offset = SpotsAndIstruct.offset;
else
    Rbkgnd = zeros(size(spots,2),1);
    Gbkgnd = zeros(size(spots,2),1);
    xlims(:,1) = zeros(size(spots,2),1);
    xlims(:,2) = ones(size(spots,2),1).*size(allRedI,2)/fps;
    ends = zeros(size(spots,2),1); % Where the end of the FRET signal should be (zero after this point)
    offset = 1;
    % Save these fields to the SpotsAndIntensities file (user-indicated
    % modifications will get re-saved below)--note that '-append' will
    % either add a field, if it doesn't exist already, or replace the field
    % with the currently saved value
    save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),...
        'Rbkgnd','Gbkgnd','xlims','ends','offset','-append')
end

% Update 3/2015: Enabling all the re-mapping to happen at the end, and then
% the user goes through only those re-mapped spots.  This vector will index
% those re-mapped spots and m will index this vector:
kk = [];
% Then these vectors will keep track of which spots need to be recalculated
% at the end. Load them from disk if the user's analysis was interrupted in
% the middle:
if isfield(SpotsAndIstruct,'RedToReCalc')
    RedToReCalc = SpotsAndIstruct.RedToReCalc;
    GrToReCalc = SpotsAndIstruct.GrToReCalc;
    AllToReCalc = SpotsAndIstruct.AllToReCalc;
else
    RedToReCalc = [];
    GrToReCalc = [];
    % and this is a placeholder that will be turned into kk to indicate it's
    % time to go through the re-mapped spots:
    AllToReCalc = [];
end

clear SpotsAndIstruct

% Get an average image of the first 10 frames to display:
[imgRinit,imgGinit] = LoadScaledMovie(PathToMovie,[1 1+params.FramesToAvg],params);
if imgRinit==-1
    return
end
imgRinitavg = mat2gray(mean(imgRinit,3));
imgGinitavg = mat2gray(mean(imgGinit,3));
imgRinit = imgRinitavg;
imgGinit = imgGinitavg;

%%%Interactive section
k = 1; % Indexes current spot being plotted

h2 = figure('Position',params.Fig2Pos);
h1 = figure('Position',params.Fig1Pos);

disp('Fig. 2 must be current.') 
disp('.=fwd; ,=back; b=background adjust; r=reset background; s=save; z = zoom; u=unzoom;')
disp('f=select frame to display; m = play movie between two points; a=show average around frame;')
disp('o=adjust black offset; l = re-_l_ocate spot now but re-calculate intensities later;')
disp('c = re-locate spot and re-_c_alculate intensities now; g=go to spot number;')
disp('e=end of trace (after this point FRET set to zero); d=done with movie')

    while k <= size(spots,2)  
       % Calculate raw intensities and raw FRET, then correct for gamma
       % and alpha as user wants:
       % Note: Before smoothing, even if user wants to smooth intensities and/or 
       % FRET, always save the raw as well (see below)
       rawRedI = allRedI(k,:)-Rbkgnd(k);
       rawGrI = allGrI(k,:)-Gbkgnd(k);
       % Correct for gamma and alpha as user wants:
       unsmoothedRedI = rawRedI - params.alpha*rawGrI;
       unsmoothedGrI = rawGrI;
       unsmoothedFRET = unsmoothedRedI./(unsmoothedGrI+params.gamma*unsmoothedRedI);% Note the raw FRET will be corrected for alpha and gamma, but
        % because the raw intensities are also saved, you can always
        % recalculate FRET from those.  rawFRET just means not smoothed.
       if params.SmoothIntensities>0
           % Update 9/2014: using a median filter if possible instead
           if exist('medfilt1') % Part of the signal processing toolbox
               RedI = medfilt1(unsmoothedRedI,round(params.SmoothIntensities));
               GrI = medfilt1(unsmoothedGrI,round(params.SmoothIntensities));
           elseif exist('medfilt2') % Part of the image processing toolbox, which is required for 
                % other parts of UserSpotsSelection, but just in case ... 
               RedI = medfilt2(unsmoothedRedI,[1,round(params.SmoothIntensities)]);
               GrI = medfilt2(unsmoothedGrI,[1,round(params.SmoothIntensities)]);
           elseif exist('smooth') % Part of the curve fitting toolbox
               RedI = smooth(unsmoothedRedI,round(params.SmoothIntensities)); % Moving average with span = SmoothIntensities
               GrI = smooth(unsmoothedGrI,round(params.SmoothIntensities));
           else
               disp('Cannot smooth data with current toolboxes; displaying unsmoothed.')
               RedI = unsmoothedRedI;
               GrI = unsmoothedGrI;
           end
       else
           RedI = unsmoothedRedI;
           GrI = unsmoothedGrI;
       end
       if params.SmoothFRET>0
           if exist('medfilt1')
               FRET = medfilt1(unsmoothedFRET,round(params.SmoothFRET));
           elseif exist('medfilt2')
               FRET = medfilt2(unsmoothedFRET,[1,round(params.SmoothFRET)]);
           elseif exist('smooth')
               FRET = smooth(unsmoothedFRET,round(params.SmoothFRET));
           else
               disp('Cannot smooth data with current toolboxes; displaying unsmoothed.')
               FRET = unsmoothedFRET;
           end
       else
           FRET = unsmoothedFRET;
       end
       
       if ends(k)~=0
           FRET(ends(k):end) = 0;
       end
       % TODO: convolve with a Gaussian instead of using a moving average
        
       xvect = ((1:length(RedI))./fps);
       
       [imgRzoom,zoomcenR] = ExtractROI(imgRinit,zoomsize,spots(:,k));
       [imgGzoom,zoomcenG] = ExtractROI(imgGinit,zoomsize,GrSpots(:,k));
       
       % In case ExtractROI had to remove some pixels because the ROI
       % center was too close to the edge:
       zoomsizeR = size(imgRzoom,1)-1;
       zoomsizeG = size(imgGzoom,1)-1;
       
       % Show plots
       figure(h2)
       % plot red channel with circle around spot
       %subplot('Position',[0.05 0.3 0.45 0.65])
       subplot('Position',[0.08 0.23 0.39 0.39*size(imgRinit,1)/size(imgRinit,2)])
       imshow(imgRinit)
       hold on
       boxfun(spots(:,k),sqrt(1./(2.*spotVars(:,k))).*5,'r');
       hold off
       title('Red','Fontsize',12)
       % plot green channel with circle around spot
       subplot('Position',[0.54 0.23 0.39 0.39*size(imgRinit,1)/size(imgRinit,2)])
       imshow(imgGinit)
       hold on
       boxfun(GrSpots(:,k),sqrt(1./(2.*spotVars(:,k))).*5,'g');
       hold off
       title('Green','Fontsize',12)
       % Show zoom in of the red spot, with the fitted Gaussian or a circle
       % around it to show how and where the spot intensity was calculated
       imgRzoom_axes = subplot('Position',[0.13 0.05 0.2 .18]);
       imshow(imgRzoom);
       hold on
       boxfun(zoomcenR,sqrt(1./(2.*spotVars(:,k))).*5,'r');
       hold off
       % Same for green
       imgGzoom_axes = subplot('Position',[0.63 0.05 0.2 .18]);
       imshow(imgGzoom)
       hold on
       boxfun(zoomcenG,sqrt(1./(2.*spotVars(:,k))).*5,'g');
       hold off
       
       figure(h1)
       trace_axes = subplot(2,1,1);
       plot(xvect,RedI,'-r',xvect,GrI,'-g',xvect,RedI+GrI+offset,'-k',...
           [xvect(1) xvect(end)],[0 0],'--k')
       if ~isempty(InjMarkToPlot)
           hold on
           for ss = 1:length(InjMarkToPlot)
               plot([InjMarkToPlot(ss), InjMarkToPlot(ss)], ...
                   [0 max(RedI+GrI+offset)+0.1],'-m')
           end
           hold off
       end
       xlabel('Time (sec)','Fontsize',12)
       ylabel('Intensity (a.u.)','Fontsize',12)
       title(strcat('Spot',int2str(k),'/',int2str(size(spots,2))),'Fontsize',12)
       xlim([xlims(k,1) xlims(k,2)])
       
       fret_axes = subplot(2,1,2);
       plot(xvect,FRET,'-b',[xvect(1) xvect(end)],[0 0],'--k',...
           [xvect(1) xvect(end)],[1 1],'--k')
       if ~isempty(InjMarkToPlot)
           hold on
           for ss = 1:length(InjMarkToPlot)
               plot([InjMarkToPlot(ss), InjMarkToPlot(ss)],[-0.2 1.2],'-m')
           end
           hold off
       end
       xlabel('Time (sec)','Fontsize',12)
       ylabel('FRET','Fontsize',12)
       xlim([xlims(k,1) xlims(k,2)])
       ylim([-0.2 1.2])
       
       % Interactive section:
       cc=1;
        while cc~=13
            ct=waitforbuttonpress;
            cc=get(gcf,'currentcharacter');
            
            if ct==1
                if ~isempty(cc)
                    % Go forward to the next spot
                    if cc=='.'
                        % If you've not gotten through the whole movie once:
                        if isempty(kk)
                            k = k+1;
                            imgRinit = imgRinitavg;
                            imgGinit = imgGinitavg;
                            % If you have now gone through the movie once and have
                            % things to re-map, re-map them now:
                            if k>size(spots,2) && sum(AllToReCalc)>0
                                RedToReCalc = sort(RedToReCalc);
                                GrToReCalc = sort(GrToReCalc);
                                if ~isempty(RedToReCalc)
                                    disp('Calculating intensities for re-mapped red spots ...')
                                    [allRedI(RedToReCalc,:), ~] = CalcIntensitiesV3(PathToMovie,...
                                        spots(:,RedToReCalc), spotVars(:,RedToReCalc),-1,params);
                                    RedI = allRedI;
                                    % Re-set RedToReCalc to empty to indicate
                                    % we've re-calculated everything here
                                    RedToReCalc = [];
                                    save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),...
                                        'RedI','-append')
                                end
                                if ~isempty(GrToReCalc)
                                    disp('Calculating intensities for re-mapped green spots ...')
                                    [~, allGrI(GrToReCalc,:)] = CalcIntensitiesV3(PathToMovie,...
                                        GrSpots(:,GrToReCalc), spotVars(:,GrToReCalc),1,params);
                                    GrI = allGrI;
                                    GrToReCalc = [];
                                    save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),...
                                        'GrI','-append')
                                end
                                clear RedI GrI
                                % Fill kk appropriately
                                kk = sort(AllToReCalc);
                                % then empty it and save that to disk
                                AllToReCalc = [];
                                save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),...
                                        'RedToReCalc','GrToReCalc','AllToReCalc','-append')
                                m = 1;
                                k = kk(m);
                            end
                        elseif m==length(kk)
                            k = size(spots,2)+1;
                            % (This will break us out of the while loop and end this analysis session)
                        else
                            m = m+1;
                            k = kk(m);
                            imgRinit = imgRinitavg;
                            imgGinit = imgGinitavg;
                        end
                        cc=13;
                    % Go back one spot
                    elseif cc==',' 
                        if isempty(kk) && k>1
                            k=k-1;
                            imgRinit = imgRinitavg;
                            imgGinit = imgGinitavg;
                        elseif ~isempty(kk) && m>1
                            m = m-1;
                            k = kk(m);
                            imgRinit = imgRinitavg;
                            imgGinit = imgGinitavg;
                        end
                        cc=13;
                    % Go to specific spot
                    elseif cc=='g'
                        newk = input('Go to bead number: ');
                        % Update 3/2015: Only allow this to work if the spot
                        % the user wants to look at is one that was re-mapped:
                        if isempty(kk)
                            if newk<=size(spots,2) && newk>=1
                                k=newk;
                                imgRinit = imgRinitavg;
                                imgGinit = imgGinitavg;
                            end
                        elseif ~isempty(kk) && ismember(newk,kk)
                            [~,indx] = ismember(newk,kk);
                            m=indx;
                            k = kk(m);
                            imgRinit = imgRinitavg;
                            imgGinit = imgGinitavg;
                        end
                        cc=13;
                    % go to the next movie:
                    elseif cc=='d'
                        k = size(spots,2)+1;
                        cc=13;
                    % Set background levels
                    elseif cc=='b'
                        [x,~] = ginput(2);
                        % Make sure user actually selected something in the
                        % correct panel
                        if isequal(trace_axes,gca)
                            x = sort(x);
                            Rbkgnd(k) = mean(RedI(round(x(1)*fps:x(2)*fps)));
                            Gbkgnd(k) = mean(GrI(round(x(1)*fps:x(2)*fps)));
                            clear x
                            save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),...
                                'Rbkgnd','Gbkgnd','-append')
                        end
                        cc=13;
                    % Reset background to zero
                    elseif cc=='r'
                        Rbkgnd(k)=0;
                        Gbkgnd(k)=0;
                        save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),...
                                'Rbkgnd','Gbkgnd','-append')
                        cc=13;
                    % Save figure.  Saves both the figure and a .mat of
                    % the red and green intensities and FRET values.
                    elseif cc=='s'
                        saveas(gca,fullfile(savedir,strcat('Spot',int2str(setnum),'_',int2str(k))),'fig')
                        print('-depsc',fullfile(savedir,strcat('Spot',int2str(setnum),'_',int2str(k))))
                        if xlims(k,2)~=size(allRedI,2)/fps
                            RedToSave = RedI(round(xlims(k,1)*fps:xlims(k,2)*fps));
                            GrToSave = GrI(round(xlims(k,1)*fps:xlims(k,2)*fps));
                            FRETtoSave = FRET(round(xlims(k,1)*fps:xlims(k,2)*fps));
                            rawRedToSave = rawRedI(round(xlims(k,1)*fps:xlims(k,2)*fps));
                            rawGrToSave = rawGrI(round(xlims(k,1)*fps:xlims(k,2)*fps));
                            rawFRETtoSave = unsmoothedFRET(round(xlims(k,1)*fps:xlims(k,2)*fps));
                            unsmthRedToSave = unsmoothedRedI(round(xlims(k,1)*fps:xlims(k,2)*fps));
                            unsmthGrToSave = unsmoothedGrI(round(xlims(k,1)*fps:xlims(k,2)*fps));
                        else
                            RedToSave = RedI;
                            GrToSave = GrI;
                            FRETtoSave = FRET;
                            rawRedToSave = rawRedI;
                            rawGrToSave = rawGrI;
                            rawFRETtoSave = unsmoothedFRET;
                            unsmthRedToSave = unsmoothedRedI;
                            unsmthGrToSave = unsmoothedGrI;
                        end 
                        clear RedI GrI FRET rawRedI rawGrI unsmoothedFRET unsmoothedRedI unsmoothedGrI
                        RedI = RedToSave;
                        GrI = GrToSave;
                        FRET = FRETtoSave;
                        rawRedI = rawRedToSave;
                        rawGrI = rawGrToSave;
                        unsmoothedRedI = unsmthRedToSave;
                        unsmoothedGrI = unsmthGrToSave;
                        unsmoothedFRET = rawFRETtoSave;
                        Rspot = spots(:,k);
                        Gspot = GrSpots(:,k);
                        variance = spotVars(:,k);
                        save(fullfile(savedir,strcat('Spot',int2str(setnum),'_',int2str(k),'.mat')),...
                            'RedI','GrI','FRET','rawRedI','rawGrI','unsmoothedRedI','unsmoothedGrI',...
                            'unsmoothedFRET','fps','Rspot','Gspot','variance','t_Inj')
                        clear RedToSave GrToSave FRETtoSave rawRedToSave rawGrToSave rawFRETtoSave
                        clear unsmthRedToSave unsmthGrToSave
                        clear Rspot Gspot variance
                        cc=13;
                    % Zoom
                    elseif cc=='z'
                        [x,~] = ginput(2);
                        if isequal(trace_axes,gca) && x(1)~=x(2)
                            x = sort(x);
                            if x(1)<0
                                x(1)=1;
                            end
                            if x(2)>size(allRedI,2)
                                x(2)=size(allRedI,2);
                            end
                            xlims(k,1) = x(1);
                            xlims(k,2) = x(2);
                        end
                        save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),...
                                'xlims','-append')
                        cc=13;
                    % Unzoom
                    elseif cc=='u'
                        xlims(k,:) = [0,size(allRedI,2)/fps];
                        save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),...
                                'xlims','-append')
                        cc=13;
                    elseif cc=='o'
                        offset = input('New offset: ');
                        save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),...
                                'offset','-append')
                        cc=13;
                    % Set end of signal (FRET set to zero after this point)
                    elseif cc=='e'
                        if isequal(trace_axes,gca) || isequal(fret_axes,gca)
                            [x,~] = ginput(1);
                            ends(k) = round(x*fps);
                        end
                        save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),...
                                'ends','-append')
                        cc=13;
                    %Show a specific frame in figure 2
                    elseif cc=='f'
                        [x,~] = ginput(1);
                        if isequal(trace_axes,gca) || isequal(fret_axes,gca)
                            % x will be in seconds, not frames. Convert to frame:
                            x = x*fps;
                            [imgRinit,imgGinit] = PlayMovie(PathToMovie,[round(x) round(x)],params,...
                                h2,strcat('subplot(',char(39),'Position',char(39),...
                                    ',[0.08 0.23 0.39 0.39*',int2str(size(imgRinit,1)/size(imgRinit,2)),'])'),...
                                strcat('subplot(',char(39),'Position',char(39),...
                                    ',[0.54 0.23 0.39 0.39*',int2str(size(imgRinit,1)/size(imgRinit,2)),'])'),...
                                spots(:,k),GrSpots(:,k),...
                                strcat('subplot(',char(39),'Position',char(39),',[0.13 0.05 0.2 .18])'),...
                                strcat('subplot(',char(39),'Position',char(39),',[0.63 0.05 0.2 .18])'),...
                                zoomsize,sqrt(1./(2.*spotVars(:,k))).*5);
                            clear x
                        end
                        cc = 13;
                    %Play a section of the movie in figure 2
                    elseif cc=='m'
                        [x,~] = ginput(2);
                        if isequal(trace_axes,gca) || isequal(fret_axes,gca)
                            % x will be in seconds, not frames. Convert to frame:
                            x = x*fps;
                            x = round(sort(x));
                            [imgRinit,imgGinit] = PlayMovie(PathToMovie,[x(1) x(2)],params,...
                                h2,strcat('subplot(',char(39),'Position',char(39),...
                                    ',[0.08 0.23 0.39 0.39*',int2str(size(imgRinit,1)/size(imgRinit,2)),'])'),...
                                strcat('subplot(',char(39),'Position',char(39),...
                                    ',[0.54 0.23 0.39 0.39*',int2str(size(imgRinit,1)/size(imgRinit,2)),'])'),...
                                spots(:,k),GrSpots(:,k),...
                                strcat('subplot(',char(39),'Position',char(39),',[0.13 0.05 0.2 .18])'),...
                                strcat('subplot(',char(39),'Position',char(39),',[0.63 0.05 0.2 .18])'),...
                                zoomsize,sqrt(1./(2.*spotVars(:,k))).*5);
                            clear x
                        end
                        cc = 13;
                    % Show an average of 10 frames around where the user clicks
                    elseif cc=='a'
                        [x,~] = ginput(1);
                        if isequal(trace_axes,gca) || isequal(fret_axes,gca)
                            % x will be in seconds, not frames. Convert to frame:
                            x = x*fps;
                            [imgRinit,imgGinit] = LoadScaledMovie(PathToMovie,...
                                [round(x)-ceil(params.FramesToAvg/2) round(x)+ceil(params.FramesToAvg/2)],...
                                params);
                            imgRinit = mat2gray(mean(imgRinit,3));
                            imgGinit = mat2gray(mean(imgGinit,3));
                            clear x
                        end
                        cc = 13;
                    % Re-locate a spot in one of the channels, if the
                    % transformation looks like it didn't get the center right:
                    % Update 3/2015: Add the option to re-find the spot as
                    % usual, but do all the re-calculating of intensities all
                    % at once at the end:
                    elseif cc=='l' || cc=='c'
                        disp('Click on trace at time around which to look for spot:')
                        [xT,~] = ginput(1);
                        if isequal(trace_axes,gca) || isequal(fret_axes,gca)
                            xT = round(xT*fps);
                            starttime = xT-ceil(params.FramesToAvg/2);
                            endtime = xT+ceil(params.FramesToAvg/2);
                            if starttime<1
                                endtime = endtime+0-starttime;
                                starttime = 1;
                            end
                            if endtime>size(allRedI,2)
                                starttime = starttime+(endtime-allRedI);
                                endtime = size(allRedI,2);
                            end
                            spottorefit = input('Relocate red spot (r) or green spot (g)?: ','s');
                            disp('Click on zoomed image of the spot where you want to look for a new one:')
                            figure(h2)
                            [xIlocal,yIlocal] = ginput(1);
                            if strcmpi(spottorefit,'r') && isequal(imgRzoom_axes,gca)
                                InitSpot = spots(:,k);
                                newcoords = GlobalToROICoords([],[yIlocal;xIlocal],spots(:,k),zoomsizeR,zoomsizeR);
                                [imgs,~] = LoadScaledMovie(PathToMovie,[starttime endtime],params);
                                [tempnewspot, ~] = FindRefinedSpotCenters(imgs,newcoords,0.02,params);
                                % Check the new spot isn't too close to the
                                % boundary
                                if ~isempty(tempnewspot)
                                    [tempnewspot,~,~,~] = CheckSpotBoundaries(tempnewspot,...
                                            [],[],[],params,PathToMovie);
                                    if ~isempty(tempnewspot)
                                        spots(:,k) = tempnewspot;
                                        NewSpot = spots(:,k);
                                        disp(sprintf('Change in spot center = %d pixels',sqrt((NewSpot(1)-InitSpot(1))^2+(NewSpot(2)-InitSpot(2))^2)))
                                        % Update 3/2015: Here's where the change happens. Don't re-calculate
                                        % till later if user entered 'c'. However, disallow the use of 'c' if
                                        % we're already looking through re-calculated spots:
                                        SpotsInR = spots;
                                        if cc=='c' || (cc=='l' && ~isempty(kk))
                                            [allRedI(k,:), ~] = CalcIntensitiesV3(PathToMovie,...
                                                spots(:,k), spotVars(:,k),-1,params);
                                            RedI = allRedI;
                                            save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),'SpotsInR',...
                                                'RedI','-append')
                                        else
                                            % Note that the new spot location is "immortalized" in spots; so we
                                            % just need to know that this spot needs its intensities
                                            % re-calculated, we already have the new spot location saved:
                                            if ~ismember(k,RedToReCalc)
                                                RedToReCalc(end+1) = k;
                                            end
                                            if ~ismember(k,AllToReCalc)
                                                AllToReCalc(end+1) = k;
                                            end
                                            save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),'SpotsInR',...
                                                'GrToReCalc','RedToReCalc','AllToReCalc','-append')
                                        end
                                        clear SpotsInR RedI NewSpot
                                    else
                                        disp('New spot center too close to edge.')
                                    end
                                else
                                    disp('Failed to find new spot center.')
                                end
                            elseif strcmpi(spottorefit,'g') && isequal(imgGzoom_axes,gca)
                                InitSpot = GrSpots(:,k);
                                newcoords = GlobalToROICoords([],[yIlocal;xIlocal],GrSpots(:,k),zoomsizeG,zoomsizeG);
                                [~,imgs] = LoadScaledMovie(PathToMovie,[starttime endtime],params);
                                [tempnewspot, ~] = FindRefinedSpotCenters(imgs,newcoords,0.02,params);
                                if ~isempty(tempnewspot)
                                    [tempnewspot,~,~,~] = CheckSpotBoundaries(tempnewspot,...
                                            [],[],[],params,PathToMovie);
                                    if ~isempty(tempnewspot)
                                        % if the fit fails, tempnewspot will be empty
                                        GrSpots(:,k) = tempnewspot;
                                        NewSpot = GrSpots(:,k);
                                        disp(sprintf('Change in spot center = %d pixels',sqrt((NewSpot(1)-InitSpot(1))^2+(NewSpot(2)-InitSpot(2))^2)))
                                        % Update 3/2015: Here's where the change happens. Don't re-calculate
                                        % till later if user entered 'c':
                                        SpotsInG = GrSpots;
                                        if cc=='c' || (cc=='l' && ~isempty(kk))
                                            [~, allGrI(k,:)] = CalcIntensitiesV3(PathToMovie,...
                                                GrSpots(:,k), spotVars(:,k),1,params);
                                            GrI = allGrI;
                                            save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),'SpotsInG',...
                                                'GrI','-append')
                                        else
                                            % Note that the new spot location is "immortalized" in GrSpots; so we
                                            % just need to know that this spot needs its intensities
                                            % re-calculated, we already have the new spot location saved:
                                            if ~ismember(k,GrToReCalc)
                                                GrToReCalc(end+1) = k;
                                            end
                                            if ~ismember(k,AllToReCalc)
                                                AllToReCalc(end+1) = k;
                                            end
                                            save(fullfile(savedir,strcat('SpotsAndIntensities',int2str(setnum),'.mat')),'SpotsInG',...
                                                'GrToReCalc','RedToReCalc','AllToReCalc','-append')
                                        end
                                        clear SpotsInG GrI NewSpot
                                    else
                                        disp('New spot center too close to edge.')
                                    end
                                else
                                    disp('Failed to find new spot center.')
                                end
                            end
                            clear imgs xT xIlocal yIlocal starttime endtime newcoords
                        end
                        cc=13;
                    % Don't let extra "enters" build up:
                    elseif isequal(cc,char(13)) %13 is the ascii code for the return key
                        cc=13;
                    end
                else
                    cc=13;
                end
            end
       % end interactive section
        end
    % End loop over all spots 
    end
close all

end