% function UserSpotSelectionV4(allRedI,allGrI,spots,PathToMovie,params,tform,savedir,fps,setnum)
%
% Iterates through all the spots in a movie and allows user to adjust
% background, keep the ones they like, etc. Note that spots are passed in in
% the frame of reference of the ACCEPTOR channel.
%
% Updated 12/2013 to allow the user to pass the figure position information
% via params, so it's easier to put this code onto different computers.
%
% Updated 2/2014 so that the intensities are passed in, instead of being
% calculated here.
%
% Steph 10/2013
% Copyright 2013 Stephanie Johnson, University of California, San Francisco

function UserSpotSelectionV4(allRedI,allGrI,spots,PathToMovie,params,tform,savedir,fps,setnum)

zoomsize = 10; % How many pixels across (and high) should the zoomed-in image around
    % a spot be
if params.DNASize>=zoomsize
    zoomsize = params.DNASize+1;
end

%%%Setting up some stuff
% subfunction for putting circles around a spot:
    function boxfun(currspot)
        % CalcSpotIntensityNoGauss puts a circle of diameter 5 around each spot:
        t = 0:pi/100:2*pi;
        plot(currspot(2)+5/2.*cos(t),currspot(1)+5/2.*sin(t),'-g')
        clear t
    end
    % subfunction for finding the local spot center in an ROI
    function localcen = FindLocalCen(ROI,currspot)
        if round(currspot(1))<=currspot(1)
            localcen(1) = size(ROI,2)/2+(currspot(1)-round(currspot(1)));
        else
            localcen(1) = size(ROI,2)/2-(round(currspot(1))-currspot(1));
        end
        if round(currspot(2))<=currspot(2)
            localcen(2) = size(ROI,1)/2+(currspot(2)-round(currspot(2)));
        else
            localcen(2) = size(ROI,1)/2-(round(currspot(2))-currspot(2));
        end
    end

% Make spots 2-by-numspots matrices, if it's not already
if size(spots,1)~=2
    spots = transpose(spots);
end

% Find spots in green channel:
GrSpots = transpose(transformPointsInverse(tform,spots'));

% Get an average image of the first 10 frames to display:
[imgRinit,imgGinit] = LoadScaledMovie(PathToMovie,[1 1+params.FramesToAvg]);
imgRinit=mat2gray(mean(imgRinit,3));
imgGinit=mat2gray(mean(imgGinit,3));

%%%Interactive section
k = 1; % Indexes current spot being plotted

% So that when you go back to a previous spot, you don't have to redo
% selections you did before:
Rbkgnd = zeros(size(spots,2),1);
Gbkgnd = zeros(size(spots,2),1);
xlims = zeros(size(spots,2),2);
ends = zeros(size(spots,2),1); % Where the end of the FRET signal should be (zero after this point)
offset = 0.1;

h2 = figure('Position',params.Fig2Pos);
h1 = figure('Position',params.Fig1Pos);

disp('Fig. 2 must be current.') 
disp('.=fwd; ,=back; b=background adjust; r=reset background; s=save; z = zoom; u=unzoom; o=adjust black offset;')
disp('f=select frame to display; m = play movie between two points;')
disp(' d=done with movie; e=end of trace (after this point FRET set to zero)')

    while k <= size(spots,2)
       RedI = allRedI(k,:)-Rbkgnd(k);
       GrI = allGrI(k,:)-Gbkgnd(k);
       FRET = RedI./(RedI+GrI);
       if params.SmoothIntensities>0
           % Another thing I'm not sure I should do: Smooth the intensities and FRET
           % values (see below)
           SmoothIntensities = round(params.SmoothIntensities); % User error handling
           tempRedI = RedI;
           tempGrI = GrI;
           clear RedI GrI
           RedI = smooth(tempRedI,SmoothIntensities); % Moving average with span = SmoothIntensities
           GrI = smooth(tempGrI,SmoothIntensities);
           clear tempRedI tempGrI SmoothIntensities
       end
       if params.SmoothFRET>0
           tempFRET = FRET;
           SmoothFRET = round(params.SmoothFRET); % User error handling
           clear FRET
           FRET = smooth(tempFRET,SmoothFRET);
           clear tempFRET SmoothFRET
       end
       
       if ends(k)~=0
           FRET(ends(k):end) = 0;
       end
       % TODO: convolve with a Gaussian instead of using a moving average
        
       xvect = ((1:length(RedI))./fps)*10^-3; % fps is actually frames per ms
       
       imgRzoom = ExtractROI(imgRinit,zoomsize,spots(:,k));
       imgGzoom = ExtractROI(imgGinit,zoomsize,GrSpots(:,k));
       
       % Show plots
       figure(h2)
       % plot red channel with circle around spot
       %subplot('Position',[0.05 0.3 0.45 0.65])
       subplot('Position',[0.08 0.23 0.39 0.39*size(imgRinit,1)/size(imgRinit,2)])
       imshow(imgRinit,[])
       hold on
       boxfun(spots(:,k));
       hold off
       title('Red','Fontsize',12)
       % plot green channel with circle around spot
       subplot('Position',[0.54 0.23 0.39 0.39*size(imgRinit,1)/size(imgRinit,2)])
       imshow(imgGinit,[])
       hold on
       boxfun(GrSpots(:,k));
       hold off
       title('Green','Fontsize',12)
       % Show zoom in of the red spot, with the fitted Gaussian or a circle
       % around it to show how and where the spot intensity was calculated
       subplot('Position',[0.13 0.05 0.2 .18])
       imshow(imgRzoom,[])
       hold on
       zoomcenR = FindLocalCen(imgRzoom,spots(:,k));
       boxfun(zoomcenR);
       hold off
       % Same for green
       subplot('Position',[0.63 0.05 0.2 .18])
       imshow(imgGzoom,[])
       hold on
       zoomcenG = FindLocalCen(imgGzoom,GrSpots(:,k));
       boxfun(zoomcenG);
       hold off
       
       figure(h1)
       subplot(2,1,1)
       plot(xvect,RedI,'-r',xvect,GrI,'-g',xvect,RedI+GrI+offset,'-k')
       xlabel('Time (sec)','Fontsize',12)
       ylabel('Intensity (a.u.)','Fontsize',12)
       title(strcat('Spot',int2str(k),'/',int2str(size(spots,2))),'Fontsize',12)
       if xlims(k,1)~=0
           xlim([xlims(k,1) xlims(k,2)])
       end
       
       subplot(2,1,2)
       plot(xvect,FRET,'-b')
       xlabel('Time (sec)','Fontsize',12)
       ylabel('FRET','Fontsize',12)
       if xlims(k,1)~=0
           xlim([xlims(k,1) xlims(k,2)])
       end
       ylim([-.2 1.2])
       
       % Interactive section:
       cc=1;
        while cc~=13
            ct=waitforbuttonpress;
            cc=get(gcf,'currentcharacter');
            
            if ct==1
                % Go forward to the next bead
                if cc=='.'
                    k = k+1;
                    cc = 13;
                % Go back one bead
                elseif cc==',' 
                    if k>1
                        k=k-1;
                    end
                    cc=13;
                % go to the next movie:
                elseif cc=='d'
                    k = size(spots,2)+1;
                    cc=13;
                % Set background levels
                elseif cc=='b'
                    [x,~] = ginput(2);
                    x = sort(x);
                    Rbkgnd(k) = mean(RedI(round(x(1)*fps/10^-3:x(2)*fps/10^-3)));
                    Gbkgnd(k) = mean(GrI(round(x(1)*fps/10^-3:x(2)*fps/10^-3)));
                    clear x
                    cc=13;
                % Reset background to zero
                elseif cc=='r'
                    Rbkgnd(k)=0;
                    Gbkgnd(k)=0;
                    cc=13;
                % Save figure.  Saves both the figure and a .mat of
                % the red and green intensities and FRET values.
                elseif cc=='s'
                    saveas(gca,fullfile(savedir,strcat('Spot',int2str(setnum),'_',int2str(k))),'fig')
                    print('-depsc',fullfile(savedir,strcat('Spot',int2str(setnum),'_',int2str(k))))
                    if xlims(k,1)~=0
                        RedToSave = RedI(round(xlims(k,1)*fps/10^-3:xlims(k,2)*fps/10^-3));
                        GrToSave = GrI(round(xlims(k,1)*fps/10^-3:xlims(k,2)*fps/10^-3));
                        FRETtoSave = FRET(round(xlims(k,1)*fps/10^-3:xlims(k,2)*fps/10^-3));
                    else
                        RedToSave = RedI;
                        GrToSave = GrI;
                        FRETtoSave = FRET;
                    end 
                    clear RedI GrI FRET
                    RedI = RedToSave;
                    GrI = GrToSave;
                    FRET = FRETtoSave;
                    save(fullfile(savedir,strcat('Spot',int2str(setnum),'_',int2str(k),'.mat')),'RedI','GrI','FRET','fps')
                    clear RedToSave GrToSave FRETtoSave
                    cc=13;
                % Zoom
                elseif cc=='z'
                    [x,~] = ginput(2);
                    x = sort(x);
                    xlims(k,1) = x(1);
                    xlims(k,2) = x(2);
                    cc=13;
                % Unzoom
                elseif cc=='u'
                    xlims(k,:) = [0,0];
                    cc=13;
                elseif cc=='o'
                    offset = input('New offset:');
                    cc=13;
                % Set end of signal (FRET set to zero after this point)
                elseif cc=='e'
                    [x,~] = ginput(1);
                    ends(k) = round(x*fps/10^-3);
                    cc=13;
                % Don't let extra "enters" build up:
                elseif isequal(cc,char(13)) %13 is the ascii code for the return key
                    cc=13;
                %Show a specific frame in figure 2
                elseif cc=='f'
                    [x,~] = ginput(1);
                    [imgRinit,imgGinit] = PlayMovie(PathToMovie,[round(x) round(x)],h2,...
                        strcat('subplot(',char(39),'Position',char(39),...
                            ',[0.08 0.23 0.39 0.39*',int2str(size(imgRinit,1)/size(imgRinit,2)),'])'),...
                        strcat('subplot(',char(39),'Position',char(39),...
                            ',[0.54 0.23 0.39 0.39*',int2str(size(imgRinit,1)/size(imgRinit,2)),'])'),...
                        spots(:,k),GrSpots(:,k),...
                        strcat('subplot(',char(39),'Position',char(39),',[0.13 0.05 0.2 .18])'),...
                        strcat('subplot(',char(39),'Position',char(39),',[0.63 0.05 0.2 .18])'),...
                        zoomsize);
                    clear x
                    cc = 13;
                %Play a section of the movie in figure 2
                elseif cc=='m'
                    [x,~] = ginput(2);
                    x = round(sort(x));
                    [imgRinit,imgGinit] = PlayMovie(PathToMovie,[x(1) x(2)],h2,...
                        strcat('subplot(',char(39),'Position',char(39),...
                            ',[0.08 0.23 0.39 0.39*',int2str(size(imgRinit,1)/size(imgRinit,2)),'])'),...
                        strcat('subplot(',char(39),'Position',char(39),...
                            ',[0.54 0.23 0.39 0.39*',int2str(size(imgRinit,1)/size(imgRinit,2)),'])'),...
                        spots(:,k),GrSpots(:,k),...
                        strcat('subplot(',char(39),'Position',char(39),',[0.13 0.05 0.2 .18])'),...
                        strcat('subplot(',char(39),'Position',char(39),',[0.63 0.05 0.2 .18])'),...
                        zoomsize);
                    clear x
                    cc = 13;
                end
            end
       % end interactive section
        end
    % End loop over all spots 
    end
close all

end