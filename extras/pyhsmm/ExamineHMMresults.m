% function ExamineHMMresults(construct)
%
% Calls auto_analysis.py by "shelling out", to perform the HMM analysis on
% data for construct, then starts a simple GUI for the user to select
% traces for further analysis.
%
% Construct can be any the subdirectory names in resultsdir in config.py.
%
% Steph 2/2015

function ExamineHMMresults(construct)

%%% Parameters: %%%
% Figure position for GUI. Run figure('Position',fig_pos') to see where
% this will show up on your screen and adjust accordingly.
fig_pos = [100,400,1100,700];

% Minimum FRET cutoff. For visualization only.
minFRET = 0.775; % Plots a horizontal line indicating min acceptable starting FRET

% Where your python installation is:
setenv('PYTHONPATH', ['/Users/Steph/code/:', getenv('PYTHONPATH')]);
setenv('PATH', ['/Users/Steph/miniconda/bin/:', getenv('PATH')]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code guts below here %%%
[~,dirnames] = system(['python config.py']);
% The two outputs come back separated by a newline character:
newlines = regexpi(dirnames,'\n');
datadir = dirnames(1:newlines(1)-1);
resultsdir = dirnames(newlines(1)+1:end-1);
resultsdir = fullfile(resultsdir,construct);

if ~exist(fullfile(resultsdir,'ResultsFigs'),'dir')
    mkdir(fullfile(resultsdir,'ResultsFigs'));
end

python_command = ['python auto_analysis.py ', construct];

% Run the HMM analysis:
if exist(fullfile(resultsdir,'modelfits.pkl'),'file')
    doanalysis = 'y';
    doanalysis = input('Re-do HMM analysis? (y/n) ','s');
    if strcmpi(doanalysis,'y')
        system(python_command);
    end
else
    system(python_command);
end

% Load results
allresults = dir(fullfile(resultsdir,'*_Results.mat'));

if exist(fullfile(resultsdir,'ToAnalyzeFurther.mat'),'file')
    tokeep = load(fullfile(resultsdir,'ToAnalyzeFurther.mat'));
    tokeep = tokeep.tokeep;
    if length(tokeep)~=length(allresults)
        disp('Current number of results files does not match previous--clearing old analysis. Continue?')
        pause
        tokeep = ones(1,length(allresults));
    end
else
    tokeep = ones(1,length(allresults));
end

k=1;
figure('Position',fig_pos)
xlims = [0,0];
disp(',: back one trace; .: forward; d: _d_iscard (do not keep for further analysis)')
disp('or un-discard; z: zoom; u: unzoom; r: _r_edo Gibbs sampling on this one trace')

while k <= length(allresults)
    tempstruct = load(fullfile(resultsdir,allresults(k).name));
    
    xvectData = ((1:length(tempstruct.RedI))./tempstruct.fps);
    xvectModel = (double(1+tempstruct.start:tempstruct.start+length(tempstruct.model_redgreenseq(:,1)))./tempstruct.fps);
   
    if xlims(2)~=0
       xlimvect = ([xlims(1) xlims(2)]);
       datastyle = 'x';
    else
       xlimvect = ([0 min(tempstruct.start+length(tempstruct.model_redgreenseq(:,2))+200,length(tempstruct.RedI))./tempstruct.fps]);
       datastyle = '-';
    end
    
    subplot(2,2,1)
    if tokeep(k)
        plot(xvectData,tempstruct.RedI,strcat(datastyle,'r'),xvectData,tempstruct.GrI,strcat(datastyle,'g'))
    else
        plot(xvectData,tempstruct.RedI,strcat(datastyle,'k'),xvectData,tempstruct.GrI,strcat(datastyle,'k'))
    end
    hold on
    if tokeep(k)
        plot(xvectModel,tempstruct.model_redgreenseq(:,2),'--m',...
            xvectModel,tempstruct.model_redgreenseq(:,1),'--c','Linewidth',2)
    else
        plot(xvectModel,tempstruct.model_redgreenseq(:,2),'--k',...
            xvectModel,tempstruct.model_redgreenseq(:,1),'--k','Linewidth',2)
    end
    hold off
    set(gca,'Fontsize',12)
    xlabel('Time (s)','Fontsize',14)
    ylabel('Smoothed intensity (a.u.)','Fontsize',14)
    title(allresults(k).name,'Interpreter','none')
    xlim(xlimvect)
    
    % Update 4/2015: Adding two additional panels to plot the unsmoothed
    % data also
    subplot(2,2,2)
    if tokeep(k)
        plot(xvectData,tempstruct.unsmoothedRedI,strcat(datastyle,'r'),xvectData,tempstruct.unsmoothedGrI,strcat(datastyle,'g'))
    else
        plot(xvectData,tempstruct.unsmoothedRedI,strcat(datastyle,'k'),xvectData,tempstruct.unsmoothedGrI,strcat(datastyle,'k'))
    end
    hold on
    if tokeep(k)
        plot(xvectModel,tempstruct.model_redgreenseq(:,2),'--m',...
            xvectModel,tempstruct.model_redgreenseq(:,1),'--c','Linewidth',2)
    else
        plot(xvectModel,tempstruct.model_redgreenseq(:,2),'--k',...
            xvectModel,tempstruct.model_redgreenseq(:,1),'--k','Linewidth',2)
    end
    hold off
    set(gca,'Fontsize',12)
    xlabel('Time (s)','Fontsize',14)
    ylabel('Unsmoothed intensity (a.u.)','Fontsize',14)
    title(allresults(k).name,'Interpreter','none')
    xlim(xlimvect)
    
    subplot(2,2,3)
    if tokeep(k)
        plot(xvectData,tempstruct.FRET,strcat(datastyle,'b'))
    else
        plot(xvectData,tempstruct.FRET,strcat(datastyle,'k'))
    end
    hold on
    plot(xvectModel,tempstruct.model_fretseq,'--k','Linewidth',2)
    plot([xvectData(1) xvectData(end)],[minFRET minFRET],'--y')
    hold off
    xlabel('Time (s)','Fontsize',14)
    ylabel('Smoothed FRET','Fontsize',14)
    set(gca,'Fontsize',12)
    ylim([-0.2, 1.2])
    xlim(xlimvect)
    
    subplot(2,2,4)
    if tokeep(k)
        plot(xvectData,tempstruct.unsmoothedFRET,strcat(datastyle,'b'))
    else
        plot(xvectData,tempstruct.unsmoothedFRET,strcat(datastyle,'k'))
    end
    hold on
    plot(xvectModel,tempstruct.model_fretseq,'--k','Linewidth',2)
    plot([xvectData(1) xvectData(end)],[minFRET minFRET],'--y')
    hold off
    xlabel('Time (s)','Fontsize',14)
    ylabel('Unsmoothed FRET','Fontsize',14)
    set(gca,'Fontsize',12)
    ylim([-0.2, 1.2])
    xlim(xlimvect)
    
    cc=1;
    while cc~=13
        ct=waitforbuttonpress;
        cc=get(gcf,'currentcharacter');

        if ct==1
            % Go forward to the next trace
            if cc=='.'
                k = k+1;
                xlims = [0,0];
                cc = 13;
            % Go back one trace
            elseif cc==',' 
                if k>1
                    k=k-1;
                    xlims = [0,0];
                end
                cc=13;
            % Discard or un-discard this trace
            elseif cc=='d'
                tokeep(k) = ~tokeep(k);
                newtempstruct = load(fullfile(resultsdir,allresults(k).name));
                if isfield(newtempstruct,'discard')
                    discard = ~newtempstruct.discard;
                else
                    discard = 1;
                end
                save(fullfile(resultsdir,allresults(k).name),'discard','-append');
                clear newtempstruct discard
                cc=13;
            % Zoom
            elseif cc=='z'
                [x,~] = ginput(2);
                x = sort(x);
                if x(1)<0
                    x(1)=1;
                end
                if x(2)>max(tempstruct.start+length(tempstruct.model_redgreenseq(:,2)),length(tempstruct.RedI))./tempstruct.fps
                    x(2)=max(tempstruct.start+length(tempstruct.model_redgreenseq(:,2)),length(tempstruct.RedI))./tempstruct.fps;
                end
                xlims(1) = x(1);
                xlims(2) = x(2);
                cc=13;
            % Unzoom
            elseif cc=='u'
                xlims = [0,0];
                cc=13;
            % Redo Gibbs sampling on this one trace
            elseif cc=='r'
                name = allresults(k).name;
                name = name(1:end-12);
                python_command = ['python auto_analysis.py ', construct, ' ', name];
                system(python_command);
                clear name
                cc=13;
            % Don't let extra "enters" build up:
            elseif isequal(cc,char(13)) %13 is the ascii code for the return key
                cc=13;
            end
        end
        
    end
    clear tempstruct xvectData xvectModel
end
close

% The tokeep vector makes subsequent rounds of data easier. This
% information is also stored in each file.
save(fullfile(resultsdir,'ToAnalyzeFurther.mat'),'tokeep')