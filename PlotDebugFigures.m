% function PlotDebugFigures(FigID,varargin)
%
% In the interest of cleaning of the smFRET wrapper function, putting calls
% to the debug figures in this separate function.
%
% If FigID==1:
%   varargin{1} = allBdImgs(:,:,i);
%   varargin{2} = spotsR_abs;
%   varargin{3} = spotsG_abs;
%   varargin(4} = params;
% If FigID==2:
%   varargin{1} = allBdImgs(:,:,i);
%   varargin{2} = matchR_abs;
%   varargin{3} = matchG_abs;
%   varargin(4} = params;
%   varargin{5} = matchR{i};
%   varargin{6} = matchG{i};
%   varargin{7} = size(imgRed);

function PlotDebugFigures(FigID,varargin)

if FigID==1%Figure with all the spots found during the channel mapping routine
    if length(varargin)==4
        % plot all the boxes for both channels on a big image:
        PutBoxesOnImageV4(varargin{1},[varargin{2},varargin{3}],varargin{4}.BeadSize,'0','w');
        pause
        close
    else
        disp('PlotDebugFigures: Not the right number of inputs with FigID = 1?')
        return
    end
elseif FigID == 2 %Box in green the ones that were matched:
    if length(varargin)==7
        PutBoxesOnImageV4(varargin{1},[varargin{2}';varargin{3}'],varargin{4}.BeadSize);
        title('Only spots that were matched')
        % Another way of plotting the matching: blue line between points
        % in the two channels, with a green dot for where the point is
        % in the green channel, and the end of the blue line where the
        % point it got matched to in the red channel is.  This is also a
        % great way of looking at the distortion between the two
        % channels
        figure('Position',[200 200 325 625])
        plot([varargin{5}(2,:);varargin{6}(2,:)],0-[varargin{5}(1,:);varargin{6}(1,:)],'-b')
        hold on
        plot(varargin{5}(2,:),0-varargin{5}(1,:),'xr')
        ylim([-varargin{7}(1), 0])
        xlim([0 varargin{7}(2)])
        title('Red x is center of point in red channel')
        pause
        close all
    else
        disp('PlotDebugFigures: Not the right number of inputs with FigID = 2?')
        return
    end
end