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
% If FigID==3
%   varargin{1} = RefinedCentersR or G;
%   varargin{2} = newspotsR or G;
%   varargin{3} = 'Red' or 'Green'
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
        % Update 10/2014 making this more universal for macs and pcs
        %figure('Position',[200 200 325 625])
        % Update 1/2015: This still doesn't work very well.  New plan:
%         figure
%         DefaultFigPos = get(gcf,'Position');
%         set(gcf,'Position',[DefaultFigPos(1)-0.5*DefaultFigPos(1), DefaultFigPos(2),...
%             DefaultFigPos(3)*0.575, DefaultFigPos(4)*1.5])
%         clear DefaultFigPos
        if ismac
            figure('Position',[200 200 325 625])
        else
            figure('Position',[200 45 325 625])
        end
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
elseif FigID==3
    if length(varargin)==3
        % Update 10/2014 making these figures show up properly on both macs
        % and pcs
        % figure('Position',[200 200 325 625])
%         figure
%         DefaultFigPos = get(gcf,'Position');
%         set(gcf,'Position',[DefaultFigPos(1)-0.5*DefaultFigPos(1), DefaultFigPos(2),...
%             DefaultFigPos(3)*0.575, DefaultFigPos(4)*1.5])
%         clear DefaultFigPos
        if ismac
            figure('Position',[200 200 325 625])
        else
            figure('Position',[200 45 325 625])
        end
        hold on
        if ~isempty(varargin{1})
            plot(varargin{1}(2,:),0-varargin{1}(1,:),'xg')
        end
        if ~isempty(varargin{2})
            plot(varargin{2}(2,:),0-varargin{2}(1,:),'xr')
        end
        hold off
        ylim([-512 0])
        xlim([0 256])
        title(strcat(varargin{3},' Channel: Newly added spots in red'))
        pause
        close
    else
        disp('PlotDebugFigures: Not the right number of inputs with FigID = 3?')
        return
    end
else
    disp('PlotDebugFigures: Unknown FigID.')
end