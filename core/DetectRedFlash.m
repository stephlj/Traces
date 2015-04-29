% function params = DetectRedFlash(AcceptorMaxes,Medians,params)
%
% Called by ScaleMovie if the DetectRedFlash parameter in TracesSetup is
% nonzero. Looks for flashes of high intensity in the acceptor channel
% which we use to mark injection points. 
%
% A movie can contain multiple red flashes. A red flash needs to meet the
% following criteria in order to be detected:
% (1) A significant spike in acceptor channel intensity
% (2) A concurrent spike in the median intensity over the whole frame
% (excludes random spikes from very brief, bright junk spots in the
% acceptor channel)
% (3) No concurrent spike in the donor channel
% (4) Spikes must last more than one frame
%
% Steph 4/2015
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

function params = DetectRedFlash(AcceptorMaxes,DonorMaxes,Medians,params)

    % Subfunction that finds potential intensity spikes.
    function candidates = FindSpikes(searchvec)
        candidates = [];
        
        % Find values in the vector that are much bigger than the values 
        % around it. Start by defining a threshold above which a value is 
        % "much bigger" than those around it:
        spikethresh = 4; %We will insist that a spike be 4 std dev's bigger
            % than the mean value of searchvec.

        maxes = searchvec > (mean(searchvec)+spikethresh*std(searchvec));

        % Assume that the first and last couple frames are never a fiducial
        % mark:
        maxes(1:3) = zeros(3,1);
        maxes(end-2:end) = zeros(3,1);
        
        maxes = find(maxes);
        
        % Lastly insist that a true candidate span more than one frame, and
        % identify the flash frame as the first in the set of maxes:
        
        k = length(maxes);
        while k > 1
            if maxes(k-1) == maxes(k)-1
                % This one spans more than one frame, keep it.
                % find where it starts:
                while k>1 && maxes(k-1) == maxes(k)-1
                    k = k-1;
                end
                candidates(end+1) = maxes(k);
                k = k-1;
            else
                k = k-1;
            end
        end
        
        candidates = sort(candidates);
    end

    % First find all potential flashes in the acceptor channel:
    RedCandidates = FindSpikes(AcceptorMaxes);
    % Check whether these are also spikes in the medians:
    MedCandidates = FindSpikes(Medians);
    % and no concurrent spikes in the donor channel:
    GrCandidates = FindSpikes(DonorMaxes);
    
    flashes = [];
    
    if isempty(RedCandidates) || isempty(MedCandidates)
        disp('No red flashes detected ... ')
        keyboard;
        
    else
        params.ManualInjectMark = 0;
        % Iterate through each candidate and ask make sure it appears in both lists:
        for j = 1:length(RedCandidates)
            if ismember(RedCandidates(j),MedCandidates) && ~ismember(RedCandidates(j),GrCandidates)
                flashes(end+1) = RedCandidates(j);
            end
        end

        % Lastly, ask the user to verify:
        figure
        isok = 'n';
        while ~isempty(isok)
            plot(AcceptorMaxes,'xr')
            hold on
            plot(DonorMaxes,'xg')
            plot(Medians,'ob')
            for f = 1:length(flashes)
                plot([flashes(f) flashes(f)],[0 max(AcceptorMaxes)],'--k')
            end
            legend('Acceptor channel','Donor channel','Medians')
            hold off
            ylabel('Intensity (a.u.)','Fontsize',14)
            xlabel('Frame','Fontsize',14)
            set(gca,'Fontsize',14)
            isok = input('Press enter if this is right; a to add a missing flash, r to remove one:','s');
            if strcmpi(isok,'a')
                disp('Click where you want to add a flash:')
                [x,~] = ginput(1);
                flashes(end+1) = round(x);
                clear x
            elseif strcmpi(isok,'r') && length(flashes)>1
                disp('Click on a flash to remove:')
                [x,~] = ginput(1);
                disttoflashes = abs(flashes-x);
                if min(disttoflashes) == disttoflashes(1)
                    flashes = flashes(2:end);
                elseif min(disttoflashes) == disttoflashes(end)
                    flashes = flashes(1:end-1);
                else
                    [~,minidx] = min(disttoflashes);
                    flashes = [flashes(1:minidx-1),flashes(minidx+1:end)];
                end
                clear x disttoflashes
            elseif strcmpi(isok,'r') 
                disp('You must end up with at least one flash; add one if you want to remove this one.')
            end
        end
    end
    params.InjectPoints = sort(flashes)/(params.fps*10^3);
    
end



