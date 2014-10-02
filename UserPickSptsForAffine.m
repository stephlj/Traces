% function [newspotsR,newspotsG] = UserPickSptsForAffine(imgR,imgG,spotsR,spotsG,boxdim)
%
% Interactive function that allows user to choose spots, from those found by 
% FindSpotsV5, that they think will form a good set of pairs from which to
% calculate a rough affine map. This will then be used to pair the rest of
% the spots so that a polynomial map can be calculated.  For use when the
% channel alignment is so off that the greedy algorithm in FindSpotMatches
% fails.
%
% Inputs are the images that spots were found in (imgR,imgG, as generated
% by SplitImg), matrices of (x;y) by numspots containing coordinates where 
% the spot centers are (spotsR, spotsG), and boxdim that determines the
% size of the circle to put around each spot.  The user should click in a
% circle in order to deselect or select a spot.  Output is a vector of 
% the newly matched spotcenters.
%
% Copyright (C) 2014 Stephanie Johnson, University of California, San Francisco
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

function [newspotsR,newspotsG] = UserPickSptsForAffine(imgR,imgG,spotsR,spotsG,boxdim)

    % Make sure the spots are listed as one spot per column:
    if size(spotsR,1)~=2
        spotsR = transpose(spotsR);
    end
    if size(spotsG,1)~=2
        spotsG = transpose(spotsG);
    end

        % Subfunction to find spot closest to where user clicked:
        function [unused,kept] = FindUserSelection(allspots,unused,kept,spotclicked)
            % Find the nearest spot center to the user's click.  Start by
            % finding the distance from user's click to every spot:
            dists = FindSpotDists(spotclicked,allspots);
            [minDist,closest] = min(dists,[],2);
            % The closest spot to the click is the one at
            % spotcenters(closest,:).
            if minDist < boxdim*2
                % If this spot is in unused, move to kept, as long as there's not
                % already a spot in kept:
                if length(find(ismember(unused,allspots(:,closest))))==2 && isempty(kept)
                    kept(:,end+1) = allspots(:,closest);
                    %kept = transpose(sortrows(kept'));
                    spotind = find(ismember(unused(1,:),allspots(1,closest)));
                    tempspots = unused;
                    clear unused
                    unused = [tempspots(:,1:spotind-1), tempspots(:,spotind+1:end)];
                    clear tempspots spotind
                    % Don't need to sort, the unused vector should still be sorted
                % if it's in kept, re-select and move to unused:
                elseif length(find(ismember(kept,allspots(:,closest))))==2
                    unused(:,end+1) = allspots(:,closest);
                    %unused = transpose(sortrows(unused'));
                    spotind = find(ismember(kept(1,:),allspots(1,closest)));
                    tempspots = kept;
                    clear kept
                    kept = [tempspots(:,1:spotind-1), tempspots(:,spotind+1:end)];
                    clear tempspots spotind
                end
            end
        end

    done = 0;
    
    while done==0
        unusedspotsR = spotsR; % Start out assuming user wants to keep all spots
        spotsRforAffine = [];
        unusedspotsG = spotsG;
        spotsGforAffine = [];

        % Plot the two images side by side, with circles around all the spots:
        % Update 10/2014: To make this more compatible with both PC's and
        % Mac's,
        figure
        DefaultFigPos = get(gcf,'Position');
        set(gcf,'Position',[DefaultFigPos(1)-0.5*DefaultFigPos(1), DefaultFigPos(2),...
            DefaultFigPos(3)*1.25, DefaultFigPos(4)*1.5])
        clear DefaultFigPos
        %figure('Position',[200 200 700 625])
        t = 0:pi/100:2*pi;

        Raxes = subplot(1,2,1);
        imshow(imgR,[])
        hold on
        for j = 1:size(spotsR,2)
            plot(spotsR(2,j)+boxdim/2*cos(t),spotsR(1,j)+boxdim/2*sin(t),'-r')
        end
        hold off
        Gaxes = subplot(1,2,2);
        imshow(imgG,[])
        hold on
        for j = 1:size(spotsG,2)
            plot(spotsG(2,j)+boxdim/2*cos(t),spotsG(1,j)+boxdim/2*sin(t),'-g')
        end

        for r = 1:3 % Three is the minimum number of points to define an affine
            z = 5;
            rspottemp = [];
            gspottemp = [];

            disp('Click on a red spot to keep; click again to un-keep it. Press enter when satisfied.')
            while ~isempty(z)
                [SelectedSpotx,SelectedSpoty,z] = ginput(1);
                if isequal(gca,Raxes) && ~isempty(z)
                    [unusedspotsR,rspottemp] = FindUserSelection(spotsR,unusedspotsR,...
                    rspottemp,[SelectedSpoty;SelectedSpotx]);
                    subplot(1,2,1)
                    imshow(imgR,[])
                    hold on
                    for j = 1:size(unusedspotsR,2)
                        plot(unusedspotsR(2,j)+boxdim/2*cos(t),unusedspotsR(1,j)+boxdim/2*sin(t),'-r')
                    end
                    for j = 1:size(spotsRforAffine,2)
                        plot(spotsRforAffine(2,j)+boxdim/2*cos(t),spotsRforAffine(1,j)+boxdim/2*sin(t),'-w')
                    end
                    for j = 1:size(rspottemp,2)
                        plot(rspottemp(2,j)+boxdim/2*cos(t),rspottemp(1,j)+boxdim/2*sin(t),'-w')
                    end
                    hold off
                    drawnow
                end
            end
            spotsRforAffine(:,end+1) = rspottemp;

            z = 5;

            disp('Now click on its match in the green channel; click again to un-keep it. Press enter when satisfied.')
            while ~isempty(z)
                [SelectedSpotx,SelectedSpoty,z] = ginput(1);
                if isequal(gca,Gaxes) && ~isempty(z)
                    [unusedspotsG,gspottemp] = FindUserSelection(spotsG,unusedspotsG,...
                        gspottemp,[SelectedSpoty;SelectedSpotx]);
                    subplot(1,2,2)
                    imshow(imgG,[])
                    hold on
                    for j = 1:size(unusedspotsG,2)
                        plot(unusedspotsG(2,j)+boxdim/2*cos(t),unusedspotsG(1,j)+boxdim/2*sin(t),'-g')
                    end
                    for j = 1:size(spotsGforAffine,2)
                        plot(spotsGforAffine(2,j)+boxdim/2*cos(t),spotsGforAffine(1,j)+boxdim/2*sin(t),'-w')
                    end
                    for j = 1:size(gspottemp,2)
                        plot(gspottemp(2,j)+boxdim/2*cos(t),gspottemp(1,j)+boxdim/2*sin(t),'-w')
                    end
                    hold off
                    drawnow
                end
            end
            spotsGforAffine(:,end+1) = gspottemp;
        end

        close

        % Make an affine transformation with these user-paired spots
        tformAffine = FRETmap(spotsGforAffine,spotsRforAffine,'Green','MatlabAffine');

        % Re-pair spots and return new matches:
        [newspotsGtemp,newspotsR] = FindSpotMatches(tformAffine.FRETmapFwd(spotsG),spotsR);
        if newspotsGtemp == -1 % Still didn't work
            disp('Affine transformation not good enough, try picking different spot-pairs.')
        else
            done = 1;
        end
    end
    newspotsG = tformAffine.FRETmapInv(newspotsGtemp);
end