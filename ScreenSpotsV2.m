%function newspots = ScreenSpots(img,spotcenters,boxdim,varargin)
%
%Interactive function that allows user to choose spots, from those found by 
%FindSpotsV3, that they don't want included in further analyses.
%
%Inputs are the image that spots were found in, a vector of (x,y)
%coordinates where the spot centers are, and boxdim that determines the
%size of the box to put around each spot.  The user should click in the box
%in order to deselect or re-select a spot.  Output is a vector of the kept
%spotcenters.
%
%Optional input is a string to put on the image, so you know what you're
%analyzing.
%
%This version uses an updated version of PutBoxesOnImage that allows figure
%data to be replaced without replotting the whole thing--makes it much faster!
%
%Steph 6/2013

function newspots = ScreenSpotsV2(img,spotcenters,boxdim,varargin)

if ~isempty(varargin) && ~ischar(varargin{1})
    disp('Spot screening: figure title must be a string.')
    newspots = -1;
    return
end

newspots = spotcenters; %Start out assuming user wants to keep all spots
badspots = [];

happy = 0;
coords=[];

h = PutBoxesOnImage(img,newspots,boxdim);
if ~isempty(varargin)
    title(strcat(varargin{1},'; Click on a spot to deselect; z to zoom; u to unzoom; return to end'));
else
    title('Click on a spot to deselect; z to zoom; u to unzoom; return to end');
end

while happy == 0     
    if ~isempty(coords)
        set(gca,'Xlim',[coords(1,1) coords(2,1)])
        set(gca,'Ylim',[coords(1,2) coords(2,2)])
    end
    [SelectedSpotx,SelectedSpoty,z] = ginput(1);
    if z==122 %User wants to zoom
        title('Click on upper left and lower right areas to zoom to:')
        coords = ginput(2);
    elseif z==117 %User wants to un-zoom
        set(gca,'Xlim',[1 size(img,2)])
        set(gca,'Ylim',[1 size(img,1)])
        coords=[];
    elseif isempty(SelectedSpotx)
        happy = 1;
    else
        %Find the nearest spot center to the user's click.  Start by
        %finding the distance from user's click to every spot:
        dists = sqrt((spotcenters(:,1)-SelectedSpotx.*ones(length(spotcenters),1)).^2+...
            (spotcenters(:,2)-SelectedSpoty.*ones(length(spotcenters),1)).^2);
        [~,closest] = min(dists);
        %The closest spot to the click is the one at
        %spotcenters(closest,:).
        %If this spot is in newspots, deselect and move to badspots:
        if ~isempty(find(ismember(newspots,spotcenters(closest,:))))
            badspots(end+1,:) = spotcenters(closest,:);
            badspots = sortrows(badspots);
            spotind = find(ismember(newspots,spotcenters(closest,:)));
            tempspots = newspots;
            clear newspots
            newspots = [tempspots(1:spotind(1)-1,:); tempspots(spotind(1)+1:end,:)];
            clear tempspots spotind
            %Don't need to sort, the newspots vector should still be sorted
        %if it's in badspots, re-select and move to newspots:
        elseif ~isempty(find(ismember(badspots,spotcenters(closest,:))))
            newspots(end+1,:) = spotcenters(closest,:);
            newspots = sortrows(newspots);
            spotind = find(ismember(badspots,spotcenters(closest,:)));
            tempspots = badspots;
            clear badspots
            badspots = [tempspots(1:spotind(1)-1,:); tempspots(spotind(1)+1:end,:)];
            clear tempspots spotind
        else
            disp('Spot is not in newspots or bad spots.  Problem!')
            keyboard
        end
       PutBoxesOnImage(img,newspots,boxdim,h);   
       
    end
end
close