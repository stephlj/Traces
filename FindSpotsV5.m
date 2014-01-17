%function [spots,n,xout] = FindSpotsV5(img,varargin)
%
%Given an image (ie. a matrix of intensity values), find the centers of all 
%the spots. Also accepts a movie, in which case it will analyze the first ten frames.
%
%Does some spot screening: spots can't be too close to another spot or to
%the edge, where "too close" is defined by maxsize.
%
%V5, compared to V3, lets spots be closer together (or further apart) than
%the scale defined by neighborhoods for ordfilt2.
%
%Inputs:
%img: image (matrix of pixel intensities: x-by-y pixels, or x by y by time, 
%   where time is in the third dimension)
%
%Optional inputs: These must be entered in pairs:
%<parameter_name1>,<value1>,<parameter_name2>,<value2>, ...
%
%'maxsize',<val>: Spots will be excluded if they are closer than maxsize.
%'NeighborhoodSize',<val>: A little bigger than the size we expect spots to be 
%   (area); for ordfilt2. Neighborhood in which to look for local maxima and
%   minima. Needs to be a perfect square.  Best if sqrt(maxsize) is odd.  Default is 9^2.
%'ShowResults',1: to show an image at the end with boxes 
%   around all the spots (you can also use ScreenSpots for this, which has
%   an interactive component).
%'UseCentroid',1: take the centroid instead of the max pixel
%   value to find the spot center.  Currently centroid
%   actually does a worse job, at least for bead maps.
%'UserThresh',<val>: Instead of calculating a threshold for "true" maxima
%   automatically, will use val.
%'ImgTitle',<string>: If ShowResults is 1 and you want to give the image a title
%   (e.g., 'Green' so you know it's the green channel currently being analyzed).
%   Default is 'Spots'.
%
%Outputs:
%Spots:  a 2-by-numspots vector where each column is the (row,col) for the center of
%   each spot found in terms of the matrix of intensities that makes the
%   image.  Note row and col don't need to be integers here--but I'm keeping
%   the output in terms of the coordinate system of the matrix, not normal
%   Cartesian (x,y) coordinates.
%n,xout: outputs of the histogram of "diffs", in case the user wants to
%   choose a better threshold for true maxima
%
%Algorithm for finding spots:
%(1) Identify pixels that are their maxes or mins in a local
%neighborhood.
%(2) Identify as true local maxes those maxes that are significantly bigger
%than the surrounding minima.
%(3) Do some light screening of spots: don't include any too close to the
%   edges of the image.
%(4) Return either the position of the max intensity or the centroid around
%   each true maximum as a spot.
%
%Steph 8/2013, updated 10/2013
%Copyright 2013 Stephanie Johnson, University of California, San Francisco

function [spots,n,xout] = FindSpotsV5(img,varargin)

%%%%TODO: 
%%%%Likelihood estimation for better threshold detection
%%%%Figure out centroid thing

%Defaults for optional inputs:
NeighborhoodSize = 9^2;
maxsize = 6;
ShowResults = 0;
UseCentroid = 0;
UserThresh = 0;
boxsize = maxsize;
ImgTitle = 'Spots';

%Other defaults:
cen_boxsize = 2;%Size of the box to use for centroid calculation will be cen_boxsize*2+1
mov_fr = 10; %Number of frames to average over, if given a movie 

%Optional inputs:
if ~isempty(varargin)
    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'ShowResults')
            ShowResults = varargin{i+1};
        elseif strcmpi(varargin{i},'UseCentroid')
            UseCentroid = varargin{i+1};
        elseif strcmpi(varargin{i},'UserThresh')
            UserThresh = 1;
            truemaxthresh = varargin{i+1};
        elseif strcmpi(varargin{i},'ImgTitle')
            ImgTitle = varargin{i+1};
        elseif strcmpi(varargin{i},'maxsize')
            maxsize = varargin{i+1};
            boxsize = maxsize;
        elseif strcmpi(varargin{i},'NeighborhoodSize')
            NeighborhoodSize = varargin{i+1};
        end
    end
end

%Error handling for inputs:
if rem(NeighborhoodSize,sqrt(NeighborhoodSize))~=0
    disp('Spot finding: NeighborhoodSize must be a perfect square.')
    spots = -1;
    n = -1;
    xout = -1;
    return
end

%Make sure the image is properly scaled, and average if given a movie:
if size(img,3)==1
    imggray = mat2gray(img);
    mov = 0;
elseif size(img,3)>=mov_fr
    imggray = mat2gray(sum(img(:,:,1:mov_fr),3)./mov_fr);
    mov = 1;
else
    imggray = mat2gray(sum(img,3)./size(img,3));
    mov = 1;
end

%May need to do some background subtraction and/or smoothing here, or maybe
%only for the non-movie inputs

%(1) Define a neighborhood as NeighborhoodSize square pixels--want the max and min
%calculations to be over areas roughly the size we expect beads to be
maxes = ordfilt2(imggray,NeighborhoodSize,ones(sqrt(NeighborhoodSize),sqrt(NeighborhoodSize))); 
    %Note about ordfilt: This tests over square neighborhoods; make a matrix that's not ones
    %everywhere to look over some other kind of domain shape.
mins  = ordfilt2(imggray,1,ones(sqrt(NeighborhoodSize),sqrt(NeighborhoodSize)));

%(2) Identify true maxes.  
%Need to define a threshold for true maxes.
%If we histogram diffs, we'll see two peaks, one for true maxes, one for
%background.  Best way to pick a threshold is to use a likelihood ratio test
%to pick where the diffs are most likely to be part of the true maxes peak
%and unlikely to be part of the background peak.  For now, just picking the
%mean value between the centerpoints of two gaussians fit to a histogram of
%the diffs (doesn't work great but works ok mostly)
diffs = maxes-mins;
[n,xout] = hist(reshape(diffs,1,numel(diffs)),50);
if ~UserThresh
    %Update 12/2013: Since not all the computers we run this code on have
    %the Curve fitting toolbox, default to threshold of 0.15 if necessary.
    %Turns out this does pretty well for most images.
    try
        opts = fitoptions('gauss2','Algorithm','Trust-Region','lower',[0 0 0 0 0.15 0],...
            'upper',[Inf 1 Inf Inf 1 Inf],'StartPoint',[1 0.1 1 1 0.5 1]);
        fitparams = fit(xout',n','gauss2',opts);
        truemaxthresh = mean([fitparams.b1,fitparams.b2]);
    catch
        truemaxthresh = 0.15;
    end
end
truemaxes = (maxes-mins) > truemaxthresh;
%Finally find where the value of the image equals the local max, and keep
%this local max only if truemax for that pixel is 1:
temppeaks = find((imggray==maxes).*truemaxes);
%Find uses linear indexing (which is done column-wise), so convert these to (row,col):
[peaks_i,peaks_j] = ind2sub(size(imggray),temppeaks);
peaks = [peaks_i,peaks_j];

disp(strcat('Number of peaks found: ',int2str(size(peaks,1))))

%For deciding if spots are too close together: first need to find all the
%distances between spots:
Dists = FindSpotDists(peaks);
%Declare that spots can't have centers closer than maxsize:
spottooclose = Dists>maxsize;
%Each row will be all 1's if the spot represented by the row is more than
%maxsize away from another spot.  If there's another spot too
%close, one or more elements will be zero.  So below, ask if the sum of a
%peak's row is equal to the length of the row (minus one, because of the
%zero element where it's too close to itself)

%Exclude any peaks that are too close to the boundaries or to another spot, 
%and calculate centroid if user desires:
spots = [];
for i = 1:size(peaks,1)
    if (peaks(i,1)-maxsize)>=1 && (peaks(i,1)+maxsize)<=size(imggray,1) && ...
            (peaks(i,2)-maxsize)>=1 && (peaks(i,2)+maxsize)<=size(imggray,2) && ...
            sum(spottooclose(i,:))==length(spottooclose(i,:))-1
        %Keep this peak.  
        if UseCentroid
            cen = FindSpotCentersV2(imggray((peaks(i,1)-cen_boxsize):(peaks(i,1)+cen_boxsize),...
                (peaks(i,2)-cen_boxsize):(peaks(i,2)+cen_boxsize)));
            %Because the center-finding algorithms were passed only a
            %submatrix (ie a cropped image of each spot), they returned the (row,col) of
            %the center of the spot in terms of submatrix.  Want to
            %return centers relative to the original imggray matrix:
            spots(end+1,:) = [peaks(i,1)-cen_boxsize-1+cen(1), peaks(i,2)-cen_boxsize-1+cen(2)];
        else
            spots(end+1,:) = [peaks(i,1),peaks(i,2)];
        end
    end
end

disp(strcat('Number of peaks kept: ',int2str(size(spots,1))))

%For checking the results: make a figure that puts a box around each spot
%found.  Also put a red x where it thinks the center of each spot is.
%Updated 1/2014 for some error handling in case no spots were found:
if ShowResults && size(spots,1) > 0 && size(spots,2) > 0
    PutBoxesOnImageV4(imggray,spots,boxsize,0);
    title(ImgTitle,'Fontsize',16)
end


