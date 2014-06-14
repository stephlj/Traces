% function [matchedR,matchedG,RedI,GrI,unmatchedR,unmatchedG,unmatchedRedI,...
%    unmatchedGrI] = CorrelateIntensities(RedSpots,VarsR,GrSpots,VarsG,PathToMovie,params,roughtform)
%
% Given a set of red spots and a set of green spots, calculates their
% intensities through a movie, and matches those that are anti-correlated
% (i.e. FRETing).
%
% Steph 6/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [matchedR,matchedG,RedI,GrI,unmatchedR,unmatchedG,unmatchedRedI,...
    unmatchedGrI] = CorrelateIntensities(RedSpots,VarsR,GrSpots,VarsG,PathToMovie,params,roughtform)

max_dist = 10; % How many pixels away to consider spots in the other channel as possible matches?

if exist('roughtform','var')
    GrSpots = roughtform.FRETmapFwd(GrSpots);
end

% Function for use later:
    function I = CalcOneSpotIntensity(spotcen,spotvars,channel,params,numtifs,path)
        I = zeros(1,numtifs);
        if rem(100,params.FramesToAvg)==0
            FramesToAvg = params.FramesToAvg;
        else
            FramesToAvg = 10;
        end
        for ff = 1:100:numtifs
            if strcmpi(channel,'green')
                [~,imgs] = LoadScaledMovie(path,[ff ff+99]);
            else
                [imgs,~] = LoadScaledMovie(path,[ff ff+99]);
            end
            [spotimg,localcen] = ExtractROI(imgs,params.DNASize,spotcen);
            
            if params.IntensityGaussWeight==1
                bkgnd = zeros(1,size(spotimg,3));
                prevbkgnd = min(min(spotimg(:,:,1)));
                prevA = max(max(spotimg(:,:,1)));
                for bb = 1:FramesToAvg:size(spotimg,3)
                    [~,~,~,~,prevbkgnd,prevA] = Fit2DGaussToSpot(mean(spotimg(:,:,bb:bb+FramesToAvg-1),3),...
                        'Background','StartParams',[localcen(1),localcen(2),spotvars(1),spotvars(2),...
                        prevbkgnd,prevA],'symGauss',params.UseSymGauss);
                    bkgnd(bb:bb+FramesToAvg-1) = prevbkgnd;
                end
                % Subtract this background value from every pixel in the ROI
                % containing the spot:
                bkgnd = repmat(bkgnd,size(spotimg,1)*size(spotimg,2),1,1);
                bkgnd = reshape(bkgnd,size(spotimg,1),size(spotimg,2),size(spotimg,3));
                I(ff:ff+size(imgs,3)-1) = CalcSpotIntensityV4(spotimg-bkgnd,...
                    localcen,spotvars);
            else
                I(ff:ff+size(imgs,3)-1) = CalcSpotIntensityNoGauss(imgs,spotcen);
            end
        end
        
    end

alltifs = dir(fullfile(PathToMovie,'img*.tif'));

matchedR = [];
matchedG = [];
RedI = [];
GrI = [];
unmatchedR = zeros(size(RedSpots));
unmatchedG = zeros(size(GrSpots));
unmatchedRedI = zeros(size(RedSpots,2),length(alltifs));
unmatchedGrI = zeros(size(GrSpots,2),length(alltifs));

% Iterate through the smallest number of spots possible, since the
% intensity calculation takes a while. Although doing it this way does mean
% that I have to load the movie over and over and over again ... 
if size(RedSpots,2)<=size(GrSpots,2)
    base_spots = RedSpots;
    spots_to_match = GrSpots;
    base_vars = VarsR;
    match_vars = VarsG;
    base_channel = 'Red';
    match_channel = 'Green';
else
    base_spots = GrSpots;
    spots_to_match = RedSpots;
    base_vars = VarsG;
    match_vars = VarsR;
    base_channel = 'Green';
    match_channel = 'Red';
end

indices = 1:size(spots_to_match,2);
Dists = FindSpotDists(base_spots,spots_to_match);

% Iterate through base_spots, looking for nearby spots_to_match
for ss = 1:size(base_spots,2)
    disp(sprintf('Spot %d of %d:',ss,size(base_spots,2)))
    CloseMatches = indices(Dists(ss,:)<=max_dist); % all spots_to_match that are
        % within max_dist of base_spot
    disp(sprintf('%d nearby spots;',length(CloseMatches)))
    if ~isempty(CloseMatches)
        % Calculate this spot's intensity as a function of time
        tic
        base_I = CalcOneSpotIntensity(base_spots(:,ss),base_vars(:,ss),...
            base_channel,params,length(alltifs),PathToMovie);
        toc
        
        match_I = zeros(length(CloseMatches),length(alltifs));
        spotcorrs = zeros(1,length(CloseMatches));
        
        for gg = 1:length(CloseMatches)
            match_I(gg,:) = CalcOneSpotIntensity(spots_to_match(:,CloseMatches(gg)),match_vars(:,CloseMatches(gg)),...
                match_channel,params,length(alltifs),PathToMovie);
            spotcorrs(gg) = sum(sum(corrcoef(base_I',match_I(gg,:)')))
            % Anticorrelated: should result in sum(sum([1 -1; -1 1]))=0
            % Correlated: should result in sum(sum([1 1; 1 1]))=4
            % Rand gives ~ sum(sum([1 -0.5 -0.5 1])) ~ 2
            
            % for debugging
                figure
                plot(base_I,'-r')
                hold on
                plot(match_I(gg,:),'-g')
                hold off
                %base_spots(:,ss)
                %spots_to_match(:,CloseMatches(gg))
                pause
                close
        end
        [closestcorr,ind] = min(spotcorrs);
        if closestcorr<1.5
            % These guys are matched
            disp('Found a match!')
            if strcmpi(base_channel,'Red')
                matchedR(:,end+1) = base_spots(:,ss);
                RedI(end+1,:) = base_I;
                matchedG(:,end+1) = spots_to_match(:,CloseMatches(ind));
                GrI(end+1,:) = match_I(CloseMatches(ind),:);
            else
                matchedG(:,end+1) = base_spots(:,ss);
                GrI(end+1,:) = base_I;
                matchedR(:,end+1) = spots_to_match(:,CloseMatches(ind));
                RedI(end+1,:) = match_I(CloseMatches(ind),:);
            end
        else
            disp('No matches found.')
            % Keep the intensity calculation we already did anyway
            if strcmpi(base_channel,'Red')
                unmatchedR(:,ss) = base_spots(:,ss);
                unmatchedRedI(ss,:) = base_I;
                unmatchedG(:,CloseMatches) = spots_to_match(:,CloseMatches);
                unmatchedGrI(CloseMatches,:) = match_I;
            else
                unmatchedG(:,ss) = base_spots(:,ss);
                unmatchedGrI(ss,:) = base_I;
                unmatchedR(:,CloseMatches) = spots_to_match(:,CloseMatches);
                unmatchedRedI(CloseMatches,:) = match_I;
            end
                
        end
    end
    
end

if exist('roughtform','var')
    matchedG = roughtform.FRETmapInv(matchedG);
    unmatchedG = roughtform.FRETmapInv(unmatchedG);
end

end