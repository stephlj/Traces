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
    unmatchedGrI] = CorrelateIntensitiesManual(RedSpots,VarsR,GrSpots,VarsG,PathToMovie,params,roughtform)

if exist('roughtform','var')
    GrSpotsInR = roughtform.FRETmapFwd(GrSpots);
    max_dist = 3; % How many pixels away to consider spots in the other channel as possible matches?
else
    GrSpotsInR = GrSpots;
    max_dist = 10;
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
            prevbkgnd = min(min(spotimg(:,:,1)));
            prevA = max(max(spotimg(:,:,1)));
            
            if params.EndInjectFrame>1 && ff<params.EndInjectFrame
                p=1;
                while p <= params.EndInjectFrame && p<=(size(spotimg,3)-FramesToAvg)
                    prevparams = [localcen(1),localcen(2),spotvars(1),spotvars(2),...
                        prevbkgnd,prevA];
                    [localcen(1),localcen(2),spotvars(1),spotvars(2),prevbkgnd,prevA] = Fit2DGaussToSpot(mean(spotimg(:,:,p:p+FramesToAvg-1),3),...
                            'Full','StartParams',prevparams,'symGauss',params.UseSymGauss);
                    bkgnd = repmat(prevbkgnd,size(spotimg,1)*size(spotimg,2),1);
                    bkgnd = reshape(bkgnd,size(spotimg,1),size(spotimg,2));
                    if params.IntensityGaussWeight==1
                        I(ff+p-1) = CalcSpotIntensityV4(spotimg(:,:,p)-bkgnd,...
                            localcen,spotvars);
                    else
                        I(ff+p-1) = CalcSpotIntensityNoGauss(spotimg(:,:,p)-bkgnd,localcen);
                    end
                    p=p+1;
                end
                if p<size(spotimg,3)
                    bkgnd = zeros(1,size(spotimg,3)-p+1);
                    for bb = p:FramesToAvg:size(spotimg,3)
                        [~,~,~,~,prevbkgnd,prevA] = Fit2DGaussToSpot(mean(spotimg(:,:,bb:min(bb+FramesToAvg-1,size(spotimg,3))),3),...
                            'Background','StartParams',[localcen(1),localcen(2),spotvars(1),spotvars(2),...
                            prevbkgnd,prevA],'symGauss',params.UseSymGauss);
                        bkgnd(bb-p+1:bb-p+1+length(prevbkgnd)) = prevbkgnd;
                    end
                    % Subtract this background value from every pixel in the ROI
                    % containing the spot:
                    bkgnd = repmat(bkgnd,size(spotimg,1)*size(spotimg,2),1);
                    bkgnd = reshape(bkgnd,size(spotimg,1),size(spotimg,2),size(spotimg,3)-p+1);
                    if params.IntensityGaussWeight==1
                        I(ff+p-1:ff+p-2+size(bkgnd,3)) = CalcSpotIntensityV4(spotimg(:,:,p:end)-bkgnd,...
                            localcen,spotvars);
        %                 I(ff:ff+size(imgs,3)-1) = CalcSpotIntensityV4(spotimg-bkgnd,...
        %                     localcen,[0.3;0.3]);
                    else
                        I(ff+p-1:ff+p-1+size(bkgnd,3)) = CalcSpotIntensityNoGauss(spotimg(:,:,p:end)-bkgnd,localcen);
                    end
                end
                spotcen = GlobalToROICoords([],localcen,spotcen,params.DNASize,params.DNASize);
            else
                %if params.IntensityGaussWeight==1
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
                if params.IntensityGaussWeight==1
                    I(ff:ff+size(imgs,3)-1) = CalcSpotIntensityV4(spotimg-bkgnd,...
                        localcen,spotvars);
    %                 I(ff:ff+size(imgs,3)-1) = CalcSpotIntensityV4(spotimg-bkgnd,...
    %                     localcen,[0.3;0.3]);
                else
                    I(ff:ff+size(imgs,3)-1) = CalcSpotIntensityNoGauss(spotimg-bkgnd,localcen);
                end
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
%if size(RedSpots,2)<=size(GrSpots,2)
    base_spots = RedSpots;
    spots_to_match = GrSpotsInR;
    base_vars = VarsR;
    match_vars = VarsG;
    base_channel = 'Red';
    match_channel = 'Green';
% else
%     base_spots = GrSpotsInR;
%     spots_to_match = RedSpots;
%     base_vars = VarsG;
%     match_vars = VarsR;
%     base_channel = 'Green';
%     match_channel = 'Red';
% end

indices = 1:size(spots_to_match,2);
Dists = FindSpotDists(base_spots,spots_to_match);

% Iterate through base_spots, looking for nearby spots_to_match
for ss = 1:size(base_spots,2)
    close all
    figure('Position',[200 200 325 625])
    plot(base_spots(2,:),0-base_spots(1,:),'xr')
    hold on
    plot([base_spots(2,ss) base_spots(2,ss)],[-512 0],'--r')
    plot([0 256],0-[base_spots(1,ss) base_spots(1,ss)],'--r')
    if exist('roughtform','var')
        plot(GrSpotsInR(2,:),0-GrSpotsInR(1,:),'xg')
    else
        plot(GrSpots(2,:),0-GrSpots(1,:),'xg')
    end
    ylim([-512 0])
    xlim([0 256])
    disp(sprintf('Spot %d of %d at (%d,%d):',ss,size(base_spots,2),base_spots(1,ss),base_spots(2,ss)))
    CloseMatches = indices(Dists(ss,:)<=max_dist); % all spots_to_match that are
        % within max_dist of base_spot
    disp(sprintf('%d nearby spots;',length(CloseMatches)))
    if ~isempty(CloseMatches)
        % Calculate this spot's intensity as a function of time
        base_I = CalcOneSpotIntensity(base_spots(:,ss),base_vars(:,ss),...
            base_channel,params,length(alltifs),PathToMovie);
        
        match_I = zeros(length(CloseMatches),length(alltifs));
        spotcorrs = zeros(1,length(CloseMatches));
        
        for gg = 1:length(CloseMatches)
            plot([spots_to_match(2,CloseMatches(gg)) spots_to_match(2,CloseMatches(gg))],[-512 0],'--g')
            plot([0 256],0-[spots_to_match(1,CloseMatches(gg)) spots_to_match(1,CloseMatches(gg))],'--g')
            disp(sprintf('Spot at (%d,%d):',spots_to_match(1,CloseMatches(gg)),...
                spots_to_match(2,CloseMatches(gg))))
            match_I(gg,:) = CalcOneSpotIntensity(GrSpots(:,CloseMatches(gg)),match_vars(:,CloseMatches(gg)),...
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
            ismatch = input('Is this a match?(y/n)','s');
            close
            if strcmpi(ismatch,'y')
                if strcmpi(base_channel,'Red')
                    matchedR(:,end+1) = base_spots(:,ss);
                    RedI(end+1,:) = base_I;
                    matchedG(:,end+1) = GrSpots(:,CloseMatches(gg));
                    GrI(end+1,:) = match_I(gg,:);
                else
                    matchedG(:,end+1) = base_spots(:,ss);
                    GrI(end+1,:) = base_I;
                    matchedR(:,end+1) = spots_to_match(:,CloseMatches(gg));
                    RedI(end+1,:) = match_I(gg,:);
                end
                break
            end
        end
        if strcmpi(ismatch,'n')
            if strcmpi(base_channel,'Red')
                unmatchedR(:,ss) = base_spots(:,ss);
                unmatchedRedI(ss,:) = base_I;
            else
                unmatchedG(:,ss) = base_spots(:,ss);
                unmatchedGrI(ss,:) = base_I;
            end
        end
            
        %[closestcorr,ind] = min(spotcorrs);
        %if closestcorr<1.5
%             % These guys are matched
%             disp('Found a match!')
%             if strcmpi(base_channel,'Red')
%                 matchedR(:,end+1) = base_spots(:,ss);
%                 RedI(end+1,:) = base_I;
%                 matchedG(:,end+1) = spots_to_match(:,CloseMatches(ind));
%                 GrI(end+1,:) = match_I(CloseMatches(ind),:);
%             else
%                 matchedG(:,end+1) = base_spots(:,ss);
%                 GrI(end+1,:) = base_I;
%                 matchedR(:,end+1) = spots_to_match(:,CloseMatches(ind));
%                 RedI(end+1,:) = match_I(CloseMatches(ind),:);
%             end
%         else
%             disp('No matches found.')
%             % Keep the intensity calculation we already did anyway
%             if strcmpi(base_channel,'Red')
%                 unmatchedR(:,ss) = base_spots(:,ss);
%                 unmatchedRedI(ss,:) = base_I;
%                 unmatchedG(:,CloseMatches) = spots_to_match(:,CloseMatches);
%                 unmatchedGrI(CloseMatches,:) = match_I;
%             else
%                 unmatchedG(:,ss) = base_spots(:,ss);
%                 unmatchedGrI(ss,:) = base_I;
%                 unmatchedR(:,CloseMatches) = spots_to_match(:,CloseMatches);
%                 unmatchedRedI(CloseMatches,:) = match_I;
%             end 
%         end
    end
    
end

if exist('roughtform','var')
    matchedG = roughtform.FRETmapInv(matchedG);
    unmatchedG = roughtform.FRETmapInv(unmatchedG);
end

end