% function [matched1,matched2] = FindSpotMatches(spots1,spots2)
%
% Given two sets of spots (e.g., a set of spots found in the red channel and
% a set of spots in the green channel), find a matching between them.
% IMPORTANT: it will match the spots in spots1 to spots in spots2, so the
% order you enter the inputs in matters!
%
% Inputs are (x,y) pairs along the rows or columns for two sets of spots.
% Returns two sets of spots, where each column of matched1 and matched 2 is
% (x,y) of a spot, and the first column of matched1 is the spot that matches
% the spot represented by the first column of matched 2, etc.
%
% Currently implemented with a "greedy" algorithm, which assumes the
% alignment between the spots (ie. the channels) isn't too far off to begin
% with.  So just find the spot in the other channel with the shortest
% distance to the current spot.
%
% Steph 9/2013
% Copyright 2013 Stephanie Johnson, University of California, San Francisco

function [matched1,matched2] = FindSpotMatches(spots1,spots2)

% global TotImg

% Make the inputs 2-by-numspots matrices, if they're not already, for
% convenience
if size(spots1,1)~=2
    spots1 = transpose(spots1);
end
if size(spots2,1)~=2
    spots2 = transpose(spots2);
end

% Get the distances between all the spots:
Dists = FindSpotDists(spots1,spots2);

[minDists,ind] = min(Dists,[],2); %Minimum along each row
% Assume that if a spot in spots1 doesn't have match that's a shorter
% distance away than the mean distance over all the spots plus the standard
% deviation, that that spot doesn't have a match:
maxDist = mean(minDists)+std(minDists,1);

matched1 = [];
matched2 = [];
matches = 1;
taken = []; %This prevents one spot in spots2 from getting matched to two 
    %spots in spots1
for i=1:size(Dists,1)
    if minDists(i)<=maxDist 
        if sum(ismember(taken,ind(i)))==0
            matched1(:,matches) = spots1(:,i);
            matched2(:,matches) = spots2(:,ind(i));
            matches = matches+1;
            taken(end+1) = ind(i);
%             %For debugging: need to uncomment the global TotImg at the top
%             matchG_abs(1,:) = matched2(1,:);
%             matchG_abs(2,:) = matched2(2,:)+size(TotImg,2)/2;
%             figure
%             PutBoxesOnImageV3(mat2gray(mean(TotImg(:,:,1:10),3)),[matched1';matchG_abs'],9);
%             pause
%             close
%             clear matchG_abs
        else %Another spot in spots1 has already claimed this spot in spots2
                %as its nearest neighbor.  Which one is closest?
            samespots2 = find(ind==ind(i));
            [~,winner] = min(minDists(samespots2));
            %samespots2(winner) is the index of the spot in spots1 that has the
            %shortest distance to the spot2 that has multiple matches.
            if samespots2(winner)==i
                %Need to remove the previous spot in spot1 that had claimed
                %this spot2.  Since this process is done iteratively, there
                %should be only one spot before i that has the same ind as
                %spot i.  So the spot to remove is spot1(samespots2(1)).
                %But that's not necessarily the same as
                %matched1(samespots2(1)).
                %We know (i-matches-1) spot1's have been skipped because they didn't
                %have a minDist to a spot2 that was less than maxDist:
                matched1temp = matched1;
                matched2temp = matched2;
                clear matched1 matched2
                if (samespots2(1)-(i-matches-1)) == 0
                    matched1 = matched1temp(:,2:end);
                    matched2 = matched2temp(:,2:end);
                elseif (samespots2(1)-(i-matches-1)+2) > length(matched1temp)
                    matched1 = matched1temp(:,1:end-1);
                    matched2 = matched2temp(:,1:end-1);
                else
                    matched1 = matched1temp(:,1:(samespots2(1)-(i-matches-1)));
                    try
                    matched1 = [matched1, matched1temp(:,(samespots2(1)-(i-matches-1)+2):end)];
                    catch
                        keyboard
                    end
                    matched2 = matched2temp(:,1:(samespots2(1)-(i-matches-1)));
                    matched2 = [matched2, matched2temp(:,(samespots2(1)-(i-matches-1)+2):end)];
                end
                matched1(:,matches) = spots1(:,i);
                matched2(:,matches) = spots2(:,ind(i));
                clear matched1temp matched2temp
            %if samespots2(winner)<i, we've already put it in matched1; if
            %samespots2(winner)>i,we'll deal with it when we get to it.  in either
            %case, do nothing with this spot for now.
%             %For debugging: need to uncomment the global TotImg at the top
%             matchG_abs(1,:) = matched2(1,:);
%             matchG_abs(2,:) = matched2(2,:)+size(TotImg,2)/2;
%             figure
%             PutBoxesOnImageV3(mat2gray(mean(TotImg(:,:,1:10),3)),[matched1';matchG_abs'],9);
%             pause
%             close
%             clear matchG_abs
            end
        end
    end
end
disp(strcat('Number of spots paired: ',int2str(matches-1)))