% function Dists = FindSpotDists(spots1,spots2)
%
% Find the pairwise distances between all the spots in spots1 and all
% the spots in spots2 (or spots1 with itself, if no spots2 is passed). Spots
% 1 and spots2 must be (x,y) pairs of positions along either the rows or the columns.
%
% Returns a matrix of distances of size length(spots1) by length(spots2),
% where each row is the distance from each spot in spots1 to each spot in
% spots2.
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

function Dists = FindSpotDists(spots1,spots2)

if ~exist('spots2','var') spots2 = spots1; end

% Make the inputs 2-by-numspots matrices, if they're not already:
if size(spots1,1)~=2
    spots1 = transpose(spots1);
end
if size(spots2,1)~=2
    spots2 = transpose(spots2);
end

Dists = zeros(size(spots1,2),size(spots2,2));

% There are two ways to find the distances.  Here's the simplest but slower way:
% for i = 1:size(spots1,2)
%     for j = 1:size(spots2,2)
%         %Dists(i,j) = sqrt((spots1(1,i)-spots2(1,j))^2 + (spots1(2,i)-spots2(2,j))^2);
%         %Using the built-in Matlab function norm is faster
%         Dists(i,j) = norm(spots1(:,i)-spots2(:,j));
%     end
% end

% Smarter way: (Thanks to Matt Johnson for this derivation and code!)
% Consider two spots vec(a) = spots1(:,i) and vec(b) = spots2(:,j).  We want
% norm(vec(a)-vec(b)). Note
% norm(vec(a)-vec(b))^2 = <vec(a)-vec(b),vec(a)-vec(b)>
%                      = norm(vec(a))^2+norm(vec(b))^2 - 2<vec(a),vec(b)>
%                      = norm(vec(a))^2+norm(vec(b))^2 - 2*(transpose(spots1)*spots2)(i,j)
% Instead of having a double for-loop, you do this in one line as a matrix
% operation, using matlab's bsxfun to take care of the fact that spots1 and
% spots2 aren't necessarily the same size, and that transpose(spots1)*spots2
% isn't the same size as spots1 and spots2

Dists = sqrt(bsxfun(@plus,sum(spots1.^2,1)',sum(spots2.^2,1))-2*spots1'*spots2);
% What the above line does:
% sum(spots1.^2,1)' will be a size(spots1,2) column that contains the norm^2
%   of each column of spots1
% sum(spots2.^2,1) will be a size(spots2,2) row that contains the norm^2 of
%   each column of spots2
% bsxfun(@plus,sum(spots1.^2,1)',sum(spots2.^2,1)) will be a size(spots1,2)
%   by size(spots2,2) matrix where each element (i,j) is
%   norm(spots1(:,i))^2+norm(spots2(:,j))^2