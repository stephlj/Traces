% function [newtestspots,newgoodspots,newtestVars,newgoodVars] = CheckSpotBoundaries(testspots,...
%    goodspots,testVars,goodVars,params,PathToMovie)
%
% It can happen that some spots, when transformed to the other channel, are
% too close to the image boundaries to be useable. This function checks all
% spots in testspots for being too close to the edge, where too close is defined
% by floor(boxdim/2); it then removes any such spots from both testspots and 
% goodspots, where it is assumed that testspots and goodspots are the same list
% of spots but in the two channels  It also removes the variances for any spots
% too close to the edge from testVars and goodVars, so that the outputs
% are all matrices of the same size. 
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

function [newtestspots,newgoodspots,newtestVars,newgoodVars] = CheckSpotBoundaries(testspots,...
    goodspots,testVars,goodVars,params,PathToMovie)

imgsize = GetInfoFromMetaData(PathToMovie,'imgsize');
if params.splitx
    imgsize(2) = imgsize(2)/2-2*params.PxlsToExclude;
else
    imgsize(1) = imgsize(1)/2-2*params.PxlsToExclude;
end

boxdim = params.DNASize;

% Check that the testspots aren't too close; remove any from both testspots and
% goodspots that are too close in the green channel
if length(find(ceil(testspots(1,:))>=1+floor(boxdim/2)))~=length(testspots(1,:))
    oldGspots = testspots;
    oldRspots = goodspots;
    oldGvars = testVars;
    oldRvars = goodVars;
    clear testspots goodspots goodVars testVars
    testspots = oldGspots(:,ceil(oldGspots(1,:))>=1+floor(boxdim/2));
    goodspots = oldRspots(:,ceil(oldGspots(1,:))>=1+floor(boxdim/2));
    testVars = oldGvars(:,ceil(oldGspots(1,:))>=1+floor(boxdim/2));
    goodVars = oldRvars(:,ceil(oldGspots(1,:))>=1+floor(boxdim/2));
    clear oldGspots oldRspots oldGvars oldRvars
end
if length(find(ceil(testspots(2,:))>=1+floor(boxdim/2)))~=length(testspots(2,:))
    oldGspots = testspots;
    oldRspots = goodspots;
    oldGvars = testVars;
    oldRvars = goodVars;
    clear testspots goodspots goodVars testVars
    testspots = oldGspots(:,ceil(oldGspots(2,:))>=1+floor(boxdim/2));
    goodspots = oldRspots(:,ceil(oldGspots(2,:))>=1+floor(boxdim/2));
    testVars = oldGvars(:,ceil(oldGspots(2,:))>=1+floor(boxdim/2));
    goodVars = oldRvars(:,ceil(oldGspots(2,:))>=1+floor(boxdim/2));
    clear oldGspots oldRspots oldGvars oldRvars
end
if length(find(ceil(testspots(1,:))<=imgsize(1)-floor(boxdim/2)))~=length(testspots(1,:))
    oldGspots = testspots;
    oldRspots = goodspots;
    oldGvars = testVars;
    oldRvars = goodVars;
    clear testspots goodspots goodVars testVars
    testspots = oldGspots(:,ceil(oldGspots(1,:))<=imgsize(1)-floor(boxdim/2));
    goodspots = oldRspots(:,ceil(oldGspots(1,:))<=imgsize(1)-floor(boxdim/2));
    testVars = oldGvars(:,ceil(oldGspots(1,:))<=imgsize(1)-floor(boxdim/2));
    goodVars = oldRvars(:,ceil(oldGspots(1,:))<=imgsize(1)-floor(boxdim/2));
    clear oldGspots oldRspots oldGvars oldRvars
end
if length(find(ceil(testspots(2,:))<=imgsize(2)-floor(boxdim/2)))~=length(testspots(2,:))
    oldGspots = testspots;
    oldRspots = goodspots;
    oldGvars = testVars;
    oldRvars = goodVars;
    clear testspots goodspots goodVars testVars
    testspots = oldGspots(:,ceil(oldGspots(2,:))<=imgsize(2)-floor(boxdim/2));
    goodspots = oldRspots(:,ceil(oldGspots(2,:))<=imgsize(2)-floor(boxdim/2));
    testVars = oldGvars(:,ceil(oldGspots(2,:))<=imgsize(2)-floor(boxdim/2));
    goodVars = oldRvars(:,ceil(oldGspots(2,:))<=imgsize(2)-floor(boxdim/2));
    clear oldGspots oldRspots oldGvars oldRvars
end

newgoodspots = goodspots;
newtestspots = testspots;
newgoodVars = goodVars;
newtestVars = testVars;