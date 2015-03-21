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