% function ComputeBackgroundImgs(PathToMovie,params)
%
% Steph 8/2014
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

function ComputeBackgroundImgs(PathToMovie,params)

% Call LoadRawImgs once to get number of frames:
[~,numframes] = LoadRawImgs(PathToMovie,'FramesToLoad',[1 1]);

if params.FrameLoadMax > numframes
    FrameLoadMax = numframes;
else
    FrameLoadMax = params.FrameLoadMax;
end

disp('Calculating background ...')
    
% Re-load everything and actually do the scaling:
% Update 8/2014: I'm not going to save the scaled version, because it's
% actually faster to load tifs and scale on the fly.  However, I do
% need to scale before I calculate background images.
for jj = 1:FrameLoadMax:numframes
    [imgR,imgG] = LoadScaledMovie(PathToMovie,...
        [jj jj+FrameLoadMax-1],params);

    % Actual background calculation is done by a separate function:
    [imgRBkgnd,imgGBkgnd] = CalcBkgnd(imgR,imgG,params);

    save(fullfile(PathToMovie,strcat('BackgroundImgs',int2str(jj),...
        'to',int2str(jj+FrameLoadMax-1),'.mat')),'imgRBkgnd','imgGBkgnd')

    clear imgR imgG imgRBkgnd imgGBkgnd

end
ScaleChannelsSeparately = params.ScaleChannelsSeparately;
save(fullfile(PathToMovie,strcat('ScalingInfo.mat')),'ScaleChannelsSeparately','-append')
clear ScaleChannelsSeparately