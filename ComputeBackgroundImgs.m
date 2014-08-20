% function ComputeBackgroundImgs(PathToMovie,params)
%
% Steph 8/2014
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