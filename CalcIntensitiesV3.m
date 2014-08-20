% function [RedI,GrI] = CalcIntensitiesV3(PathToMovie, Rspots, spotVars,...
%       tform,params)
%
% Calculates intensities for all spots in a movie.
%
% Inputs:
% PathToMovie: a path to a directory with the movie
% Rspots: locations of the spots in the acceptor channel
% spotVars: x,y variances to use for each spot for the Gaussian weighting
% tform: mapping information for finding spots in the donor channel. If you
%   only want to calculate the intensities for the acceptor channel, pass -1 for
%   tform; if you want to only calculate intensities for the donor, pass 1.
%   Otherwise tform must be a FRETmap class object.
% params: file saved by smFRETsetup
%
% Outputs:
% RedI, GrI: Intensity-vs-time information for each spot
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

function [RedI, GrI] = CalcIntensitiesV3(PathToMovie, Rspots, spotVars, ...
    tform,params)

% Figure out the total number of image files in this movie:
[~,totframes] = LoadRawImgs(PathToMovie,'FramesToLoad',[1 1]);

if params.FrameLoadMax > totframes
    FrameLoadMax = totframes;
else
    FrameLoadMax = params.FrameLoadMax;
end

RedI = zeros(size(Rspots,2),totframes);
GrI = zeros(size(Rspots,2),totframes);
% Find the spots in the coordinate system of the other (green) channel:
if isa(tform,'FRETmap')
    Gspots = tform.FRETmapInv(Rspots);
elseif tform == 1
    Gspots = Rspots;
    Rspots = [];
elseif tform == -1
    Gspots = [];
end    

for jj = 1:FrameLoadMax:totframes
    
    if tform == 1
        [~,imgG,~,bkgndG] = LoadScaledMovie(PathToMovie,[jj jj+FrameLoadMax-1],params,'gbkgnd');
        framesdone = jj+size(imgG,3)-1;
    elseif tform == -1
        [imgR,~,bkgndR,~] = LoadScaledMovie(PathToMovie,[jj jj+FrameLoadMax-1],params,'rbkgnd');
        framesdone = jj+size(imgR,3)-1;
    else
        [imgR,imgG,bkgndR,bkgndG] = LoadScaledMovie(PathToMovie,[jj jj+FrameLoadMax-1],params,'bkgnd');
        framesdone = jj+size(imgR,3)-1;
    end
    
    for kk = 1:size(Rspots,2)
        
       if params.IntensityGaussWeight==1

           if ~isempty(Rspots)
               [spotimgR,localcenR] = ExtractROI(imgR,params.DNASize,Rspots(:,kk));
               [spotRbkgnd,~] = ExtractROI(bkgndR,params.DNASize,Rspots(:,kk));
               RedI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensity('Gauss',...
                   spotimgR-spotRbkgnd,localcenR,spotVars(:,kk),params);
           end
           
           if ~isempty(Gspots)
               [spotimgG,localcenG] = ExtractROI(imgG,params.DNASize,Gspots(:,kk));
               [spotGbkgnd,~] = ExtractROI(bkgndG,params.DNASize,Rspots(:,kk));
               GrI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensity('Gauss',...
                   spotimgG-spotGbkgnd,localcenG,spotVars(:,kk),params);
           end

       else
           if ~isempty(Rspots)
               [spotimgR,localcenR] = ExtractROI(imgR,5,Rspots(:,kk));
               [spotRbkgnd,~] = ExtractROI(bkgndR,params.DNASize,Rspots(:,kk));
               RedI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensity('NoGauss',...
                   spotimgR-spotRbkgnd,localcenR,[],params);
           end
           % Get ROI in green channel:
           if ~isempty(Gspots)
               [spotimgG,localcenG] = ExtractROI(imgG,5,Gspots(:,kk));
               [spotGbkgnd,~] = ExtractROI(bkgndG,params.DNASize,Rspots(:,kk));
               GrI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensity('NoGauss',...
                   spotimgG-spotGbkgnd,localcenG,[],params);
           end
       end
       clear spotimgG spotimgR spotRbkgnd spotRbgknd
    end
    
   disp(sprintf('Calculated intensity for frames %d to %d of %d', jj, framesdone,totframes))
   clear imgR imgG bkgndR bkgndG framesdone
end
end
