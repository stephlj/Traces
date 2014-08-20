% function [RedI,GrI] = CalcIntensitiesV3(PathToMovie, Rspots, spotVars,...
%       tform,params)
%
% Calculates intensities for all spots in a movie.
%
% Inputs:
% PathToMovie: a path to a directory with the movie
% Rspots: locations of the spots in the acceptor channel
% spotVars: x,y variances to use for each spot for the Gaussian weighting
% tform: mapping information for finding spots in the donor channel
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
if ~isempty(tform)
    Gspots = tform.FRETmapInv(Rspots);
else
    Gspots = Rspots;
end    

for jj = 1:FrameLoadMax:totframes
    [imgR,imgG,bkgndR,bkgndG] = LoadScaledMovie(PathToMovie,[jj jj+FrameLoadMax-1],params,'bkgnd');
    
    for kk = 1:size(Rspots,2)
        
       if params.IntensityGaussWeight==1
            % Get ROI in red channel
           [spotimgR,localcenR] = ExtractROI(imgR,params.DNASize,Rspots(:,kk));
           [spotRbkgnd,~] = ExtractROI(bkgndR,params.DNASize,Rspots(:,kk));
           % Get ROI in green channel:
           [spotimgG,localcenG] = ExtractROI(imgG,params.DNASize,Gspots(:,kk));
           [spotGbkgnd,~] = ExtractROI(bkgndG,params.DNASize,Rspots(:,kk));
           
           RedI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensity('Gauss',...
               spotimgR-spotRbkgnd,localcenR,spotVars(:,kk),params);
           GrI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensity('Gauss',...
               spotimgG-spotGbkgnd,localcenG,spotVars(:,kk),params);

       else
           % Get ROI in red channel
           [spotimgR,localcenR] = ExtractROI(imgR,5,Rspots(:,kk));
           % Get ROI in green channel:
           [spotimgG,localcenG] = ExtractROI(imgG,5,Gspots(:,kk));

           RedI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensity('NoGauss',...
               spotimgR,localcenR,[],params);
           GrI(kk,jj:jj+size(imgR,3)-1) = CalcSpotIntensity('NoGauss',...
               spotimgG,localcenG,[],params);
       end
       clear spotimgG spotimgR
    end
   
   disp(sprintf('Calculated intensity for frames %d to %d of %d', jj, jj+size(imgR,3)-1,totframes))
   clear imgR imgG bkgndR bkgndG
end
end
