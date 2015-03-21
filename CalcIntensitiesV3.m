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
    numspots = size(Rspots,2);
elseif tform == 1
    Gspots = Rspots;
    Rspots = [];
    numspots = size(Gspots,2);
elseif tform == -1
    Gspots = [];
    numspots = size(Rspots,2);
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
    
    for kk = 1:numspots
        
       if params.IntensityGaussWeight==1

           if ~isempty(Rspots)
               [spotimgR,localcenR] = ExtractROI(imgR,params.DNASize,Rspots(:,kk));
               [spotRbkgnd,~] = ExtractROI(bkgndR,params.DNASize,Rspots(:,kk));
               RedI(kk,jj:framesdone) = CalcSpotIntensity('Gauss',...
                   spotimgR-spotRbkgnd,localcenR,spotVars(:,kk),params);
           end
           
           if ~isempty(Gspots)
               [spotimgG,localcenG] = ExtractROI(imgG,params.DNASize,Gspots(:,kk));
               [spotGbkgnd,~] = ExtractROI(bkgndG,params.DNASize,Gspots(:,kk));
               GrI(kk,jj:framesdone) = CalcSpotIntensity('Gauss',...
                   spotimgG-spotGbkgnd,localcenG,spotVars(:,kk),params);
           end

       else
           if ~isempty(Rspots)
               [spotimgR,localcenR] = ExtractROI(imgR,5,Rspots(:,kk));
               [spotRbkgnd,~] = ExtractROI(bkgndR,params.DNASize,Rspots(:,kk));
               RedI(kk,jj:framesdone) = CalcSpotIntensity('NoGauss',...
                   spotimgR-spotRbkgnd,localcenR,[],params);
           end
           % Get ROI in green channel:
           if ~isempty(Gspots)
               [spotimgG,localcenG] = ExtractROI(imgG,5,Gspots(:,kk));
               [spotGbkgnd,~] = ExtractROI(bkgndG,params.DNASize,Gspots(:,kk));
               GrI(kk,jj:framesdone) = CalcSpotIntensity('NoGauss',...
                   spotimgG-spotGbkgnd,localcenG,[],params);
           end
       end
       clear spotimgG spotimgR spotRbkgnd spotRbgknd
    end
    
   disp(sprintf('Calculated intensity for frames %d to %d of %d', jj, framesdone,totframes))
   clear imgR imgG bkgndR bkgndG framesdone
end
end
