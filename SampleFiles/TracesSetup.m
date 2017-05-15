% function TracesSetup()
%
% Called by Traces.m; sets up some basic parameters. This is the only code 
% that the user should have to change.
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

function TracesSetup()

%%%%%%%% Directory defaults: %%%%%%%%
% Where to save analyzed data:
defaultsavedir = '/Volumes/smFRET/smFRET data analysis';
% Where to load data from:
defaultdatadir = '/Volumes/smFRET/smFRET data';
% Where the code is (which is where it saves this parameter file)
codedir = '/Users/Steph/Documents/UCSF/Narlikar lab/smFRET analysis code';

%%%%%%%% Display defaults: %%%%%%%
% These control where the two figures that show traces and a field of view
% show up on your monitor. Numbers are [top left corner x position, top left y,
% width, height].  Run figure('Position',Fig2Pos) to see where the figure
% will show up; adjust the Fig2Pos vector elements as necessary.  Note: Macs
% and PC's have a different (0,0) for the monitor.
% Fig2Pos = [650,800,700,550];
Fig2Pos = [650,800,500,650];
% Fig1Pos = [0,400,600,500];
Fig1Pos = [25,400,600,600];

%%%%%%%% Microscope-specific parameters: %%%%%%%%
% Note: This code hasn't really been de-bugged for settings other than
% splitx = 1, Acceptor = 0.
splitx = 1; % If this is 1, red-green channels are left and right (not top and bottom)
Acceptor = 0; % If this is 1, the acceptor channel is the one on the right 
    % (or on the bottom if splitx is 0)
% (Our acquisition procedure automatically saves other necessary variables
% like how many pixels on each side our images our, the frame rate, etc.
% But those could be coded in here and then Traces.m rewritten to use values
% from the AnalysisParameters.mat file rather than from the metadata file.)  
    
%%%%%%%% Analysis parameters: %%%%%%%%
SmoothIntensities = 5; % If this is zero (or negative), don't do any smoothing 
    % of the acceptor and donor intensities; if greater than zero, median
    % filter of width specified by this variable.  Must be an integer.
SmoothFRET = 5; % Same as SmoothIntensities but for the FRET signal. 
EndInjectFrame = 1;% If this is bigger than 1, spotfinding will start after
    % this frame (instead of the first 1:FramesToAvg frames). (Relic from
    % when we did manual injections, which bumped the stage.)
DetectRedFlash = 0; % If this is greater than 0, Traces will look for a "flash"
    % in the acceptor channel, which we use to mark injection via an
    % automated syringe pump.
InjectDelay = 2.6; % If you are using an automated syringe pump to inject,
    % and you have measured the delay, put that information here.  This is 
    % only used if red laser flashes are detected (see above). Our delay is
    % 2.6+/0.3 sec.
ManualInjectMark = 11.8; %If this is greater than 0, UserSpotSelection 
    % will plot a vertical line at this number of seconds, and this value 
    % will be saved as the time of injection. This isn't used if
    % DetectRedFLash is nonzero and acceptor channel flashes are found.
FramesToAvg = 20; % How many frames to average over for spotfinding. 10-20 is a good value
    % for me.
FindSpotsEveryXFrames = 0; % If this is greater than 0,
    % spots will be found every this many frames (but still averaging over
    % FramesToAvg frames)
CheckSpotFindingEveryXFrames = 0; % If this is greater than zero, Traces will ask
    % the user to check the fidelity of the spotfinding threshhold every
    % this many frames. I recommend if FindSpotsEveryXFrames is greater
    % than 0, that this is set to something like 5*FindSpotsEveryXFrames.
UseCombinedImage = 0; % If this is 1, use an image of one (transformed) channel
    % overlaid on the other to find spots in real data. Otherwise, find
    % spots separately in each channel. While using a combined image has
    % the advantage of capturing mid-FRET spots more easily, it depends heavily 
    % on the quality of the channel mapping. NOTE this does not work very well
    % right now, depending on the fidelity of the transform. Note also that
    % this currently uses an affine transformation only, since polynomial
    % always does worse at creating a combined image. It's probably better
    % to just adjust the spotfinding threshold to capture mid-FRET spots.
    % Finally, only Matlab 2013 and later has imwarp, the overlay
    % functionality; if you don't have 2013 or later, you can't use this
    % functionality (you'll get a warning if you try and UseCombinedImage
    % will be automatically reset to 0).
TransformToCalc = 'MatlabPoly'; % Options are Affine, Poly, MatlabAffine, MatlabPoly
    % (caps insensitive; the Matlab* versions use built-in Matlab functions
    % instead of my hand-written code)
TformMaxDeg = 4; % If TransformToCalc is Poly or MatlabPoly, max degree of the polynomial
    % (note if using built-in Matlab functions, this should equal TformTotDeg)
TformTotDeg = 4; % If TransformToCalc is Poly or MatlabPoly, max degree of the polynomial
    % (note if using built-in Matlab functions, this must be 2, 3, or 4)
ResidTolerance = 0.008; % When calculating channel mapping: what's the maximum residual, 
    % divided by total number of spots, allowable.
Refine_Bd_Cen = 1; % If this is 1, use a 2D gaussian fit to refine the bead center
    % position.  Highly recommended.
IntensityGaussWeight = 1; % If this is 1, weight the intensity of each spot  
    % in each frame by a Gaussian whose center (and possibly variance) are determined
    % from a fit to the spot's first FramesToAvg frames. Note that if this is 0, it
    % will calculate intensities over a 5 pxl diameter disk.  That's
    % hard-wired into the code--see CalcIntensitiesV3.m to change
    % the size. Obviously it's best to do the Gaussian weighting, but it
    % works surprsingly well not to, if that's preferable. (There's not 
    % really a difference in speed or anything.)
GaussWeightAmp = 2; % If IntensityGaussWeight is 1, this determines the amplitude
    % of the Gaussian used to weight each spot's intensity. Effectively
    % this just scales the intensity values--the Ha lab code uses an
    % amplitude of 2, probably so that changes in intensity are more
    % noticeable. I haven't noticed much of a difference bewteeen 1 vs 2.
FixSpotVar = []; % If IntensityGaussWeight is 1, and FixSpotVar is not
    % the empty matrix, all spot intensities will be weighted by a Gaussian
    % with [xvar;yvar] = 1./(2.*FixSpotVar). The Ha lab always uses a fixed
    % FixSpotVar = [0.3; 0.3]; I haven't found it to make a huge
    % difference, but 0.3 seems to be good.  Most of my spots have values
    % between 0.3 and 1, with an average at about 0.7. Err on the smaller
    % side for these values, as that will mean you don't run the risk of
    % under-weighting real intensity.
UseSymGauss = 0; % If this is 1, insist that the Gaussian used if IntensityGaussWeight=1
    % is symmetric (same variances in x and y). Does not affect the use of
    % Gaussian fitting in spotfinding. Given the distortions at the bottom
    % of our images I do not use a symmetric Gaussian (though of course this parameter
    % does not matter if you use FixSpotVar with identical x and y values).
BeadSize = 8; % Diameter of a circle that defines a bead (used for the channel
    % mapping procedure); beads whose centers are closer than BeadSize will 
    % not be included, and found beads will be circled by a circle of radius BeadSize.  
BeadNeighborhood = 9^2; % Our spot-finding algorithm looks for local maxima in local
    % "neighborhoods". This defines the size of a neighborhood (area, in square
    % pixels) for the beads.  Needs to be a perfect square, and best if
    % sqrt(BeadNeighborhood) is odd.  Should be a little bigger than we
    % expect beads to be.
DNASize = 8; % Same as BeadSize but for DNA: diameter of expected spots.  Note that
    % if IntensityGaussWeight=1, this is also the side of a square over which
    % a Gaussian is fit and the intensity calculated. However, if IntensityGaussWeight=0,
    % the intensity is summed over a 5-pixel diameter circle and this parameter
    % has no effect.  In both cases DNASize also determines how close two 
    % spots can be and still be included in the analysis.
    % I have found that 6 or 8 is a good number.
DNANeighborhood = 9^2; % Same as BeadNeighborhood but for DNA.
alpha = 0.1; % Crosstalk between channels: Corrects for bleed-through of donor intensity
    % into acceptor channel.  Corrects raw acceptor intensities I_A,raw
    % according to the formula I_A = I_A,raw - alpha*I_D, where I_D is the
    % donor intensity. Set to 0 to not correct for channel cross-talk.
    % On our setup, alpha should be about 0.1 at ~10 mW laser power, maybe
    % as high as 0.14 for 30 mW laser power.
    % TODO: Not clear to me if this should be done before or after
    % background subtraction? Doing it after, consistent with Ha lab IDL
    % code.
gamma = 1; % Detection inequality (and quantum yield inequalities, etc) between dyes. 
    % Corrects FRET values according to FRET = I_A / (I_D + gamma*I_A). Set
    % to 1 to not correct for detection inequality (Ha and Selvin 2007 say
    % for Cy3 and Cy5, gamma is roughly 1).
PxlsToExclude = 10; % How many pixels on each side of the image, along the axis that
    % contains both channels, to exclude from analysis.  On our system with
    % a decent channel alignment this is about 10 pixels.  This avoids
    % finding spots in areas of the image where the channels might overlap.
    % Set to zero to use the whole image. NOTE: This parameter MUST be the
    % same for the data used to create a channel map, as for any data you
    % analyze with that map. Traces will insist on using the value for
    % PxlsToExclude that comes with any map you load. Re-do a map with a
    % different PxlsToExclude value in order to change the pixels excluded
    % with real data!
FrameLoadMax = 500; % Maximum number of frames to load into Matlab's memory 
    % at one time. Making this number bigger will make some parts of the
    % analysis process marginally faster, but not hugely significantly so, at
    % least on my laptop.  Note: you can change this number as much
    % as you like, and the rest of the analysis suite will run fine 
    % (for example, you can scale movies with it set to 500, and then change 
    % it, making it either bigger or smaller, and calculate intensities, etc.).
    % It will always use whatever the current value of FrameLoadMax is (so
    % if, for example, you scaled movies in 500-frame chunks, but then change 
    % this to 100 and rerun the same data set, the play movie feature in
    % the GUI will only load 100 frames at a time).
ScaleChannelsSeparately = 1; % If this is 0, ScaleMovieV2 will scale every intensity
    % in every frame between the maximum and minimum values across the
    % entire movie.  If this is 1, it will scale the donor channel to the
    % max and min of the donor side only, and the acceptor to the max and
    % min of the acceptor side only.  I believe the convention in the field
    % has been to scale both channels together, but it appears for our
    % system that that will lead to donor intensities that are consistently
    % half those of the acceptor channel.
NormImage = 0; % If this is 1, ScaleMovieV2 will normalize each pixel's intensity, 
    % in each frame, to the median intensity of the (512x512) image at that
    % frame. We've been observing large fluctuations in total image
    % intensity over time, which may be due to laser power fluctuations;
    % this is an attempt to correct for that.
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PLEASE DO NOT CHANGE CODE BELOW HERE %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Check parameters are reasonable %%%%%%%%%%%%%%%%%%%%
    PxlsToExclude = round(PxlsToExclude);
    EndInjectFrame = round(EndInjectFrame);
    FindSpotsEveryXFrames = round(FindSpotsEveryXFrames);
    if EndInjectFrame<=0
        EndInjectFrame = 1;
    end
    if InjectDelay<0
        InjectDelay=0;
    end
    if FindSpotsEveryXFrames<0
        FindSpotsEveryXFrames = 0;
    end
    if ~isempty(FixSpotVar) && (FixSpotVar(1)<=0 || ...
            FixSpotVar(2)<=0)
        FixSpotVar = [];
    end
    MatlabVer = ver;
    MatlabDate = MatlabVer(1).Date;
    if UseCombinedImage == 1 && str2double(MatlabDate(end-1:end))<=11 % Testing for Matlab versions older than 2012
        disp('Warning: This version of Matlab does not support creation of a combined image.')
        UseCombinedImage = 0;
    end
    clear MatlabVer MatlabDate

%%%%%%%%%%%%% Save the paramters %%%%%%%%%%%%%%%%%%%%
save(fullfile(codedir,'AnalysisParameters.mat'),'defaultsavedir',...
    'defaultdatadir','splitx','Acceptor','BeadSize','BeadNeighborhood',...
    'DNASize','DNANeighborhood','SmoothIntensities','SmoothFRET',...
    'Fig1Pos','Fig2Pos','FramesToAvg','PxlsToExclude','Refine_Bd_Cen',...
    'FrameLoadMax','UseCombinedImage','IntensityGaussWeight','NormImage',...
    'TransformToCalc','TformMaxDeg','TformTotDeg','ResidTolerance',...
    'UseSymGauss','EndInjectFrame','DetectRedFlash','InjectDelay',...
    'ManualInjectMark','FindSpotsEveryXFrames','alpha','gamma',...
    'CheckSpotFindingEveryXFrames','GaussWeightAmp','FixSpotVar',...
    'ScaleChannelsSeparately');
