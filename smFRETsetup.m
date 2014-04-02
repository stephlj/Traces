% function smFRETsetup
%
% Called by smFRET; sets up some basic parameters. This is the only code 
% that a new user might have to change.
%
% Steph 9/2013
% Copyright 2013 Stephanie Johnson, University of California, San Francisco

function smFRETsetup

%%%%%%%% Directory defaults: %%%%%%%%
% Where to save analyzed data:
% defaultsavedir = '/Volumes/smFRET/smFRET data analysis';
% defaultsavedir = '/Users/Steph/Documents/UCSF/Narlikar lab/smFRET analysis code/SampleData/RealDNA/analysis';
defaultsavedir = '/Users/Steph/Dropbox/Steph Dropbox/Narlikar Lab DB/smFRET data analysis';
% Where to load data from:
defaultdatadir = '/Users/Steph/Documents/UCSF/Narlikar lab/smFRET analysis code/SampleData';
% defaultdatadir = '/Users/Steph/Dropbox/Steph Dropbox/Narlikar Lab DB/smFRET data';
% Where the code is (which is where it saves this parameter file)
% codedir = '/Volumes/smFRET/smFRET analysis code';
codedir = '/Users/Steph/Documents/UCSF/Narlikar lab/smFRET analysis code';

%%%%%%%% Display defaults: %%%%%%%
% These control where the two figures that show traces and a field of view
% show up on your monitor. Numbers are [top left corner x position, top left y,
% width, height].  Run figure('Position',Fig2Pos) to see where the figure
% will show up; adjust the Fig2Pos vector elements as necessary.  Note: Macs
% and PC's have a different (0,0) for the monitor I think.
% Fig2Pos = [650,800,700,550];
Fig2Pos = [650,800,500,650];
% Fig1Pos = [0,400,600,500];
Fig1Pos = [25,400,600,500];

%%%%%%%% Microscope-specific parameters: %%%%%%%%
% Note: This code hasn't really been de-bugged for settings other than
% splitx = 1, Acceptor = 0.
splitx = 1; % This means red-green channels are left and right (not top and bottom)
Acceptor = 0; % This means the acceptor channel is the one on the left 
    % (or on the top if splitx is 0)
% Our acquisition procedure automatically saves other necessary variables
% like how many pixels on each side our images our, the frame rate, etc.
% But those could be coded in here and then smFRET rewritten to use values
% from smFRETsetup rather than from the acquisition file.
    
%%%%%%%% Analysis parameters: %%%%%%%%
FramesToAvg = 10; % How many frames to average over for spotfinding and calculating
    % the local background that will be subtracted
PxlsToExclude = 10; % How many pixels on each side of the image, along the axis that
    % contains both channels, to exclude from analysis.  On our system with
    % a decent channel alignment this is about 10 pixels.  This avoids
    % finding spots in areas of the image where the channels might overlap.
    % Set to zero to use the whole image.
UseCombinedImage = 0; % If this is 1, use an image of one (transformed) channel
    % overlaid on the other to find spots in real data. Otherwise, find
    % spots separately in each channel. While using a combined image has
    % the advantage of capturing mid-FRET spots, it depends heavily on the
    % quality of the channel mapping.
Refine_Bd_Cen = 1; % If this is 1, use a 2D gaussian fit to refine the bead center
    % position.  This will increase computational time for the channel
    % mapping by about a factor of 2, for roughly a factor of 2 improvement
    % in bead center localization (that is, the mean error will go down by
    % a factor of 2).
IntensityGaussWeight = 0; % If this is 1, weight the intensity of each spot  
    % in each frame by a Gaussian whose center and variance are determined
    % from a fit to the spot's first 10 frames. Note that if this is 0, it
    % will calculate intensities over a 5 pxl diameter circle.  That's
    % hard-wired into the code--see CalcSpotIntensityNoGauss.m to change
    % the size.
BeadSize = 8; % Diameter of a circle that defines a bead (used for the channel
    % mapping procedure); beads whose centers are closer than BeadSize will 
    % not be included, and found beads will be circled by a circle of radius BeadSize.  
BeadNeighborhood = 9^2; % Our spot-finding algorithm looks for local maxima in 
    % "neighborhoods" (see FindSpotsV4 comments). This defines the size of a
    % neighborhood (area, in square pixels) for the beads.  Needs to be a perfect
    % square, and best if sqrt(BeadNeighborhood) is odd.  Should be a little 
    % bigger than we expect beads to be.
DNASize = 8; % Same as BeadSize but for DNA: diameter of expected spots.  Note that
    % the code that actually calculates the intensity of a DNA spot is
    % currently hard-coded to integrate over a circle of diameter 5 pixels;
    % DNASize only determines how close two spots can be and still be
    % included in the analysis.
DNANeighborhood = 9^2; % Same as BeadNeighborhood but for DNA.
BkgndSubSigma = 4; % For background subtraction: variance of the Gaussian filter that is applied
SmoothIntensities = 0; % If this is zero (or negative), don't do any smoothing 
    % of the acceptor and donor intensities; if greater than zero, moving
    % average smoothing filter of width specified by this variable.  Must be
    % an integer.
SmoothFRET = 0; % Same as SmoothIntensities but for the FRET signal.  At some
    % point should implement a Gauss filter instead
    
%%%%%%%%%%%%% Save the paramters %%%%%%%%%%%%%%%%%%%%
save(fullfile(codedir,'AnalysisParameters.mat'),'defaultsavedir',...
    'defaultdatadir','splitx','Acceptor','BeadSize','BeadNeighborhood',...
    'DNASize','DNANeighborhood','SmoothIntensities','SmoothFRET',...
    'Fig1Pos','Fig2Pos','FramesToAvg','PxlsToExclude','Refine_Bd_Cen',...
    'BkgndSubSigma','UseCombinedImage','IntensityGaussWeight');