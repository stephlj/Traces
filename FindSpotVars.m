% function Vars = FindSpotVars(imgs,spots,params,debug)
%
% Given a time series (or one frame) of a movie and the locations of spots
% in that movie, fit a 2D Gaussian to each spot to obtain the variance of
% each spot. This will then be used in CalcIntensities to weight a spot's
% intensity by a Gaussian. Note that this function is called by smFRET
% after spots have already been found so the (x,y) center of the spot is
% taken to be fixed and only the variances, amplitude and background are
% fit for each spot.
%
% The assumption is that this is called separately on spots in each
% channel.
%
% Inputs:
% imgs: x-by-y or x-by-y-by-frame total fields of view in which spots have
%   already been located
% spots: (x,y)-by-number locations of spots in imgs, which have already
%   been fit to a Gaussin in FindSpotsV5 to refine their center location
% params: file saved by smFRETsetup
% debug: If this is 1, plots the figures in Fit2DGaussToSpot. Default is 0.
%
% Outputs:
% Vars: (xvar,yvar)-by-spots variances, which are the mean values of these
%   variances across however many frames are in imgs
% optional output: amplitudes from the fits, in case that information is
%   used to choose which set of variances (red channel vs green channel) to
%   keep in the function that calls this one.
%
% Stephanie 5/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [Vars,varargout] = FindSpotVars(imgs,spots,params,debug)

if ~exist('debug','var') debug = 0; end

ROIsize = params.DNASize;
symGauss = params.UseSymGauss;

defaultXvar = 0.3; %If the fit fails, set the spot variance to this value
    % For now taking this value from the Ha lab IDL code, but it's probably
    % system (camera, microscope, etc) specific
defaultYvar = defaultXvar;

Vars = zeros(2,size(spots,2),size(imgs,3)); % Vars will be:
    % First row will contain xvars, second row yvars
    % Columns will be spots
    % Third dimension is frame
Amps = zeros(size(spots,2),size(imgs,3));

for i = 1:size(imgs,3)
    for ss = 1:size(spots,2)
       % Extract an ROI over which to calculate the fit:
       [spotimg,localcen] = ExtractROI(imgs(:,:,i),ROIsize,spots(:,ss));
       % Note if the fit fails, Fit2DGaussToSpot will return the
       % StartParams, so make sure the StartParams for xvar and yvar are
       % the defaults hardcoded above
       [~,~,Vars(1,ss,i),Vars(2,ss,i),~,Amps(ss,i)] = Fit2DGaussToSpot(spotimg,'Vars',...
           'StartParams',[localcen(1),localcen(2),defaultXvar,defaultYvar,min(spotimg(:)),max(spotimg(:))],...
           'Debug',debug,'symGauss',symGauss);
       clear xvar yvar A
    end
end

% Not sure yet if I want to do this:
Vars = mean(Vars,3);
varargout = mean(Amps,3);