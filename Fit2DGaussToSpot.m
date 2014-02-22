% function [Xcen, Ycen, Xvar, Yvar, bkgnd, A] = Fit2DGaussToSpot(spotimg,varargin)
%
% Given a region-of-interest (i.e. part of a frame from a FRET movie) that
% nominally contains a single fluorescent spot, find the parameters of the
% 2D Gaussian that best fits that spot.
%
% Input:
% spotimg: image of a single spot
%
% Optional inputs: enter these as pairs ('<paraname>',<value>)
% 'Debug', 1: display a set of images of the fit
% 'symGauss', 1: force the variances in x and y to be the same
%
% Outputs: best fit values for:
% Xcen, Ycen: location of the center of the spot in x and y
% Xvar, Yvar: variance of the Gaussian in x and y
% bkgnd: background
% A: amplitude
% 
% Note: if you're wondering how well a 2D Gaussian approximates the
% point-spread function of a single fluorophore, look at the Wikipedia
% article on Airy disks (http://en.wikipedia.org/wiki/Airy_disc#Approximation_using_a_Gaussian_profile).
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function [Xcen, Ycen, Xvar, Yvar, bkgnd, A] = Fit2DGaussToSpot(spotimg,varargin)

debug = 0;
symGauss = 0;

if ~isempty(varargin)
    for k = 1:2:length(varargin)
        if strcmpi(varargin{k},'debug')
            debug = varargin{k+1};
        elseif strcmpi(varargin{k},'symGauss')
            symGauss = varargin{k+1};
        end
    end
end

% Start with some intelligent guesses for initial parameters:
A_init = max(spotimg(:)); %Guess that the amplitude is the intensity of the brightest pixel
bkgnd_init = min(spotimg(:)); %Guess that the background is the minimum pixel intensity
Xcen_init = size(spotimg,2)/2; %Guess that the spot is roughly in the center.
Ycen_init = size(spotimg,1)/2; 
    % Alternatively, you could use the location of the brightest pixel to
    % set Xcen_init, Ycen_init.
Xvar_init = 1/(size(spotimg,2)/4); % Assume the user didn't give you a huge ROI for a tiny spot ... 
Yvar_init = 1/(size(spotimg,1)/4);

% Find parameters that minimize the difference between a 2D Gaussian and
% the actual image.  This "minimize the difference" problem is encapsulated
% in the Gauss2DCost function.
if symGauss
    fitparams = fminsearch(@(params)Gauss2DCostSym(params,spotimg),...
        [A_init,bkgnd_init,Xcen_init,Ycen_init,Xvar_init]);
    A = fitparams(1);
    bkgnd = fitparams(2);
    Xcen = fitparams(3);
    Ycen = fitparams(4);
    Xvar = fitparams(5);
    Yvar = Xvar;
else
    fitparams = fminsearch(@(params)Gauss2DCost(params,spotimg),...
        [A_init,bkgnd_init,Xcen_init,Ycen_init,Xvar_init,Yvar_init]);
    A = fitparams(1);
    bkgnd = fitparams(2);
    Xcen = fitparams(3);
    Ycen = fitparams(4);
    Xvar = fitparams(5);
    Yvar = fitparams(6);
end

if debug
    % Debugging: plot a surface map of the spot versus the fit:
    figure('Position',[200,0,900,700])
    subplot(2,2,1)
    surf(spotimg)
    colormap gray
    title('Original image','Fontsize',14)
    zlim([0 1])
    
    subplot(2,2,2)
    surf(PlotGauss2D(size(spotimg),fitparams))
    colormap jet
    title('Best-fit Gaussian','Fontsize',14)
    zlim([0 1])
    
    subplot(2,2,3)
    surf(spotimg)
    hold on
    mesh(PlotGauss2D(size(spotimg),fitparams))
    colormap pink
    title('Overlay','Fontsize',14)
    zlim([0 1])
    
    subplot(2,2,4)
    surf(spotimg-PlotGauss2D(size(spotimg),fitparams))
    colormap hot
    title('Difference','Fontsize',14)
    zlim([0 1])
    
    pause
    close
end
