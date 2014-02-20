% function [Xcen, Ycen, Xvar, Yvar, bkgnd, A] = Fit2DGaussToSpot(spotimg)
%
% Given a region-of-interest (i.e. part of a frame from a FRET movie) that
% nominally contains a single fluorescent spot, find the parameters of the
% 2D Gaussian that best fits that spot.
%
% Input:
% spotimg: image of a single spot
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

function [Xcen, Ycen, Xvar, Yvar, bkgnd, A] = Fit2DGaussToSpot(spotimg)

% Start with some intelligent guesses for initial parameters:
A_init = max(spotimg(:)); %Guess that the amplitude is the intensity of the brightest pixel
bkgnd_init = min(spotimg(:)); %Guess that the background is the minimum pixel intensity
Xcen_init = size(spotimg,2); %Guess that the spot is roughly in the center.
Ycen_init = size(spotimg,1); 
    % Alternatively, you could use the location of the brightest pixel to
    % set Xcen_init, Ycen_init.
Xvar_init = size(spotimg,2)/2; % Assume the user didn't give you a huge ROI for a tiny spot ... 
Yvar_init = size(spotimg,1)/2;

% Find parameters that minimize the difference between a 2D Gaussian and
% the actual image.  This "minimize the difference" problem is encapsulated
% in the Gauss2DCost function.
fitparams = fminsearch(@(params)Gauss2DCost(params,spotimg),...
    [A_init,bkgnd_init,Xcen_init,Ycen_init,Xvar_init,Yvar_init]);

A = fitparams(1);
bkgnd = fitparams(2);
Xcen = fitparams(3);
Ycen = fitparams(4);
Xvar = fitparams(5);
Yvar = fitparams(6);

    % Debugging: plot a surface map of the spot versus the fit:
    %figure
    %imagesc(spotimg)
    %colormap gray
    %hold on
    %plot([outparams(5) outparams(5)],[outparams(4) outparams(4)],'xr')
    
    figure
    surf(spotimg)
    colormap gray
    figure
    mesh(PlotGauss2D(size(spotimg),outparams))
    colormap pink
    pause
    close
