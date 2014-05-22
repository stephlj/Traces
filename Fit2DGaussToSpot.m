% function [Xcen, Ycen, Xvar, Yvar, bkgnd, A] = Fit2DGaussToSpot(spotimg,varargin)
%
% Given a region-of-interest (i.e. part of a frame from a FRET movie) that
% nominally contains a single fluorescent spot, find the parameters of the
% 2D Gaussian that best fits that spot.
%
% Inputs:
% spotimg: image of a single spot
% mode: 'full','vars',background' (case insensitive)
%   Use mode = 'full' to fit all 6 variables (5, for a symmetric Gauss)
%   Use mode = 'vars' to fit only the variances and the background (and
%       amplitude), but not the center position
%   Use mode = 'background' to fit only the background (and amplitude), but
%       not the variances or center position.
%   Default is 'full'. If mode is 'vars' or 'background', the optional
%       'StartParams' input must also be passed, and the values for the
%       parameters not to be fit will be used from StartParams.
%
% Optional inputs: enter these as pairs ('<paraname>',<value>)
% 'Debug', 1: display a set of images of the fit
% 'symGauss', 1: force the variances in x and y to be the same
% 'StartParams',[Xcen, Ycen, Xvar, Yvar, bkgnd, A]: start parameters to use
%       for the fit
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

function [Xcen, Ycen, Xvar, Yvar, bkgnd, A] = Fit2DGaussToSpot(spotimg,mode,varargin)

% Error handling
if size(spotimg,3)>1
    disp('Fit2DGaussToSpot: Must pass a 2D image only, not a movie.')
    return
end

% Set mode to default of 'full' if only a spotimg was passed
if ~exist('mode','var') mode = 'full'; end

% Handle optional inputs before dealing further with mode
debug = 0;
symGauss = 0;

% Start with some intelligent guesses for initial parameters:
A_init = max(spotimg(:)); %Guess that the amplitude is the intensity of the brightest pixel
bkgnd_init = min(spotimg(:)); %Guess that the background is the minimum pixel intensity
Xcen_init = size(spotimg,2)/2; %Guess that the spot is roughly in the center.
Ycen_init = size(spotimg,1)/2; 
    % Alternatively, you could use the location of the brightest pixel to
    % set Xcen_init, Ycen_init.
Xvar_init = 1/(size(spotimg,2)/4); % Assume the user didn't give you a huge ROI for a tiny spot ... 
Yvar_init = 1/(size(spotimg,1)/4);

if ~isempty(varargin)
    for k = 1:2:length(varargin)
        if strcmpi(varargin{k},'debug')
            debug = varargin{k+1};
        elseif strcmpi(varargin{k},'symGauss')
            symGauss = varargin{k+1};
        elseif strcmpi(varargin{k},'StartParams') && ~isempty(varargin{k+1})
            A_init = varargin{k+1}(6);
            bkgnd_init = varargin{k+1}(5);
            Xcen_init = varargin{k+1}(1);
            Ycen_init = varargin{k+1}(2); 
            Xvar_init = varargin{k+1}(3);
            Yvar_init = varargin{k+1}(4);
        end
    end
end

switch mode
    case {'full','Full','FULL'}
        fixedparams = [];
        if symGauss
            startparams = [Xcen_init,Ycen_init,Xvar_init,bkgnd_init,A_init];
        else
            startparams = [Xcen_init,Ycen_init,Xvar_init,Yvar_init,bkgnd_init,A_init];
        end
    case {'vars','Vars','VARS'}
        Xcen = Xcen_init;
        Ycen = Ycen_init;
        fixedparams = [Xcen,Ycen];
        if symGauss
            startparams = [Xvar_init,bkgnd_init,A_init];
        else
            startparams = [Xvar_init,Yvar_init,bkgnd_init,A_init];
        end
    case {'background','Background','BACKGROUND'}
        Xcen = Xcen_init;
        Ycen = Ycen_init;
        Xvar = Xvar_init;
        startparams = [bkgnd_init,A_init];
        if symGauss
            Yvar = Xvar;
            fixedparams = [Xcen,Ycen,Xvar];
        else
            Yvar = Yvar_init;
            fixedparams = [Xcen,Ycen,Xvar,Yvar];
        end
    otherwise
        disp('Fit2DGaussToSpot: Mode not recognized.')
        return
end

% Find parameters that minimize the difference between a 2D Gaussian and
% the actual image.  This "minimize the difference" problem is encapsulated
% in the Gauss2DCost function.
% Note if you have the Optimization toolbox, it would be better (faster and
% more accurate) to use lsqnonlin rather than fminsearch. However,
% fminsearch is used here to reduce this software suite's dependance on
% Matlab toolboxes.
%%%% TODO: use a try/catch block to use lsqnonlin if the Optimization
%%%% toolbox is installed.

opts = optimset('Display','off'); % Don't display a warning if the fit doesn't converge

if symGauss
    [fitparams,~,exitflag] = fminsearch(@(params)Gauss2DCostSym(params,spotimg,fixedparams),...
        startparams,opts);
else
    [fitparams,~,exitflag] = fminsearch(@(params)Gauss2DCost(params,spotimg,fixedparams),...
        startparams,opts);
end

% Define output parameters:
if exitflag==0
    % If the fit fails, return the default parameters. Sometimes the
    % parameters it gets from a failed fit are pretty whacky. Return a
    % very low amplitude as indication that this is not a good Gaussian.
    A = 0.0001;
    bkgnd = bkgnd_init;
    Xcen = Xcen_init;
    Ycen = Ycen_init;
    Xvar = Xvar_init;
    Yvar = Xvar_init; 
else
    A = fitparams(end);
    bkgnd = fitparams(end-1);
    switch mode
        case {'full','Full','FULL'}
            Xcen = fitparams(1);
            Ycen = fitparams(2);
            Xvar = fitparams(3);
            if symGauss
                Yvar = Xvar;
            else
                Yvar = fitparams(4);
            end
        case {'vars','Vars','VARS'}
            Xvar = fitparams(1);
            if symGauss
                Yvar = Xvar;
            else
                Yvar = fitparams(2);
            end
    end

end

if debug
    clear fitparams
    if symGauss
        fitparams = [Xcen,Ycen,Xvar,bkgnd,A];
    else
        fitparams = [Xcen,Ycen,Xvar,Yvar,bkgnd,A];
    end
    
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
