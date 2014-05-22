% function cost = Gauss2DCostSym(params,data)
%
% Cost function for fitting a 2D Gaussian to an image.  The cost is the sum
% of the differences between the computed Gaussian, given the parameters in
% params, and the actual image, in data. Called by Fit2DGaussToSpot.m
%
% The difference between this and Gauss2DCost is that this insists the
% Gaussian be rotationally symmetric (i.e. have the same xvar as yvar).
%
% Params must be a vector of: [x0,y0,var,B,A], where
% Gauss2D = A*exp(-var*(x-x0)^2-var*(y-y0)^2)+B
% in which case varargin must be empty; 
% OR
% params must be [xvar,B,A], and varargin must be [x0,y0];
% OR
% params must be [B,A] and varargin must be [x0,y0,xvar].
%
% Thanks to David Wu for showing me how to write this (and the original
% versions of the other associated functions, Fit2DGaussToSpot and
% PlotGauss2D!)
%
% Steph 2/2014, updated 5/2014 to allow for only some parameters to be fit
% and others (those in varargin) to remain constant.
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function cost = Gauss2DCostSym(params,data,varargin)

    if length(params)==5
        A = params(5);
        B = params(4);
        x0 = params(1);
        y0 = params(2);
        var = params(3);
    elseif length(params)==3 && length(varargin{1})==2
        A = params(3);
        B = params(2);
        x0 = varargin{1}(1);
        y0 = varargin{1}(2);
        var = params(1);
    elseif length(params)==2 && length(varargin{1})==3
        A = params(2);
        B = params(1);
        x0 = varargin{1}(1);
        y0 = varargin{1}(2);
        var = varargin{1}(3);
    else
        disp('Gauss2DCostSym: Input parameters unusable.')
        return
    end
    
    difference = zeros(size(data));
    
    for i=1:size(data,1)
        for j=1:size(data,2);
            difference(i,j) = (data(i,j)-(A*exp(-var*(i-x0)^2-var*(j-y0)^2) + B))^2;
        end
    end
    
    cost = sum(sum(difference));

end