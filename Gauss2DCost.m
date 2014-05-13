% function cost = Gauss2DCost(params,data,varargin)
%
% Cost function for fitting a 2D Gaussian to an image.  The cost is the sum
% of the differences between the computed Gaussian, given the parameters in
% params (and optionally in varargin), and the actual image, in data. 
% Called by Fit2DGaussToSpot.m
%
% Params must be a vector of: [x0,y0,xvar,yvar,B,A], where
% Gauss2D = A*exp(-xvar*(x-x0)^2-yvar*(y-y0)^2)+B
% in which case varargin must be empty; 
% OR
% params must be [xvar,yvar,B,A], and varargin must be [x0,y0];
% OR
% params must be [B,A] and varargin must be [x0,y0,xvar,yvar].
%
% Thanks to David Wu for showing me how to write this (and the original
% versions of the other associated functions, Fit2DGaussToSpot and
% PlotGauss2D!)
%
% Steph 2/2014, updated 5/2014 to allow for only some parameters to be fit
% and others (those in varargin) to remain constant.
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function cost = Gauss2DCost(params,data,varargin)

    if length(params)==6 && isempty(varargin)
        A = params(6);
        B = params(5);
        x0 = params(1);
        y0 = params(2);
        xvar = params(3);
        yvar = params(4);
    elseif length(params)==4 && length(varargin{1})==2
        A = params(4);
        B = params(3);
        x0 = varargin{1}(1);
        y0 = varargin{1}(2);
        xvar = params(1);
        yvar = params(2);
    elseif length(params)==2 && length(varargin{1})==4
        A = params(2);
        B = params(1);
        x0 = varargin{1}(1);
        y0 = varargin{1}(2);
        xvar = varargin{1}(3);
        yvar = varargin{1}(4);
    else
        disp('Gauss2DCost: Input parameters unusable.')
        return
    end
    
    difference = zeros(size(data));
    
    for i=1:size(data,1)
        for j=1:size(data,2);
            difference(i,j) = (data(i,j)-(A*exp(-xvar*(i-x0)^2-yvar*(j-y0)^2) + B))^2;
        end
    end
    
    cost = sum(sum(difference));

end