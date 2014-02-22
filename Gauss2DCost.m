% function cost = Gauss2DCost(params,data)
%
% Cost function for fitting a 2D Gaussian to an image.  The cost is the sum
% of the differences between the computed Gaussian, given the parameters in
% params, and the actual image, in data. Called by Fit2DGaussToSpot.m
%
% Params must be a vector of: [x0,y0,xvar,yvar,B,A], where
% Gauss2D = A*exp(-xvar*(x-x0)^2-yvar*(y-y0)^2)+B
%
% Thanks to David Wu for showing me how to write this (and the original
% versions of the other associated functions, Fit2DGaussToSpot and
% PlotGauss2D!)
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function cost = Gauss2DCost(params,data)

    A = params(6);
    B = params(5);
    x0 = params(1);
    y0 = params(2);
    xvar = params(3);
    yvar = params(4);
    
    difference = zeros(size(data));
    
    for i=1:size(data,1)
        for j=1:size(data,2);
            difference(i,j) = (data(i,j)-(A*exp(-xvar*(i-x0)^2-yvar*(j-y0)^2) + B))^2;
        end
    end
    
    cost = sum(sum(difference));

end