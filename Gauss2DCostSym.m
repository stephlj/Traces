% function cost = Gauss2DCostSym(params,data)
%
% Cost function for fitting a 2D Gaussian to an image.  The cost is the sum
% of the differences between the computed Gaussian, given the parameters in
% params, and the actual image, in data. Called by Fit2DGaussToSpot.m
%
% The difference between this and Gauss2DCost is that this insists the
% Gaussian be rotationally symmetric (i.e. have the same xvar as yvar).
%
% Params must be a vector of: [A,B,x0,y0,var], where
% Gauss2D = A*exp(-var*(x-x0)^2-var*(y-y0)^2)+C
%
% Thanks to David Wu for showing me how to write this (and the original
% versions of the other associated functions, Fit2DGaussToSpot and
% PlotGauss2D!)
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function cost = Gauss2DCostSym(params,data)

    A = params(1);
    B = params(2);
    x0 = params(3);
    y0 = params(4);
    var = params(5);
    
    difference = zeros(size(data));
    
    for i=1:size(data,1)
        for j=1:size(data,2);
            difference(i,j) = (data(i,j)-(A*exp(-var*(i-x0)^2-var*(j-y0)^2) + B))^2;
        end
    end
    
    cost = sum(sum(difference));

end