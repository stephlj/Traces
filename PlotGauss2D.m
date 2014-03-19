% function PlotGauss2D(plotsize,params)
%
% Creates a discretized 2D Gaussian for checking the fit performed by
% Fit2DGaussToSpot.
% 
% Inputs:
% plotsize: 2 element vector; extent of the Gaussian (nominally the size of
%   the image that this Gaussian will be plotted over)
% params: the six-element output of the call to fminsearch in
%   Fit2DGaussToSpot. Params must be a vector of: [x0,y0,xvar,yvar,B,A], 
%   where Gauss2D = A*exp(-xvar*(x-x0)^2-yvar*(y-y0)^2)+B.
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function Discrete2DGauss = PlotGauss2D(plotsize,params)

    Discrete2DGauss = zeros(plotsize);
    
    if length(params)==5
        params(6) = params(5);
        params(5) = params(4);
        params(4) = params(3);
    end

    for i=1:plotsize(1)
        for j=1:plotsize(2)
            Discrete2DGauss(i,j)=params(6)*exp(-params(3)*(i-params(1))^2-params(4)*(j-params(2))^2) + params(5);
        end
    end
end