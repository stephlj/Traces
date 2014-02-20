% function PlotGauss2D(plotsize,params)
%
% Creates a discretized 2D Gaussian for checking the fit performed by
% Fit2DGaussToSpot.
% 
% Inputs:
% plotsize: 2 element vector; extent of the Gaussian (nominally the size of
%   the image that this Gaussian will be plotted over)
% params: the six-element output of the call to fminsearch in
%   Fit2DGaussToSpot. Params must be a vector of: [A,B,x0,y0,xvar,yvar], 
%   where Gauss2D = A*exp(-xvar*(x-x0)^2-yvar*(y-y0)^2)+C.
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function Discrete2DGauss = PlotGauss2D(plotsize,params)

    Discrete2DGauss = zeros(plotsize);
    
    for i=1:plotsize(1)
        for j=1:plotsize(2)
            Discrete2DGauss(i,j)=params(1)*exp(-params(5)*(i-params(3))^2-params(6)*(j-params(4))^2) + params(2);
        end
    end
end