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
% Copyright (C) 2014 Stephanie Johnson, University of California, San Francisco
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% A copy of the GNU General Public License can be found in the LICENSE.txt 
% file that accompanies this software; it can also be found at 
% <http://www.gnu.org/licenses/>.

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