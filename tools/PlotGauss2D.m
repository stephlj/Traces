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
% The MIT License (MIT)
% 
% Copyright (c) 2014 Stephanie Johnson, University of California, San Francisco
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

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