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
% Mode determines whether this function is used in conjunction with
% lsqnonlin (requires just the differences between the data and the equation
% to fit; pass mode = 'diffonly') or fminsearch (requires sum of squares 
% of differences; pass mode = 'sumsquares').
%
% Thanks to David Wu for showing me how to write this (and the original
% versions of the other associated functions, Fit2DGaussToSpot and
% PlotGauss2D!)
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

function cost = Gauss2DCost(params,data,mode,varargin)

    if ~strcmpi(mode,'sumsquares') && ~strcmpi(mode,'diffonly')
        disp('Gauss2DCost: Mode not recognized, should be either sumsquares or diffonly.')
        cost = [];
        return
    end

    if length(params)==6
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
        cost = [];
        return
    end
    
    difference = zeros(size(data));
    
    for i=1:size(data,1)
        for j=1:size(data,2);
            if strcmpi(mode,'sumsquares')
                difference(i,j) = (data(i,j)-(A*exp(-xvar*(i-x0)^2-yvar*(j-y0)^2) + B))^2;
            else
                difference(i,j) = data(i,j) - (A*exp(-xvar*(i-x0)^2-yvar*(j-y0)^2) + B);
            end
        end
    end
    
    if strcmpi(mode,'sumsquares')
        cost = sum(sum(difference));
    else
        cost = difference;
    end

end