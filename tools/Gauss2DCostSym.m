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

function cost = Gauss2DCostSym(params,data,mode,varargin)

    if ~strcmpi(mode,'sumsquares') && ~strcmpi(mode,'diffonly')
        disp('Gauss2DCost: Mode not recognized, should be either sumsquares or diffonly.')
        cost = [];
        return
    end

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
        cost = [];
        return
    end
    
    difference = zeros(size(data));
    
    for i=1:size(data,1)
        for j=1:size(data,2);
            if strcmpi(mode,'sumsquares')
                difference(i,j) = (data(i,j)-(A*exp(-var*(i-x0)^2-var*(j-y0)^2) + B))^2;
            else
                difference(i,j) = data(i,j)-(A*exp(-var*(i-x0)^2-var*(j-y0)^2) + B);
            end
        end
    end
    
    if strcmpi(mode,'sumsquares')
        cost = sum(sum(difference));
    else
        cost = difference;
    end

end