% function whichPolyVersion = CheckPoly()
%
% Uses a set of points with a known transformation to check whether the
% Matlab polynomial transformation calculator in FRETmap or FRETmapR2017
% will work, or whether Traces should use its built-in function. 
%
% Outputs are 'preR2017a', 'R2017a', or 'Traces'.
%
% Called by TracesSetup.
%
% By extension this will also tell Traces which version of affine transform
% code to use.
%
% The MIT License (MIT)
% 
% Copyright (c) 2018 Stephanie Johnson, University of California, San Francisco
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

function whichPolyVersion = CheckPoly()

% Pairs of points to use to check transform calculators:
load('PolyPoints.mat');

MatlabPoly = FRETmap(StartData,EndData,'Green','MatlabPoly',4,4);
MatlabPolyR2017a = FRETmapR2017a(StartData,EndData,'Green','MatlabPoly',4);
Poly = FRETmap(StartData,EndData,'Green','Poly',4,4);

% Traces internal polynomial transformation should work no mattter what, so
% use the residuals on that transform to test whether the different Matlab
% poly versions work:
tolerance = round(Poly.ResidualsFwd,1); % Round the known good residuals to the first decimal point

if tolerance > 0.2 % then something is really wrong!
    disp('CheckPoly: Residuals too high even with Traces built-in polynomial transform calculator;')
    disp('something is wrong!')
    keyboard
end

if MatlabPoly.ResidualsFwd < tolerance
    whichPolyVersion = 'preR2017a';
elseif MatlabPolyR2017a.ResidualsFwd < tolerance
    whichPolyVersion = 'R2017a';
else
    whichPolyVersion = 'Traces';
end