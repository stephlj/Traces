% function I = CalcSpotIntensity(kind,img,spotcen,spotvar,params,bkgnd)
%
% Calculates the intensity in a single spot.  Called by CalcIntensities.
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

function I = CalcSpotIntensity(kind,img,spotcen,spotvar,params,bkgnd)

    if ~exist('bkgnd','var')
        bkgnd = 0;
    end
    if ~isempty(params.FixSpotVar)
        spotvar = params.FixSpotVar;
    end

    if strcmpi(kind,'Gauss')
        spotsize = size(img(:,:,1));
        mask = PlotGauss2D(spotsize,...
           [spotcen(1), spotcen(2), spotvar(1), spotvar(2), bkgnd, params.GaussWeightAmp]);
    else
        SpotDiam = 5;
        mask = zeros(SpotDiam,SpotDiam,class(img));
        mask(2:4,:) = ones(3,SpotDiam,class(img));
        mask(:,2:4) = ones(SpotDiam,3,class(img));
    end

    % Instead of using a for-loop here, reshape so that we can take
    % advantage of Matlab's fast matrix multiplication:
    newmask = reshape(mask,1,size(mask,1)*size(mask,2));
    data = reshape(img,size(img,1)*size(img,2),size(img,3));
    I = newmask*data;

end