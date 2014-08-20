% function I = CalcSpotIntensity(kind,img,spotcen,spotvar,params,bkgnd)
%
% Calculates the intensity in a single spot.  Called by CalcIntensities.
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