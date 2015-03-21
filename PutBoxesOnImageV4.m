% function PutBoxesOnImageV4(img,spots,boxdim,boxcolor)
%
% Same as PutBoxesOnImageV3, but puts circles on the image instead of boxes.
% Boxdim is the DIAMETER of the circle.
%
% Given a matrix of pixel intensities (ie. an image) and spot centers (in terms
% of the matrix coordinates--that is, (row,col) for spot centers), and a circle 
% dimension (diameter), put circles on top of the image 
% that bound the spots, and an x where it thinks the center is.  Note that 
% spots doesn't have to contain integers.
%
% Combines plot and imshow, so the centers can show up as
% non-integer pixels (V2 rounds to the nearest pixel).
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

function PutBoxesOnImageV4(img,spots,boxdim,NoX,boxcolor)

% Decided not to force the image to be scaled
% if min(img(:))~=0 || max(img(:))~= 1
%     img = mat2gray(img);
% end

if isempty(spots)
    return
end

if ~exist('NoX','var') NoX = 0; end
if ~exist('boxcolor','var') boxcolor = 'g'; end

% Make sure the spots are listed as one spot per column:
if size(spots,1)~=2
    spots = transpose(spots);
end

figure
imshow(img)
hold on

t = 0:pi/100:2*pi;

for j = 1:size(spots,2)
    plot(spots(2,j)+boxdim/2*cos(t),spots(1,j)+boxdim/2*sin(t),strcat('-',boxcolor))
end

% Lastly add a red dot in the center:
% Update 9/2013: a red dot instead of an x
if ~NoX
    % plot(spots(:,2),spots(:,1),'xr')
    plot(spots(2,:),spots(1,:),'.r','MarkerSize',3)
end
hold off

% if ~isempty(varargin)
%     set(varargin{1},'CData',imgwboxes)
% else
%     figure
%     varargout{1} = imshow(imgwboxes);
% end