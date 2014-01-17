%function PutBoxesOnImageV4(img,spots,boxdim,boxcolor)
%
%Same as PutBoxesOnImageV3, but puts circles on the image instead of boxes.
%Boxdim is the DIAMETER of the circle.
%
%Given a matrix of pixel intensities (ie. an image) and spot centers (in terms
%of the matrix coordinates--that is, (row,col) for spot centers), and a circle 
%dimension (diameter), put circles on top of the image 
%that bound the spots, and an x where it thinks the center is.  Note that 
%spots doesn't have to contain integers.
%
%Combines plot and imshow, so the centers can show up as
%non-integer pixels (V2 rounds to the nearest pixel).
%
%Steph 10/2013
%Copyright 2013 Stephanie Johnson, University of California, San Francisco

function PutBoxesOnImageV4(img,spots,boxdim,NoX,boxcolor)

%Decided not to force the image to be scaled
% if min(img(:))~=0 || max(img(:))~= 1
%     img = mat2gray(img);
% end

if ~exist('NoX','var') NoX = 0; end
if ~exist('boxcolor','var') boxcolor = 'g'; end

%Make sure the spots are listed as one spot per row:
if size(spots,2)~=2
    spots = transpose(spots);
end

figure
imshow(img)
hold on

t = 0:pi/100:2*pi;

for j = 1:length(spots)
    plot(spots(j,2)+boxdim/2*cos(t),spots(j,1)+boxdim/2*sin(t),strcat('-',boxcolor))
end

%Lastly add a red dot in the center:
%Update 9/2013: a red dot instead of an x
if ~NoX
    %plot(spots(:,2),spots(:,1),'xr')
    plot(spots(:,2),spots(:,1),'.r','MarkerSize',3)
end
hold off

% if ~isempty(varargin)
%     set(varargin{1},'CData',imgwboxes)
% else
%     figure
%     varargout{1} = imshow(imgwboxes);
% end