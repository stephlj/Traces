% function mappedspots = CalcSpotTransform(Gspots,Rspots,A,b)
%
% Given a mapping FROM GREEN CHANNEL TO RED, parameterized by matrix A and
% vector b (the outputs of CalcChannelMapping or of fitgeotrans), and a spot 
% or set of spots in the green channel (Gspots) OR in the red channel (Rspots),
% calculates where the spot or spots are in the other channel.  Note one of
% Gspots or Rspots must be an empty vector, so it knows which way to
% calculate the mapping.
%
% Note if your map goes from red channel to green, then enter green spots
% in the Rspots parameter, or vice versa.  The point is, the map
% parameterized by A and b takes the first input and maps to the channel
% represented by the second input, or the second input and maps to the
% channel represented by the first input.
%
% Steph 2/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

function mappedspots = CalcSpotTransform(Gspots,Rspots,A,b)

%Input error handling
if ~isempty(Gspots) && ~isempty(Rspots)
    disp('CalcSpotTransform accepts only one spot or spots to map!')
    mappedspots = -1;
    return
end

if size(A,1)~=2 || size(A,2)~=2
    disp('CalcSpotTransform: A must be a 2x2 matrix.')
    mappedspots = -1;
    return
end

if ~(size(b,1)==2 && size(b,2)==1) && ~(size(b,1)==1 && size(b,2)==2)
    disp('CalcSpotTransform: b must be a 2x1 vector.')
    mappedspots = -1;
    return
end

% Make sure spots are 2xnumspots:
if ~isempty(Gspots) && size(Gspots,1)~=2
    Gspots = transpose(Gspots);
end
if ~isempty(Rspots) && size(Rspots,1)~=2
    Rspots = transpose(Rspots);
end

% Do the channel mapping
% Figure out which direction the user wants to map
if ~isempty(Gspots) %Map from green channel to red
    % Is this a single spot or a set of spots?
    if size(Gspots,2)>1
        b = repmat(b,1,size(Gspots,2));
    end
    mappedspots = A*Gspots+b;
else %Map from red to green
    % Again is this a single spot or set of spots?
    if size(Rspots,2)>1
        b = repmat(b,1,size(Rspots,2));
    end
    mappedspots = inv(A)*(Rspots-b);
end

