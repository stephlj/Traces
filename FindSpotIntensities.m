%function SpotIntensities = FindSpotIntensities(spots)
%
%Takes as its input the cell array of ROIs that's the output of
%FindSpotsV2, and finds the total intensity in each ROI, and returns those
%intensities in a vector of length length(spots).
%
%Steph 6/2013

function SpotIntensities = FindSpotIntensities(spots)

SpotIntensities = zeros(length(spots),1);

for i = 1:length(spots)
    SpotIntensities(i) = sum(sum(spots{i}));
end

%That almost wasn't worth a separate function ... maybe I'll turn this into
%a function that also plots traces or something ... 