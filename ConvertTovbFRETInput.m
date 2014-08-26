% function ConvertTovbFRETInput()
%
% Converts all the Spot*_*.mat files in a directory to a suitable input
% format for vbFRET.
%
% This currently converts the ``unsmoothed'' intensities and FRET values,
% which will have gamma and alpha applied (but are, obviously, unsmoothed).
%
% Steph 8/2014
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

function ConvertTovbFRETInput()

ToAnalyze = uigetdir('','Directory with Spot files to convert?:');

AllSpots = dir(fullfile(ToAnalyze,'Spot*_*.mat'));

FRET = cell(1,length(AllSpots));
data = cell(1,length(AllSpots));
labels = cell(1,length(AllSpots));

for k = 1:length(AllSpots)
    labels{k} = AllSpots(k).name;
    tempspot = load(fullfile(ToAnalyze,AllSpots(k).name));
    % vbFRET wants the first column to be donor
    data{k}(:,1) = tempspot.unsmoothedGrI;
    data{k}(:,2) = tempspot.unsmoothedRedI;
    FRET{k} = tempspot.unsmoothedFRET';
    clear tempspot
end

mkdir(fullfile(ToAnalyze,'vbFRETformat'))
save(fullfile(ToAnalyze,'vbFRETformat','AllSpots.mat'),'FRET','data','labels')

