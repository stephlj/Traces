% function ConvertTovbFRETInput()
%
% Converts all the Spot*_*.mat files in a directory to a suitable input
% format for vbFRET.
%
% This currently converts the ``unsmoothed'' intensities and FRET values,
% which will have gamma and alpha applied (but are, obviously, unsmoothed).
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

