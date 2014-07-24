% function val = GetInfoFromMetaData(dirname,paramname)
%
% Paramname can be: 
% imgsize ->output will be a vector of xpxls,ypxls
% fps -> output will be frames per second
%
% Update 7/2014: To (hopefully) make loading large data sets faster, this
% function now saves the result to disk so that the next time you have to
% get this information out of the metadata file, it'll be faster.
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

function val = GetInfoFromMetaData(dirname,paramname)

if ~exist(fullfile(dirname,'metadata.mat'),'file')
    fid = fopen(fullfile(dirname,'metadata.txt'));
    alllines = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    %if strcmpi(paramname,'imgsize')
        % The way our metadata files are written, the size of the image
        % will be in different places in the file depending on whether it
        % was collected with Mult-D acquisition or Snap:
        if strcmpi(alllines{1}(19),'"ROI": [')
            tempx = alllines{1}(22);
            tempy = alllines{1}(23);
        else
            tempx = alllines{1}(12);
            tempy = alllines{1}(22);
        end
        tempx2 = tempx{1};
        tempy2 = tempy{1};
        [~,~,~,xmatch] = regexpi(tempx2,'\d\d\d');
        [~,~,~,ymatch] = regexpi(tempy2,'\d\d\d');
        %val(1) = str2double(xmatch);
        %val(2) = str2double(ymatch);
        xsize = str2double(xmatch);
        ysize = str2double(ymatch);
        clear tempx tempy tempx2 tempxy xmatch ymatch fid
        %%%%TODO: the above part needs some error handling ... 
    %elseif strcmpi(paramname,'fps')
        temp = regexpi(alllines{1}(4),'Interval_ms');
        if ~isempty(temp{1})
            thisline = char(alllines{1}(4));
            temp2 = regexpi(thisline,'\d');
            valstr = '';
            for i = 1:length(temp2)
                valstr = strcat(valstr,thisline(temp2(i)));
            end
            % val = str2double(valstr);
            fps = str2double(valstr);
        else
            disp('GetInfoFromMetaData: Unexpected metadata format.')
            val = -1;
            return
        end
    save(fullfile(dirname,'metadata.mat'),'xsize','ysize','fps')
    if strcmpi(paramname,'imgsize')
        val(1) = xsize;
        val(2) = ysize;
    elseif strcmpi(paramname,'fps')
        val = fps;
    else
        disp('GetInfoFromMetaData: Invalid paramname.')
        val = -1;
        return
    end
else
    info = load(fullfile(dirname,'metadata.mat'));
    if strcmpi(paramname,'imgsize')
        val(1) = info.xsize;
        val(2) = info.ysize;
    elseif strcmpi(paramname,'fps')
        val = info.fps;
    else
        disp('GetInfoFromMetaData: Invalid paramname.')
        val = -1;
        return
    end
end