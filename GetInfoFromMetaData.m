% function val = GetInfoFromMetaData(dirname,paramname)
%
% Gets information from the metadata.txt file that MicroManager saves for
% each acquisition (or alternatively, the metadata.txt file that the user
% created for each acquisition--see below).
% 
% Inputs:
% dirname: path to the directory that contains the metadata.txt file
% paramname: can be
%    imgsize ->output ("val") will be a vector of xpxls,ypxls
%    fps -> output ("val") will be frames per MILLIsecond not frames per
%       second (I know, bad nomenclature ... )
%    precision -> output ("val") will be the numeric type of the raw
%       images. THIS IS ONLY AN OPTION IF YOU CREATED YOUR OWN METADATA FILE,
%       and only necessary if you're loading pma's instead of tif's.
%
% Update 7/2014: To (hopefully) make loading large data sets faster, this
% function now saves the result to disk as a .mat file, so that the next 
% time you have to get this information out of the metadata file, it'll be 
% faster.
%
% IF YOU DON'T USE MICROMANAGER: Here's how to create a metadata file that
% GetInfoFromMetaData can read:
% Create a text file that has the following lines: (copy them exactly as
% here, except change the "150" to be whatever your frame interval is, in 
% milliseconds, change the last two numbers to be the pixel size of 
% your camera's field of view, and change the precision to whatever your 
% raw images' numeric type is. NOTE you only need that last line if you're
% going to load pma files! The lines below are for 150 ms frame interval,
% a 512x512 camera, and for raw image files saved in pma's as uint8's. 
% See also the sample metadata.txt file that accompanies this software package):
%
%%%%% SAMPLE METADATA.TXT CONTENT %%%%%%
% "Interval_ms": 150 
% "ROI": [
%    0,
%    0,
%    512,
%    512
% ]
% "Precision": uint8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the above as "metadata.txt".  Each data-containing directory will 
% need its own metadata file.
%
% Alternatively, you can create a metadata.mat file that has fields fps,
% xsize and ysize, and precision (if you need it), since GetInfoFromMetaData 
% checks first for a metadata.mat file before looking for a .txt file. 
% Create such a metadata.mat file as follows, from the Matlab command line:
%
%%%%% SAMPLE METADATA.MAT CREATION COMMANDS %%%%%%
% fps = 150;
% xsize = 512;
% ysize = 512;
% save('metadata.mat','fps','xsize','ysize')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Again each data-containing directory will need its own metadata file.
%
% Lastly, you could also edit smFRETsetup to have fps, xsize and ysize
% included in the params structure that gets passed around between
% functions, and bypass GetInfoFromMetaData entirely. I would recommend
% against this, though, because having a file the user has to generate for
% each acquisition ensures that they think about what fps should be. If fps
% is one of many parameters in smFRETsetup, users might forget to modify it
% so that it has the correct frame interval rate.
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

if ~exist(fullfile(dirname,'metadata.mat'),'file') && ...
        ~exist(fullfile(dirname,'metadata.txt'),'file')
    disp('GetInfoFromMetaData: File does not exist?')
    val = -1;
    return
end

if ~exist(fullfile(dirname,'metadata.mat'),'file')
    fid = fopen(fullfile(dirname,'metadata.txt'));
    alllines = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    clear fid
    % First extract image size:
        % The way our metadata files are written, the size of the image
        % will be in different places in the file depending on whether it
        % was collected with Mult-D acquisition or Snap:
        % Update 7/2014 to also allow for a simple user-created metadata
        % file:
        if strcmpi(alllines{1}(2),'"ROI": [')
            tempx = alllines{1}(5);
            tempy = alllines{1}(6);
        elseif strcmpi(alllines{1}(19),'"ROI": [')
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
        clear tempx tempy tempx2 tempxy xmatch ymatch
        %%%%TODO: the above part needs some error handling ... 
    %elseif strcmpi(paramname,'fps')
    
    % Now extract the fps
        temp = regexpi(alllines{1}(1),'Interval_ms');
        if ~isempty(temp{1})
            thisline = char(alllines{1}(1));
            temp2 = regexpi(thisline,'\d');
            valstr = '';
            for i = 1:length(temp2)
                valstr = strcat(valstr,thisline(temp2(i)));
            end
            fps = str2double(valstr);
        else
            temp = regexpi(alllines{1}(4),'Interval_ms');
            if ~isempty(temp{1})
                thisline = char(alllines{1}(4));
                temp2 = regexpi(thisline,'\d');
                valstr = '';
                for i = 1:length(temp2)
                    valstr = strcat(valstr,thisline(temp2(i)));
                end
                fps = str2double(valstr);
            else
                disp('GetInfoFromMetaData: Unexpected metadata format.')
                val = -1;
                return
            end
        end
        clear temp temp2
        
    % Lastly, get numeric type, if available:
        try
            temp = regexpi(alllines{1}(8),'Precision');
            if ~isempty(temp{1})
                thisline = char(alllines{1}(8));
                % Doing the following so it doesn't matter if the user put
                % a space after the : or not:
                temp2 = regexpi(thisline(13:end),'[uids]*');
                precision = thisline(13+temp2-1:end);
            else
                disp('GetInfoFromMetaData: Unexpected metadata format for precision.')
                val = -1;
                return
            end
            clear temp
        catch
        end
    
    % Save all this to a .mat file for future use:
    if ~exist('precision','var')
        save(fullfile(dirname,'metadata.mat'),'xsize','ysize','fps')
    else
        save(fullfile(dirname,'metadata.mat'),'xsize','ysize','fps','precision')
    end
    
    % Finally, return the value that the user asked for:
    if strcmpi(paramname,'imgsize')
        val(1) = xsize;
        val(2) = ysize;
    elseif strcmpi(paramname,'fps')
        val = fps;
    elseif strcmpi(paramname,'precision') && exist('precision','var')
        val = precision;
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
    elseif strcmpi(paramname,'precision') && isfield(info,'precision')
        val = info.precision;
    else
        disp('GetInfoFromMetaData: Invalid paramname.')
        val = -1;
        return
    end
end