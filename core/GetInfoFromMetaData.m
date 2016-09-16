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
%    fps -> output ("val") will be frame rate at which the data were
%       acquired (frames per second)
%    precision -> output ("val") will be the numeric type of the raw
%       images. THIS IS ONLY AN OPTION IF YOU CREATED YOUR OWN METADATA FILE,
%       and only necessary if you're loading pma's instead of tif's.
%
% Update 7/2014: To (hopefully) make loading large data sets faster, this
% function now saves the result to disk as a .mat file, so that the next 
% time you have to get this information out of the metadata file, it won't 
% have to parse a text file.
%
% IF YOU DON'T USE MICROMANAGER: Here's how to create a metadata file that
% GetInfoFromMetaData can read:
% Create a text file that has the following lines: (copy them exactly as
% here, except change the "150" to the interval, in MILLIseconds, to whatever
% 1/<your frame rate> is, change the last two numbers to be the pixel size of 
% your camera's field of view, and change the precision to whatever your 
% raw images' numeric type is. NOTE you only need that last line if you're
% going to load pma files! The lines below are for 150 ms frame interval,
% a 512x512 camera, and for raw image files saved in pma's as uint8's. 
% See also the sample metadata.txt file that accompanies this software package):
%
%%%%% SAMPLE METADATA.TXT CONTENT %%%%%%
% "Andor-ActualInterval-ms": 150 
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
% fps = 0.15;
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

function val = GetInfoFromMetaData(dirname,paramname)

    % Update 9/2016: Allow loading via either old uManager or new uManager
    % directory structure:    
    if exist(fullfile(dirname,'Pos0','metadata.txt'),'file')
        dirname = fullfile(dirname,'Pos0');
    end

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
        elseif strcmpi(alllines{1}(19),'"ROI": [') %THis is for previous metadata verison
            tempx = alllines{1}(22);
            tempy = alllines{1}(23);
        elseif strcmpi(alllines{1}(20),'"ROI": [') %THis is for new metadata verison
            tempx = alllines{1}(23);
            tempy = alllines{1}(24);
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
        temp = regexpi(alllines{1}(1),'Andor-ActualInterval-ms');
        if ~isempty(temp{1})
            thisline = char(alllines{1}(1));
            temp2 = regexpi(thisline,'\d');
            valstr = '';
            for i = temp2(1):temp2(end)
                    valstr = strcat(valstr,thisline(i));
            end
            fps = 1/(str2double(valstr)*10^-3);
        else
            temp = regexpi(alllines{1}(61),'Andor-ActualInterval-ms'); %For old metadata version
            if ~isempty(temp{1})
                thisline = char(alllines{1}(61));
                temp2 = regexpi(thisline,'\d');
                valstr = '';
                for i = temp2(1):temp2(end)
                    valstr = strcat(valstr,thisline(i));
                end
                fps = 1/(str2double(valstr)*10^-3);
            else
                temp = regexpi(alllines{1}(62),'Andor-ActualInterval-ms'); %For new metadata version
                if ~isempty(temp{1})
                    thisline = char(alllines{1}(62));
                    temp2 = regexpi(thisline,'\d');
                    valstr = '';
                    for i = temp2(1):temp2(end)
                        valstr = strcat(valstr,thisline(i));
                    end
                    fps = 1/(str2double(valstr)*10^-3);
                else
                    disp('GetInfoFromMetaData: Unexpected metadata format.')
                    val = -1;
                    return
                end
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