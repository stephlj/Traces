% function [allimgs,numframes] = LoadPMA(D,numtype,varargin)
%
% Loads frames from a .pma file (the binary filetype used in the Ha lab's 
% IDL code) into a 3-d matrix.  It is assumed that the first four bytes of 
% the file encode the camera chip size in pixels (e.g. 512 by 512). More
% specifically, the first two bytes should encode the x-width, and the 
% second two bytes the y-width. If you created a .pma some other way, 
% you'll need to modify this function accordingly.
%
% NOTE: You will also have to create a "metadata" file that 
% GetInfoFromMetaData.m can read. Type "help GetInfoFromMetaData" in the 
% command line for more information.
% 
% Inputs:
% D: the full path, including filename, to the .pma file (but does not 
%   have to include the .pma extension)
% numtype: Numeric type of the data in the .pma file. This should be 'uint8' 
%   if you're using the Ha lab's acquisition code or the UCSF mrc to pma
%   converter.
% 
% Optional inputs: Pass as a '<name>',<val> pair
% 'FramesToLoad',[start end]: allows the user to specify how many images to
%    load, in the form of [start end] vector.
%
% Outputs:
% allimgs: an image matrix. Images returned are NOT scaled between 0 and 1 
%    (or scaled at all), and are returned as the same integer type in the pma file.
%    Acceptor channel ends up on the left (at least if you created a pma using
%    the Ha lab code or UCSF's mrc to pma converter).
% numframes: total number of frames in this pma (regardless of how many you
%    load)
%
% Note that this does not make use of the FrameLoadMax parameter in 
% smFRETsetup!  You can load as many frames as you want with this function,
% including so many it'll crash Matlab if Matlab doesn't have enough
% memory ...
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

function [allimgs,totframes] = LoadPMA(D,numtype,varargin)

    % Deal with inputs
    if ~strcmpi(D(end-2:end),'pma')
        D = strcat(D,'.pma');
    end
    if ~exist(D,'file')
        disp('LoadPMA: File does not exist? Type "help LoadPMA" for more information on input format.');
        allimgs = -1;
        return
    end
    if ~strcmpi(numtype,'uint8') && ~strcmpi(numtype,'int8') && ... 
            ~strcmpi(numtype,'uint16') && strcmpi(numtype,'int16') && ...
            ~strcmpi(numtype,'double') && strcmpi(numtype,'single') && ...
            ~strcmpi(numtype,'uint32') && ~strcmpi(numtype,'int32') && ... 
            ~strcmpi(numtype,'uint64') && strcmpi(numtype,'int64') && ...
        disp('LoadPMA: Numeric type not recognized.')
        allimgs = -1;
        return
    end
    if ~isempty(varargin) && strcmpi(varargin{1},'FramesToLoad')
        StartStop = sort(varargin{1+1});
        if StartStop(1)<=0
            StartStop(1)=1;
        end
    end
    
    % Open the pma as a binary file:
    pmafile = fopen(D);
    
    if pmafile==-1
        disp('LoadPMA: File unreadable.')
        allimgs = -1;
        return
    end
    
    % The first four BYTES of the .pma are the size of each frame. Note
    % that this has to be read as a 16-bit integer, because you can't get
    % the number 512 out of 8 bits of information in binary. (And I doubt
    % you have a camera that's smaller than 512 pixels on a side, in this
    % day and age ... )
    xpxls = fread(pmafile,1,'uint16');
    ypxls = fread(pmafile,1,'uint16');
    
    % Figure out how many frames are in this pma, based on the file size:
    fileinfo = dir(D);
    if round(fileinfo.bytes/2)==fileinfo.bytes/2
        totbytes = fileinfo.bytes-4;
    else
        totbytes = fileinfo.bytes-5;
    end
    totframes = totbytes/(xpxls*ypxls);
    clear totbytes fileinfo
    
    if round(totframes)~=totframes
        disp('LoadPMA: Cannot get a complete set of frames out of this file. Something is wrong.')
        allimgs = -1;
        return
    end
    if StartStop(2)>totframes
        StartStop(2) = totframes;
    end

    if isempty(varargin) || (StartStop(2)-StartStop(1)+1)>=totframes
        allimgs = zeros(xpxls,ypxls,totframes,numtype);
    else
        allimgs = zeros(xpxls,ypxls,(StartStop(2)-StartStop(1)+1),numtype);
    end

    % Load all the images into a 3d matrix:
    if isempty(varargin)
        for i = 1:totframes
            img = zeros(xpxls,ypxls,numtype);
            img = flipud(transpose(fread(pmafile,size(img),strcat('*',numtype)))); 
                % img is filled COLUMN-WISE, hence the transpose (the pma's
                % are created row-wise); flipud is because they're also
                % filled bottom to top, whereas fread pulls them out top to
                % bottom; the asterisk is so the output is returned in the same
                % format as the input
            allimgs(:,:,i) = img;
            clear img
        end
    else
        incr = 1;
        % Get the pointer to the right place in the file:
        if StartStop(1) > 1
            junk = fread(pmafile,xpxls*ypxls*(StartStop(1)-1),numtype);
            clear junk
        end
        for i = StartStop(1):StartStop(2);
            img = zeros(xpxls,ypxls,numtype);
            img = flipud(transpose(fread(pmafile,size(img),strcat('*',numtype)))); 
            allimgs(:,:,incr) = img;
            clear img
            incr = incr+1;
        end
    end
    
    fclose(pmafile);
end