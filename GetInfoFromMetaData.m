% function val = GetInfoFromMetaData(dirname,paramname)
%
% Paramname can be: 
% imgsize ->output will be a vector of xpxls,ypxls
% fps -> output will be frames per second
%
% Steph 10/2013
% Copyright 2013 Stephanie Johnson, University of California, San Francisco

function val = GetInfoFromMetaData(dirname,paramname)

    fid = fopen(fullfile(dirname,'metadata.txt'));
    alllines = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    if strcmpi(paramname,'imgsize')
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
        val(1) = str2double(xmatch);
        val(2) = str2double(ymatch);
        clear tempx tempy tempx2 tempxy xmatch ymatch fid
        %%%%TODO: the above part needs some error handling ... 
    elseif strcmpi(paramname,'fps')
        temp = regexpi(alllines{1}(4),'Interval_ms');
        if ~isempty(temp{1})
            thisline = char(alllines{1}(4));
            temp2 = regexpi(thisline,'\d');
            valstr = '';
            for i = 1:length(temp2)
                valstr = strcat(valstr,thisline(temp2(i)));
            end
            val = str2double(valstr);
        else
            disp('Unexpected metadata format.')
            val = -1;
            return
        end
    else
        disp('Invalid paramname.')
        val = -1;
        return
    end
end