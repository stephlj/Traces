% function [movRed,movGreen,movRedBkgnd,movGrBkgnd,lastframe] = LoadScaledMovie(PathToMovie,frames,params,varargin)
%
% Loads the scaled and background-subtracted movie frames saved by
% ScaleMovie.m. 
%
% Inputs:
% PathToMovie: Full path to folder with the "ScaledMovieFrames..." files
% frames: [start end] vector of frames to return. 
% params: parameters file saved by smFRETsetup.
% Optional: you must also pass 'bkgnd' as the last input, if you want it
%   to return the third and fourth outputs (see below). Update 8/2014:
%   pass 'rbkgnd' to only load the acceptor background image, 'gbkgnd' to
%   only load the donor bkgnd image, 'bkgnd' to load both.
%
% Outputs:
% movRed, movGreen: x-by-y-by-frames matrix of intensities for the red and
%   green channels
% movRedBkgnd, movGrBkgnd: images used to compute background for intensity
%   calculation (currently only called by CalcIntensitiesV3)
% lastframe: The last frame returned.  This is frames(2) if the user asked
%   for fewer than 100 frames to be returned; else this is start+99.
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

function [movRed,movGreen,movRedBkgnd,movGrBkgnd,lastframe] = LoadScaledMovie(PathToMovie,frames,params,varargin)

    % Input error handling
    frames = sort(frames);
    if frames(1)<=0
        frames(1)=1;
    end
    if frames(1)==frames(2)
        frames(2)=frames(1)+1;
    end
    % Make sure all files and frames exist
    [~,totframes] = LoadRawImgs(PathToMovie,'FramesToLoad',[1 1]);
    if frames(1) > totframes
        disp('LoadScaledMovie: Movie not that long?')
        movRed = -1;
        movGreen = -1;
        lastframe = -1;
        movRedBkgnd = -1;
        movGrBkgnd = -1;
        return
    end
    if frames(2) > totframes
        frames(2) = totframes;
        if frames(1) == frames(2)
            frames(2) = frames(2)+1;
        end
    end
    
    % Make sure the scaling information is available:
    if ~exist(fullfile(PathToMovie,strcat('ScalingInfo.mat')),'file')
        disp('LoadScaledMovie: Scaling info not saved?')
        movRed = -1;
        movGreen = -1;
        lastframe = -1;
        movRedBkgnd = -1;
        movGrBkgnd = -1;
        return
    else
        load(fullfile(PathToMovie,strcat('ScalingInfo.mat')));
    end
    
    % For some reason, if you initialize movRed to be
    % movRed = [];
    % size(movRed,3) = 1 ... 
    movRed = zeros(0,0,0);
    movGreen = zeros(0,0,0);
    movRedBkgnd = zeros(0,0,0);
    movGrBkgnd = zeros(0,0,0);
    
    % First, load and scale images:
    firstframe = frames(1);
    lastframe = firstframe;
    while lastframe<frames(2)
        lastframe = min(firstframe+params.FrameLoadMax-1, frames(2));
        [moviebit,~] = LoadRawImgs(PathToMovie,'FramesToLoad',[firstframe lastframe]);
        moviebit = double(moviebit);
        
        if params.NormImage
            tempMeds(1,1,:) = allMedians(firstframe:lastframe);
            % moviebit = moviebit./repmat(tempMeds,size(moviebit,1),size(moviebit,2),1);
            % For backwards compatibility with older Matlab versions, using
            % the other syntax for repmat:
            moviebit = moviebit./repmat(tempMeds,[size(moviebit,1),size(moviebit,2),1]);
            clear tempMeds
        end
        
        if ~params.ScaleChannelsSeparately
            ScaledMovie = mat2gray(moviebit,[MovieMin MovieMax]); % This also converts it to double precision,
                % but need to explicitly do so earlier in case NormImage is 1
            [movRed,movGreen] = SplitImg(ScaledMovie,params);
        else
            [movRed,movGreen] = SplitImg(moviebit,params);
            movRed = mat2gray(movRed,[MovieMinRed,MovieMaxRed]);
            movGreen = mat2gray(movGreen,[MovieMinGr,MovieMaxGr]);
        end
        clear ScaledMovie
    end
    
    % Then, load and return background images if requested:
    if ~isempty(varargin) %User wants to load background images, which are saved
            % in separate files
        % Figure out how many frames were saved per file:
        allfiles = dir(fullfile(PathToMovie,'BackgroundImgs*.mat'));
        % Update 11/2018: Since I delete background images (but not scaling
        % info) after I'm done analyzing a movie, if I then want to re-load
        % the movie, Traces crashes here because there aren't any
        % background files to load!
        if isempty(allfiles)
            ComputeBackgroundImgs(PathToMovie,params);
            allfiles = dir(fullfile(PathToMovie,'BackgroundImgs*.mat'));
        end
        sample_filename = allfiles(1).name; % This isn't necessarily ScaledMovieFrames1to100!!
        first_num = regexpi(sample_filename,'Imgs\d+to','match');
        second_num = regexpi(sample_filename,'to\d+.mat','match');
        FrameSaveIncr = str2double(second_num{1}(3:end-4))-str2double(first_num{1}(5:end-2))+1;
        clear sample_filename first_num second_num
    
        % Converts global frame number to per-file frame number
        FileFrame = @(frame) rem(frame - 1, FrameSaveIncr) + 1;
        
        firstframe = frames(1);
        lastframe = firstframe;

        while lastframe<frames(2)
            FileNum = firstframe-FileFrame(firstframe) + 1;

            if strcmpi(varargin{1},'rbkgnd')
                temp = load(fullfile(PathToMovie,sprintf('BackgroundImgs%dto%d.mat',...
                        FileNum,FileNum+FrameSaveIncr-1)),'imgRBkgnd');
            elseif strcmpi(varargin{1},'gbkgnd')
                temp = load(fullfile(PathToMovie,sprintf('BackgroundImgs%dto%d.mat',...
                        FileNum,FileNum+FrameSaveIncr-1)),'imgGBkgnd');
            else
                temp = load(fullfile(PathToMovie,sprintf('BackgroundImgs%dto%d.mat',...
                        FileNum,FileNum+FrameSaveIncr-1)));
            end

            lastframe = min(FileNum + FrameSaveIncr-1, frames(2));

            if strcmpi(varargin{1},'rbkgnd')
                movRedBkgnd(:,:,end+1:end+1+lastframe-firstframe) = temp.imgRBkgnd(:,:,...
                    FileFrame(firstframe):FileFrame(lastframe));
            elseif strcmpi(varargin{1},'gbkgnd')
                movGrBkgnd (:,:,end+1:end+1+lastframe-firstframe) = temp.imgGBkgnd(:,:,...
                    FileFrame(firstframe):FileFrame(lastframe));
            else
                movRedBkgnd(:,:,end+1:end+1+lastframe-firstframe) = temp.imgRBkgnd(:,:,...
                    FileFrame(firstframe):FileFrame(lastframe));
                movGrBkgnd (:,:,end+1:end+1+lastframe-firstframe) = temp.imgGBkgnd(:,:,...
                    FileFrame(firstframe):FileFrame(lastframe));
            end
            
            firstframe = lastframe+1;

            clear temp
        end
    end
end
    