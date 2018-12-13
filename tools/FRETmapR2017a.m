% FRETmapR2017a class
%
% Mostly the same as FRETmap class, but compatible with Matlab's polynomial 
% transform functions in R2017a and later versions. Around version R2017a, 
% Matlab's built-in code switched to normalizing StartData and EndData around 
% (0,0), and then computing the transformation for these normalized points. 
% All of the functions that manipulate the resulting tform object know about 
% this normalization, but it means that the tform is not compatible with the
% rest of Traces in this normalized form.
%
% So the biggest difference for the user is that FRETmapR2017a will NOT
% return the A and Ainv matrices to the user.
%
% Object containing information about an smFRET channel mapping. 
% Uses Matlab's handle class rather than value class.
%
% Properties of this class:
% Kind: MatlabAffine/MatlabPoly
% ResidualsFwd,ResidualsInv: Sum of the squared residuals of the fitted transformation
% TotalDegree: total degree for polynomial transformations. Note that
%   Matlab defines the degree of the transformation in terms of total
%   degree, always. (Not max.) This can be 2, 3, or 4.
% StartData: in the transformation y=Ax, this is x. Forward transformation
%       maps from StartData to EndData
% EndData
% StartChannel: 'Green' or 'Red'. Never used in any consequential way by 
%       the object (the exception is the HistResiduals method); just for the
%       user's benefit so they know if they defined forward as green->red or vice
%       versa
%
% Methods of this class:
% MapObj = FRETmap(StartData,EndData,StartChannel,Kind,TotalDegree)
%       Constructor method
% newspots = MapObj.FRETmapFwd(oldspots)
%       Performs forward transformation (StartChannel->EndChannel)
% newspots = MapObj.FRETmapInv(oldspots)
%       Performs inverse transformation (EndChannel->StartChannel)
% tform = MapObj.ReturnMatlabTform(direction)
%       Returns the MapObj information as a Matlab tform.
%       Direction is a string and is either 'fwd' or 'inv'
% MapObj.HistResiduals(direction)
%       Histograms the residuals from the transformation. Direction is a
%       string and is either 'fwd' or 'inv'
% MapObj.TformResiduals(Data1,Data2,direction)
%       Same as HistResiduals but done for new data (not the
%       data that generated the object's transformation). Note that this
%       method always transforms Data1 in the direction
%       indicated by the string direction (either 'fwd' or 'inv'), and histograms
%       the differences from Data2.
% MapObj.PlotTform(Data1,Data2)
%       Static method. Plots the transformation in various ways. 
%       Calls to this function could look like:
%           MapObj.PlotTform(MapObj.StartData,MapObj.EndData): plots the
%               untransformed StartData vs. untransformed EndData of this
%               object
%           MapObj.PlotTform(MapObj.StartData,MapObj.FRETmapFwd(MapObj.StartData))
%               This would plot the results of performing the forward
%               transformation on the StartData of this object versus the
%               StartData. The first input could instead be EndData, which
%               would check the fidelity of the transformation.
%           MapObj.PlotTform(newData,MapObj.FRETmapFwd(newData))
%               This would plot a comparison between a user-inputted data
%               set and the transform of that data set.
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

classdef FRETmapR2017a < handle
    properties (SetAccess=private)
        Kind % MatlabAffine, MatlabPoly
        tform %Matlab's tform object, for NORMALIZED points
        tformInv
        ResidualsFwd % of the fitted forward transformation
        ResidualsInv % of the fitted inverse transformation
        TotalDegree % Total degree for the polynomial transformation (Matlab is always TOTAL degree)
        StartData % The data used to generate the tform
        EndData % The data used to generate the tform
        StartChannel % 'green' or 'red', whichever corresponds to StartData.
    end
    methods
        % Constructor method
        function self = FRETmapR2017a(StartData,EndData,StartChannel,Kind,TotalDegree)
            % Set up the data fields:
                % Make sure StartData and EndData are 2-by-numspots matrices
                    if size(StartData,1)~=2
                        StartData = transpose(StartData);
                    end
                    if size(EndData,1)~=2
                        EndData = transpose(EndData);
                    end
                    if size(StartData,2)~=size(EndData,2)
                        disp('FRETmap class: StartData and EndData must contain the same numbers of points.')
                        self = -1;
                        return
                    end
                self.StartData = StartData;
                self.EndData = EndData;
                self.StartChannel = StartChannel;
            % Error handling (some additional error handling below)
                % Check that Kind is supported
                if strcmpi(Kind,'MatlabAffine') || strcmpi(Kind,'MatlabPoly')
                    self.Kind = Kind;
                elseif strcmpi(Kind,'Affine') % This class forces the use of Matlab built-in functions only!
                    self.Kind = 'MatlabAffine';
                elseif strcmpi(Kind,'Poly')
                    self.Kind = 'MatlabPoly';
                else
                    disp(strcat('Kind:',Kind,' not supported for class FRETmapR2017a.'))
                    self.Kind = -1;
                    return
                end
            
            % If Kind is MatlabAffine:
            if strcmpi(self.Kind,'MatlabAffine')
                % Check that there are at least three points:
                if size(StartData,2)<3
                    disp('FRETmapR2017a class: Minimum three pairs of points required for affine transformation')
                    self.tform = -1;
                    self.tformInv = -1;
                    return
                end
                % Calculate the transformation:
                % No longer need try/catch here because CheckPoly already
                % found out this will work
                self.tform = fitgeotrans(StartData',EndData','Affine');
                % This is kind of stupid because affine is easily
                % invertible mathematically but whatever:
                self.tformInv = fitgeotrans(EndData',StartData','Affine');
                
            % Kind is MatlabPoly  
            else % More error handling
                % Check that if Kind is not affine, TotalDegree has been given
                if ~exist('TotalDegree','var')
                    self.TotalDegree = 4;
                else
                    self.TotalDegree = TotalDegree;
                end
                % Check that enough points were supplied: this is an upper
                % bound
                if size(StartData,2) < self.TotalDegree^2
                    disp('FRETmapR2017a class: Too few spots for chosen polynomial transform.')
                    disp('Reduce degree of desired polynomial or use affine.')
                    self.tform = -1;
                    self.tformInv = -1;
                    return
                end
                % Calculate the transformation
                % No longer need try/catch here because CheckPoly already
                % found out this will work
                % Note that this is technically the tform for
                % StartData->EndData; but because of Matlab weirdness, to
                % actually compute EndData from StartData below, I'll end
                % up calling transformPointsInverse(tformInv,StartData')
                self.tform = fitgeotrans(StartData',EndData','polynomial',self.TotalDegree);
                self.tformInv = fitgeotrans(EndData',StartData','polynomial',self.TotalDegree);
            end
            
            % Lastly, calculate the residuals:
            self.ResidualsFwd = sum(sum((EndData-self.FRETmapFwd(StartData)).^2));
            self.ResidualsInv = sum(sum((StartData-self.FRETmapInv(EndData)).^2));
        end
        
        % Given a spot or set of spots in StartChannel, find it/them in the
        % other channel
        function newspots = FRETmapFwd(self,oldspots)
            flipoutput = 0;
            if size(oldspots,1)~=2
                oldspots = transpose(oldspots);
                flipoutput = 1;
            end
            newspots = transformPointsInverse(self.tformInv, oldspots')';
            if flipoutput
                newspots = transpose(newspots);
            end
        end
        
        % Given a spot or set of spots in EndChannel, find it/them in the
        % StartChannel
        function newspots = FRETmapInv(self,oldspots)
            flipoutput = 0;
            if size(oldspots,1)~=2
                oldspots = transpose(oldspots);
                flipoutput = 1;
            end
            newspots = transformPointsInverse(self.tform, oldspots')';
            if flipoutput
                newspots = transpose(newspots);
            end
        end
      
        function MatlabTform = ReturnMatlabTform(self,direction)
            if strcmpi(direction, 'fwd')
                % empirically:
                MatlabTform = self.tformInv;
            elseif strcmpi(direction,'inv')
                MatlabTform = self.tform;
            else
                disp(strcat('Direction: ',direction,' not defined for FRETmapR22017a class method ReturnMatlabTform.'))
                return
            end
        end
        
        % Plot the residuals from the fit that produced the mapping
        % contained in this object
        function HistResiduals(self,direction)
            figure
            if strcmpi(direction,'fwd')
                errs = FindSpotDists(self.EndData,self.FRETmapFwd(self.StartData));
            elseif strcmpi(direction,'inv')
                errs = FindSpotDists(self.StartData,self.FRETmapInv(self.EndData));
            else 
                disp(strcat('Direction: ',direction,' not defined for FRETmap class method HistResiduals.'))
                return
            end
            hist(min(errs,[],2),[0:0.1:10])
            hold on
            plot([mean(min(errs,[],2)) mean(min(errs,[],2))], [0 size(errs,1)/4],'--k');
            hold off
            ylabel('Counts','Fontsize',12)
            if strcmpi(direction,'fwd')
                xlabel('Distance between mapped end-channel spot and real end-channel spot','Fontsize',12)
            else
                xlabel('Distance between mapped start-channel spot and real start-channel spot','Fontsize',12)
            end
        end
        % Histogram the differences between a transformed set of points and
        % their real positions
        function TformResiduals(self,Data1,Data2,direction)
            figure
            if strcmpi(direction,'fwd')
                newspots = self.FRETmapFwd(Data1);
            elseif strcmpi(direction,'inv')
                newspots = self.FRETmapInv(Data1);
            else 
                disp(strcat('Direction: ',direction,' not defined for FRETmap class method TformResiduals.'))
                return
            end
            errs = FindSpotDists(Data2,newspots);
            hist(min(errs,[],2),[0:0.1:10])
            hold on
            plot([mean(min(errs,[],2)) mean(min(errs,[],2))], [0 size(errs,1)/4],'--k');
            hold off
            ylabel('Counts','Fontsize',12)
        end
    end
    
    methods(Static)
        % Plot, as a scatter plot, a set of data points and their
        % transformation
        function PlotTform(Data1,Data2)
            % Update 10/2014: Making this show up on both macs and pc's
            % better
            %figure('Position',[200 200 325 625])
            % Update 1/2015: This doesn't work well. It's actually the
            % second element of the Position vector that's the problem.
            if ismac
                figure('Position',[200 200 325 625])
            else
                figure('Position',[200 45 325 625])
            end
%             figure
%             DefaultFigPos = get(gcf,'Position');
%             set(gcf,'Position',[DefaultFigPos(1)-0.5*DefaultFigPos(1), DefaultFigPos(2),...
%                 DefaultFigPos(3)*0.575, DefaultFigPos(4)*1.5])
%             clear DefaultFigPos
            plot(Data1(2,:),0-Data1(1,:),'xg')
            hold on
            plot(Data2(2,:),0-Data2(1,:),'xr')
            hold off
            % ylim([-512 0])
            % xlim([0 256])
        end
    end
end