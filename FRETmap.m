% FRETmap class
%
% Object containing information about an smFRET channel mapping. 
% Uses Matlab's handle class rather than value class.
%
% Properties of this class:
% Kind: Affine/Poly/MatlabAffine/MatlabPoly
% A: Forward transformation matrix
% Ainv: Reverse transformation matrix. For polynomial transformations this
%       is really the transformation matrix; for affine, this is an artificial
%       construction of [inv(M) b] where the inverse transformation is 
%       x = inv(M)(y-b).
% ResidualsFwd,ResidualsInv: Sum of the squared residuals of the fitted transformation
% MaxDegree: maximum degree to use for polynomial transformations
% TotalDegree: total degree for polynomial transformations. Not required
%       for constructor method unless you explicitly want to set it 
%       (MaxDegree, though, is required).
% StartData: in the transformation y=Ax, this is x. Forward transformation
%       maps from StartData to EndData
% EndData
% StartChannel: 'Green' or 'Red'. Never used in any consequential way by 
%       the object (the exception is the HistResiduals method); just for the
%       user's benefit so they know if they defined forward as green->red or vice
%       versa
%
% Methods of this class:
% MapObj = FRETmap(StartData,EndData,StartChannel,Kind,MaxDegree,TotalDegree)
%       Constructor method
% newspots = MapObj.FRETmapFwd(oldspots)
%       Performs forward transformation
% newspots = MapObj.FRETmapInv(oldspots)
%       Performs inverse transformation
% tform = MapObj.ReturnMatlabTform(direction)
%       Returns the MapObj information as a Matlab tform (either a tform
%       object, if using a version of Matlab with fitgeotrans; or a tform
%       structure for older Matlab versions, which goes with cp2tform).
%       Direction is a string and is either 'fwd' or 'inv'
% MapObj.HistResiduals(direction)
%       Histograms the residuals from the transformation. Direction is a
%       string and is either 'fwd' or 'inv'
% MapObj.TformResiduals(self,Data1,Data2,direction)
%       Same as HistResiduals but done for new data (not the
%       data that generated the object's transformation). Note that this
%       method always transforms Data1 in the direction
%       indicated by the string direction (either 'fwd' or 'inv'), and histograms
%       the differences from Data2.
% TformMatrix = CalcAffineTform(pointsStart,pointsEnd)
%       Static method. Calculates an affine transformation without calling
%       any of Matlab's built-in functions.
% TformInv = CalcAffineInv(TformFwd)
%       Static method. Given an affine transformation matrix, calculates
%       the artificial inverse stored in a FRETmap object as Ainv. (See
%       above for precise definition.)
% newStart = EnumerateMonomials(StartVec,MaxDeg,TotDeg)
%       Static method. For a 2xN matrix StartVec, computes monomials of orders 
%       up to and including MaxDeg and TotDeg for each column. The result
%       will be num_monomials by N.
% TformPoly = CalcPolyTform(StartMonomials,EndPts)
%       Static method. Computes polynomial transformation without calling
%       any of Matlab's built-in functions for doing so. StartMonomials
%       MUST be an output of EnumerateMonomials.
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
% Stephanie 4/2014
% Copyright 2014 Stephanie Johnson, University of California, San Francisco

classdef FRETmap < handle
    properties (SetAccess=private)
        Kind % Affine, poly, MatlabAffine, MatlabPoly
        A % Forward transformation matrix (ie from StartData to EndData)
        Ainv % Reverse transformation matrix
        ResidualsFwd % of the fitted forward transformation
        ResidualsInv % of the fitted inverse transformation
        MaxDegree % Max degree, for polynomial tforms
        TotalDegree % Total degree, for polynomial tforms
        StartData % The data used to generate the tform
        EndData % The data used to generate the tform
        StartChannel % 'green' or 'red', whichever corresponds to StartData.
    end
    methods
        % Constructor method
        function self = FRETmap(StartData,EndData,StartChannel,Kind,MaxDegree,TotalDegree)
            % Set up the data fields:
                % Make sure StartData and EndData are 2-by-numspots matrices
                    if size(StartData,1)~=2
                        StartData = transpose(StartData);
                    end
                    if size(EndData,1)~=2
                        EndData = transpose(EndData);
                    end
                self.StartData = StartData;
                self.EndData = EndData;
                self.StartChannel = StartChannel;
            % Error handling (some additional error handling below)
                % Check that Kind is supported
                if strcmpi(Kind,'affine') || strcmpi(Kind,'poly') || ...
                        strcmpi(Kind,'MatlabAffine') || strcmpi(Kind,'MatlabPoly')
                    self.Kind = Kind;
                else
                    disp(strcat('Kind:',Kind,' not supported for class FRETmap.'))
                    self.Kind = -1;
                    return
                end
            
            % If Kind is Affine or MatlabAffine:
            if strcmpi(Kind,'affine') || strcmpi(Kind,'MatlabAffine')
                % Calculate the transformation:
                % Embed in a linear space
                if strcmpi(Kind,'affine')
                    self.A = self.CalcAffineTform(StartData,EndData);
                else % Use Matlab's built-in function IF AVAILABLE
                    try
                        tform = fitgeotrans(EndData',StartData','Affine');
                        self.A = transpose(inv(tform.T)); % Putting Matlab's output into my format
                        clear tform
                    catch
                        try 
                            tform = cp2tform(EndData',StartData','Affine');
                            self.A = transpose(inv(tform.tdata.T)); % Putting Matlab's output into my format
                            % Ooh cp2tform gives the inverse automatically
                            % too. Oh but that's not consistent with my
                            % notation ... not sure what Tinv is so not
                            % using it.
                            %self.Ainv = transpose(inv(tform.tdata.Tinv));
                            clear tform
                        catch 
                            disp('FRETmap class: No built-in Matlab affine transformation function available; computing manually.')
                            self.Kind = 'Affine';
                            % Just do it manually
                            self.A = self.CalcAffineTform(StartData,EndData);
                        end
                    end
                end
                self.Ainv = self.CalcAffineInv(self.A);
                
            % Kind is Poly or MatlabPoly  
            else % More error handling
                % Check that if Kind is not affine, Max and/or TotalDegree
                % have been given
                if ~exist('MaxDegree','var') && ~exist('TotalDegree','var')
                    self.MaxDegree = 1;
                    self.TotalDegree = 1;
%                 elseif ~exist('MaxDegree','var') % User passed only a TotalDegree
%                     self.TotalDegree = TotalDegree;
%                     self.MaxDegree = TotalDegree; %Max can be as great as Total
                elseif ~exist('TotalDegree','var') %User passed only MaxDegree
                    self.MaxDegree = MaxDegree;
                    self.TotalDegree = 2*MaxDegree; % For every MaxDegree, the total degree (for
                        % a polynomial of 2 variables, like those
                        % implemented here) could be as great as twice
                        % MaxDegree. Assume that maximum.
                else
                    self.MaxDegree = MaxDegree;
                    self.TotalDegree = TotalDegree;
                end
                % Calculate the transformation
                if strcmpi(Kind,'poly')
                    x = self.EnumerateMonomials(StartData,self.MaxDegree,self.TotalDegree);
                    self.A = self.CalcPolyTform(x,EndData);
                    clear x
                    x = self.EnumerateMonomials(EndData,self.MaxDegree,self.TotalDegree);
                    self.Ainv = self.CalcPolyTform(x,StartData);
                    clear x
                else % Use Matlab's built-in functionality IF AVAILABLE
                    try
                        tform = fitgeotrans(EndData',StartData','polynomial',self.TotalDegree);
                        % Degree for Matlab is total degree, always
                        self.A = [tform.A;tform.B]; % Putting Matlab's output into my format
                        % NOTE I order the monomials different, so while my
                        % code and Matlab's generally give the same output,
                        % the order of the coefficients will be different!
                        % The A matrices are not interchangeable unless you
                        % reorder the columns.
                        clear tform
                        tform = fitgeotrans(StartData',EndData','polynomial',self.TotalDegree);
                        self.Ainv = [tform.A;tform.B];
                        clear tform
                    catch
                        try 
                            tform = cp2tform(EndData',StartData','polynomial',self.TotalDegree);
                            self.A = tform.tdata'; % Putting Matlab's output into my format
                            clear tform
                            tform = cp2tform(StartData',EndData','polynomial',self.TotalDegree);
                            self.Ainv = tform.tdata'; % Putting Matlab's output into my format
                            clear tform
                        catch 
                            disp('FRETmap class: No built-in Matlab polynomial transformation function available; computing manually.')
                            self.Kind = 'Poly';
                            % Just do it manually
                            x = self.EnumerateMonomials(StartData,self.MaxDegree,self.TotalDegree);
                            self.A = self.CalcPolyTform(x,EndData);
                            clear x
                            x = self.EnumerateMonomials(EndData,self.MaxDegree,self.TotalDegree);
                            self.Ainv = self.CalcPolyTform(x,StartData);
                            clear x
                        end
                    end
                end
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
            if strcmpi(self.Kind,'affine') || strcmpi(self.Kind,'AffineMatlab')
                x = [oldspots; ones(1,size(oldspots,2))];
                y = self.A*x;
                newspots = y(1:2,:);
            elseif strcmpi(self.Kind,'poly')
                OldMonomials = self.EnumerateMonomials(oldspots,self.MaxDegree,self.TotalDegree);
                newspots = self.A*OldMonomials;
            else % Matlab poly
                try 
                    tform = images.geotrans.PolynomialTransformation2D(self.A(1,:),self.A(2,:));
                    newspots = transpose(transformPointsInverse(tform,oldspots'));
                    clear tform
                catch
                    % This is so silly, but it's the only way I could
                    % figure out how to use maketform with a polynomial:
                    tempdata1 = rand(2,10);
                    tempdata2 = rand(2,10);
                    temptform = cp2tform(tempdata1',tempdata2','polynomial',2);
                    tform = maketform('custom',2,2,[],temptform.inverse_fcn,transpose([self.A(1,:);self.A(2,:)]));
                    newspots = transpose(tforminv(tform,oldspots'));
                    clear tform tempdata1 tempdata2 temptform
                end
                        
            end
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
            if strcmpi(self.Kind,'affine') || strcmpi(self.Kind,'AffineMatlab')
                M = self.Ainv(:,1:2);
                b = self.Ainv(:,3);
                b = repmat(b,1,size(oldspots,2));
                newspots = M*(oldspots-b);
            elseif strcmpi(self.Kind,'Poly') 
                OldMonomials = self.EnumerateMonomials(oldspots,self.MaxDegree,self.TotalDegree);
                newspots = self.Ainv*OldMonomials;
            else
                try 
                    tform = images.geotrans.PolynomialTransformation2D(self.Ainv(1,:),self.Ainv(2,:));
                    newspots = transpose(transformPointsInverse(tform,oldspots'));
                    clear tform
                catch
                    % This is so silly, but it's the only way I could
                    % figure out how to use maketform with a polynomial:
                    tempdata1 = rand(2,10);
                    tempdata2 = rand(2,10);
                    temptform = cp2tform(tempdata1',tempdata2','polynomial',2);
                    tform = maketform('custom',2,2,[],temptform.inverse_fcn,transpose([self.Ainv(1,:);self.Ainv(2,:)]));
                    newspots = transpose(tforminv(tform,oldspots'));
                    clear tform tempdata1 tempdata2 temptform
                end
            end
            if flipoutput
                newspots = transpose(newspots);
            end
        end
        
        % Return a tform structure or affine2d/polynomial class object in
        % Matlab's format, for use with other Matlab functions like imwarp
        function MatlabTform = ReturnMatlabTform(self,direction)
            if strcmpi(direction,'fwd') || (strcmpi(direction,'inv') && ...
                    (strcmpi(self.Kind,'Affine') || strcmpi(self.Kind,'MatlabAffine')))
                Atouse = self.A;
            elseif strcmpi(direction,'inv')
                Atouse = self.Ainv;
            else
                disp(strcat('Direction: ',direction,' not defined for FRETmap class method HistResiduals.'))
                return
            end
            if strcmpi(self.Kind,'Affine') || strcmpi(self.Kind,'MatlabAffine')
                % Due to small numerical errors, the last column of
                % Atouse might not be [0;0;1] exactly like it needs to
                % be.
                Atouse = inv(transpose(Atouse));
                Atouse(:,3) = [0;0;1];
                try 
                    temptform = affine2d(Atouse);
                    if strcmpi(direction,'fwd')
                        MatlabTform = temptform;
                        clear temptform
                    else
                        MatlabTform = invert(temptform);
                    end
                catch
                    temptform = maketform('affine',Atouse);
                    if strcmpi(direction,'fwd')
                        MatlabTform = temptform;
                        clear temptform
                    else
                        MatlabTform = fliptform(temptform);
                    end
                end
            elseif strcmpi(self.Kind,'MatlabPoly')
                try 
                    MatlabTform = images.geotrans.PolynomialTransformation2D(Atouse(1,:),Atouse(2,:));
                catch
                    % This is so silly, but it's the only way I could
                    % figure out how to use maketform with a polynomial:
                    tempdata1 = rand(2,10);
                    tempdata2 = rand(2,10);
                    temptform = cp2tform(tempdata1',tempdata2','polynomial',2);
                    MatlabTform = maketform('custom',2,2,[],temptform.inverse_fcn,transpose([Atouse(1,:);Atouse(2,:)]));
                    clear tempdata1 tempdata2 temptform
                end
            else
                disp('Class FRETmap: method ReturnMatlabTform not defined for Kind = Poly.')
                return
            end
        end
        
        % Plot the residuals from the fit that produced the mapping
        % contained in this object
        function HistResdiuals(self,direction)
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
        % Homemade affine transformation calculator:
        function TformMatrix = CalcAffineTform(pointsStart,pointsEnd)
            % Redundant, but static methods can be called without an
            % instance of the class object:
            if size(pointsStart,1)~=2
                pointsStart = transpose(pointsStart);
            end
            if size(pointsEnd,1)~=2
                pointsEnd = transpose(pointsEnd);
            end
            x = [pointsStart; ones(1,size(pointsStart,2))];
            y = [pointsEnd; ones(1,size(pointsEnd,2))];
            TformMatrix = y*x'/(x*x');
        end
        % Homemade affine inverter, with my own notation that differs from
        % Matlab's
        function TformInv = CalcAffineInv(TformFwd)
            % pretty sure I have to extract M and b, and can't do something
            % fancy with A?
            M(:,1) = TformFwd(1:2,1);
            M(:,2) = TformFwd(1:2,2);
            b = TformFwd(1:2,3);
            TformInv = [inv(M) b];
        end
        % For making homemade polynomials: enumerate all the monomials the
        % user wants in the polynomial fit
        function newStart = EnumerateMonomials(StartVec,MaxDeg,TotDeg)
            i=0; % Indexes the power of x_1 (first variable of the polynomial)
            j=0; % Indexes the power of x_2 (second variable of the polynomial)
                % (The idea here is that the Cartesian coordinates of a
                % spot x are (x_1,x_2).)
            newStart = []; % Can I figure out ahead of time how big x will be,
                % based on MaxDeg and TotDeg?
            while i<=MaxDeg && i+j<=TotDeg
                while j<=MaxDeg && i+j<=TotDeg
                    % Uncomment the next line to display the monomials
                    % being included. NOTE my code enumerates monomials in
                    % a different order than Matlab does, so the parameters
                    % returned will be in a different order than Matlab's
                    % built-in fitgeotrans or cp2tform!
                    %disp(sprintf('x_1^%d x_2^%d',i,j))
                    newStart(end+1,:) = StartVec(1,:).^i.*StartVec(2,:).^j;
                    j = j+1;
                end
                i=i+1;
                j=0;
            end
        end
        % Homemade polynomial transformation calculator
        function TformPoly = CalcPolyTform(StartMonomials,EndPts)
            x = StartMonomials;
            y = EndPts; % EndPts doesn't change to embed polynomial in a linear space
            TformPoly = y*x'/(x*x');
        end
        % Plot, as a scatter plot, a set of data points and their
        % transformation
        function PlotTform(Data1,Data2)
            figure('Position',[200 200 325 625])
            plot(Data1(2,:),0-Data1(1,:),'xg')
            hold on
            plot(Data2(2,:),0-Data2(1,:),'xr')
            hold off
            % ylim([-512 0])
            % xlim([0 256])
        end
    end
end