% function Map = CalcChannelMapping(points1,points2)
%
% Given a set of points in one channel and a matched set of points in the
% other channel, calculated the matrix that maps one channel to the other.
% IMPORTANT the input order matters because this will calculate the
% transformation from points1 to points2!
%
% Inputs:
% 2xnumspots or numspotsx2 matrices of spots
% NOTE this code forces points1 and points2 to end up as 2xnumspots
% matrices, regardless of their sizes when entered as inputs.  Therefore to
% use the outputted mapping, make sure that points1 and points2 are
% 2xnumspots!
%
% Output:
% A matrix A and a vector b, such that points2 = A*points1 + b
%
% Algorithm:
% We want points1 = A*points1 + b = points2
% if we want to allow translations, reflection, rotation, scaling
% Want to embed an affine transformation (Ax+b) in a linear space (Ax)
% So define 
% x = [points1;1]
% M = [A, b; 0 1]
% Then
% M*x = [A*points1+b; 1] = y
% Now, as with a linear transformation, find
% argmin of M(||y - M*x||_2)^2
% That is, minimize the error of our estimate of y = [points2;1] given an M
% Use the normal equations as conditions of optimality:
% (y-M*x)*x' = 0
% Therefore
% argmin of M = y*x'*(x*x')^(-1)
%
% Steph 9/2013
% Copyright 2013 Stephanie Johnson, University of California, San Francisco

function [A,b] = CalcChannelMapping(points1,points2)

% Make the inputs 2-by-numspots matrices, if they're not already, for
% convenience
if size(points1,1)~=2
    points1 = transpose(points1);
end
if size(points2,1)~=2
    points2 = transpose(points2);
end

x = [points1; ones(1,size(points1,2))];
y = [points2; ones(1,size(points1,2))];

M = y*x'*inv(x*x');

A(:,1) = M(1:2,1);
A(:,2) = M(1:2,2);
b = M(1:2,3);

