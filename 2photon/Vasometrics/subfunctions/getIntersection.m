function [val] = getIntersection(centerX, centerY, radius, slope, yInt)

%Vasometric function get coordinates of crossline intersection with main
%line. Input the last coordinates and function finds the next start point
%at a distance of 'radius'
% Inputs :
%   - centerX : coordinates of start point of x-axis
%   - centerY : coordinates of start point of y-axis
%   - radius : displacement (conjugate of x and y)
%   - slope : slope of main line segment
%   - yInt : intercept of main line segment
%
% Outputs : 
%   - val : array of two possible values

%% Script
a = power(slope,2)+1;
b = 2 * ((slope * yInt) - (slope * centerY) - centerX);
c = power(centerY,2) - power(radius,2) + power(centerX,2) - (2 * yInt * centerY) + power (yInt,2);
x1 = (-1 * b + sqrt(power(b,2) - 4 * a * c)) / (2 * a);
x2 = (-1 * b - sqrt(power(b,2) - 4 * a * c)) / (2 * a);
val = [x1 x2];
