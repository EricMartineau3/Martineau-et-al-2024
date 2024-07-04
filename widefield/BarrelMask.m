function [mask] = BarrelMask(data, Threshold, Start, End, SeriesID)
% User assisted function to find regions of activity (ROI) based on signal intensity
% during a specified time window
%
%   INPUTS:
%       - data : xyt matrix of the channel used to find the ROI. Usually
%       F2corr or F1corr. Consider excluding areas where the ciment is
%       visible using drawassisted if signals are detected in the ciment.
%       - Threshold: threshold over which signal is included in the ROI. Value between 0 and 1.
%       - Start : frame number corresponding to the start of the time window
%       - End : frame number corresponding to the end of the time window
%       - SeriesID: string identifier of the recording being analyzed
%
%   OUTPUTS:
%       -mask : xy mask of the identified ROI
%
% Written by Eric Martineau and Antoine Malescot - Universite de Montreal


%%
data = data(:,:,(Start+1):End);
data = mean(data,3,'omitnan');

maxVal = max(data,[],'all');
mask = data > (maxVal*Threshold);


%Prompt user to select right ROI
imshow(~mask);
title(SeriesID);
[column, row] = ginput(1);
labelNumber = mask(row, column);
extractedObject = ismember(mask, labelNumber);
mask = bwselect(mask,column,row); 

clear max 
close all


