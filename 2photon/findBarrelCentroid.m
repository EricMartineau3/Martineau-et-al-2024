%% DESCRIPTION %%
% Analysis script to :
% 1 - find the barrel centroid from each mask
% 2 - align a 2P xyz stacks onto the centroid map 
% 3 - identify vessels in stack and extract xy coordinates
% 4 -  report coordinates and measure the distance between a vessel and the C2 and D2
% centroids. Barrel identity can be modified by changing lines 39, 40. 

% Will prompt users to import widefield barrel ROI from a 'parameters' structure and 2P stack of the
% barrel
%   Parameters: structure containing all analysis parameters. The relevant ones for this function are :
%             - Parameters.windowMap : xy grayscale image of the window (usually from the green reflectance to see surface vasculature)
%             - Parameters.ROI_IDX: {nMasks x 1} containing the string identifiers of the recording the ROI has been pulled from
%             - Parameters.Masks: xy-nMasks matrix of masks extracted using BarrelMask or the standard ROI alignement%
%               Will prompt user to open if empty.

%   stack2P = 3D stack of the region imaged in 2P, must include surface
%             vessels. Will open a ui to open if empty. The prompt assumes
%             the stack is saved as a multi-dimension .tif with the number
%             of channels and slices encoded in the file infos. 


% Written by Éric Martineau - Universite de Montreal

%% Parameters
vesselChan = 2; %2P channel containing the vessels for alignement with widefield
standardROI = 1; %1 to use standard ROIs, 0 to use signal-based ROIs
depth = 60; %depth to include in mip of surface vessels
wf_pixelSize = 16.86; %um per pixel
stack_pixelSize = 0.522; %um per pixel
sf = stack_pixelSize/wf_pixelSize; %approximation of the scale factor

outpath = uigetdir('',"Select output path");
cd(outpath)
mkdir Centroids;

%% Find the barrel centroids
close all
%Import ROIs and window map
mapPath = uigetdir('',"Select the parameters file");
load(fullfile(mapPath,"Parameters.mat"),'Parameters');

%Find centroid
if standardROI == 0
    D2 = find(contains(Parameters.ROI_IDX,"D2")==1);
    C2 = find(contains(Parameters.ROI_IDX,"C2")==1);
    centroids = zeros(size(Parameters.Masks,4),2);
    for k = 1:size(Parameters.Masks,3)
        [row, col] = ind2sub(size(Parameters.Masks(:,:,k)),find(Parameters.Masks(:,:,k) == 1)); %find all pixels in ROI
        centroids(k,:) = [mean(col),mean(row)]; %calculate centroids (mean of xy coordinates) and store them
    end
elseif standardROI == 1
    D2 = find(contains(Parameters.StdROI_order(1,:),"D2")==1);
    C2 = find(contains(Parameters.StdROI_order(1,:),"C2")==1);
    centroids = zeros(size(Parameters.StandardMasks,4),2);
    for k = 1:size(Parameters.StandardMasks,3)
        [row, col] = ind2sub(size(Parameters.StandardMasks(:,:,k)),find(Parameters.StandardMasks(:,:,k) == 1)); %find all pixels in ROI
        centroids(k,:) = [mean(col),mean(row)]; %calculate centroids (mean of xy coordinates) and store them
    end
end

clear row col map_rgb
%Import 2P stack if not passed to 2P function
[file,path] = uigetfile('',"Select 2P stack");
stack_info = imfinfo(fullfile(path,file));
chan = str2double(extractAfter(extract(stack_info(1).ImageDescription,"channels=" + digitsPattern(1)),"channels="));
z = str2double(extractAfter(extract(stack_info(1).ImageDescription,"slices=" + digitsPattern(2,3)),"slices="));
if isempty(chan)
    chan = 1;
end

stack2P = zeros(stack_info(1).Height,stack_info(1).Width,z,chan);
k = 1;
for i = 1 : z %read tiff and rebuild xyzc stack
    for j = 1 : chan
        stack2P(:,:,i,j) = imread(fullfile(path,file), k); %vessels should always be 1st index in 4th dimension
        k = k+1;
    end
end


%% Align 2P stack to window map and create figures
[xWF, yWF, xy, ~, fig1,fig2,fig3,fig4] = register2PtoWidefield(Parameters, [], vesselChan, depth, centroids,stack_pixelSize, wf_pixelSize,2000, standardROI);

% Save images
savefig(fig1,fullfile(outpath,'Centroids', 'CentroidMap'));
savefig(fig2,fullfile(outpath,'Centroids', 'AlignmentMap'));
savefig(fig3,fullfile(outpath,'Centroids', 'VesselID'));
savefig(fig4,fullfile(outpath,'Centroids', 'VesselMap'));

save(fullfile(outpath,'Centroids',"coordinates.mat"),'xWF','yWF','xy')

%% Calculate distance from D2 and C2 and ratio for each point 
d_xy = centroids(D2,:);
c_xy = centroids(C2,:);
vesselDist = zeros(size(xWF,1),2);

for i = 1:size(xWF,1)
    vesselDist(i,1) = pdist([c_xy;[xWF(i) yWF(i)]])*wf_pixelSize;
    vesselDist(i,2) = pdist([d_xy;[xWF(i) yWF(i)]])*wf_pixelSize;
end
CDDist = pdist([c_xy;d_xy])*wf_pixelSize;

save(fullfile(outpath, 'Centroids', 'vesselDistance.mat'), 'vesselDist','centroids','xWF','yWF','CDDist', '-mat');