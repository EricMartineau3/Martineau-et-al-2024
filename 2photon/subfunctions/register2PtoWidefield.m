function [xWF, yWF,xy, z, fig1,fig2,fig3,fig4] = register2PtoWidefield(wfParameters, stack2P, vesselChan, depth, barrelCentroids,stack_pixelSize, wf_pixelSize, maxHist,standardROI)
%% DESCRIPTION %%
% Script to align 2P z-stack to widefield surface vessel and map features
% (ROIs,vessels,etc.) onto widefield image.
%
%%%INPUTS%%%
%  wfParameters: structure containing all analysis parameters. The relevant
%             ones for this function are :
%             - Parameters.windowMap : xy grayscale image of the window (usually from the green reflectance to see surface vasculature)
%             - Parameters.ROI_IDX: {nMasks x 1} containing the string identifiers of the recording the ROI has been pulled from
%             - Parameters.Masks: xy-nMasks matrix of masks extracted using BarrelMask or the standard ROI alignement%
%              

%   stack2P = XYZ stack of the region imaged in 2P, must include surface
%             vessels.

%   vesselChan = Index of vessel channel in the stack.
%
%
%   depth = # of slices in 2P stack to make projection and map onto
%           widefield map.
%   
%   barrelCentroids = coordinates of barrelCentroids, with barrelCentroids(:,1) 
%                     corresponding to x-axis values and  barrelCentroids(:,2) 
%                     corresponding to y-axis values. Optional, will
%                     not plot if empty
%
%   stack_pixelSize = size of pixels in stack2P
%
%   wf_pixelSize = size of pixels in wf images

%%%OUTPUTS%%%
%   xWF = x-axis coordinates of identified features (ROIs, vessels, etc) in
%       the widefield reference frame
%
%   yWF = y-axis coordinates of identified features (ROIs, vessels, etc) in
%       the widefield reference frame
%
%   xy = xy-coordinates of identified features (ROIs, vessels, etc) in
%       the 2P reference frame
%
%   zsize = depth of stack in um;
%
%   fig1 = Widefield image with barrel masks and centroids (optional)
%   fig2 = Alignement of 2P stack to surface vessel image
%   fig3 = location of features (vessels or ROIs) in the 2P stack
%   fig4 = location of features (vessels or ROIs) in the widefield image

% Written by Éric Martineau - Universite de Montreal

%% Parameters
sf = stack_pixelSize/wf_pixelSize; %approximation of the scale factor

%% Plot WF map with masks if they exists
wfMap = wfParameters.windowMap*2;

%Show barrel masks
MaskMerge = zeros(size(wfParameters.Masks,1),size(wfParameters.Masks,2),3);
if standardROI == 1
    for i = 1:size(wfParameters.StandardMasks,3)
        im = uint8(255 * mat2gray(wfParameters.StandardMasks(:,:,i)));
        im = ind2rgb(im,custom_cmap(i));
        MaskMerge = max(MaskMerge, im);
        clear img
    end
else
    for i = 1:size(wfParameters.Masks,3)
        im = uint8(255 * mat2gray(wfParameters.Masks(:,:,i)));
        im = ind2rgb(im,custom_cmap(i));
        MaskMerge = max(MaskMerge, im);
        clear img
    end
end
map_rgb = mat2gray(wfParameters.windowMap)*0.75;
map_rgb = cat(3,map_rgb, map_rgb, map_rgb);
MaskMerge = max(MaskMerge,map_rgb);
fig1 = figure();
wfIm = imshow(imresize(MaskMerge,2));

if isempty(barrelCentroids) == 0 %plot centroids if passed to function
    hold(wfIm.Parent,'On')
    for k = 1:size(wfParameters.Masks,3)  
        plot(wfIm.Parent,barrelCentroids(k,1)*2,barrelCentroids(k,2)*2,'y+','MarkerSize',6,'LineWidth',1);
    end
    hold(wfIm.Parent,'Off')
end
clear map_rgb

%% Align 2P stack to window map
% Make surface projection
substack = squeeze(max(stack2P(:,:,1:depth,vesselChan),[],3));
substack = imcomplement(imresize(substack,sf)); %scale down image and invert scale

% Register substack to widefield map
answer = 'no'; %initialize answer
fig2 = figure();
while strcmp(answer,'Continue')~= 1
    [mp,fp] = cpselect(mat2gray(substack),wfMap,'Wait',true); %control point selection
    t = fitgeotrans(mp,fp,'similarity'); %transformation calculation
    Rfixed = imref2d(size(wfMap)); %reference frame %%see if necessary later
    reg_substack = imwarp(substack,t,'OutputView',imref2d(size(wfMap)));    
    alignIm = imshowpair(imresize(wfMap,2), imresize(reg_substack,2),'Scaling','independent');      
    if isempty(barrelCentroids) == 0
        hold(alignIm.Parent,'On')
        plot(alignIm.Parent,barrelCentroids(:,1)*2,barrelCentroids(:,2)*2,'y+','MarkerSize',6,'LineWidth',1);
        hold(alignIm.Parent,'Off')
    end
    
    answer = questdlg("Is the alignement okay?","Alignement check","Continue","Redo alignement",1);
end

%%calculate true scale factor from registration
u = [0 1]; 
v = [0 0]; 
[sX, sY] = transformPointsForward(t, u, v); 
dx = sX(2) - sX(1); 
dy = sY(2) - sY(1); 
angle = (180/pi) * atan2(dy, dx);
scale = 1 / sqrt(dx^2 + dy^2);

true_sf = sf/scale

clear answer u v x y dx dy angle scale substack depth mp fp chan 

%% Select vessels
%open sliceviewer and obtain handle
fig3 = figure();
if size(stack2P,4)>1 %if multichannel, convert to rgb 
    stack2P = cat(4,stack2P(:,:,:,vesselChan),stack2P(:,:,:,3-vesselChan),zeros(size(stack2P,1,2,3)));
end
s = sliceViewer(mat2gray(stack2P,[0 maxHist]));
sAx = getAxesHandle(s);

%obtain xy coordinates of vessels
xy = [];
answer = 'no'; %initialize answer
while strcmp(answer,'Done')~= 1 
    point = drawpoint(sAx);
    xy = [xy; point.Position]; 
    answer = questdlg("Add another feature?","Feature selector","Add feature","Done",1);    
end
%Plot vessel makers on MIP and label them
fig4 = figure();
mip = imshow(squeeze(mat2gray(max(stack2P,[],3))));
hold(mip.Parent,'On')
plot(mip.Parent,xy(:,1),xy(:,2),'bo','MarkerSize',6,'LineWidth',2);
label = strcat("Feature ", cellstr(num2str((1:size(xy,1))')));
text(xy(:,1)+3, xy(:,2)-3, label, 'Fontsize', 10,'Color', 'w');
hold(mip.Parent,'Off')

clear answer point label

%transpose vessel xy-coordinates into the widefield reference frame
sX = xy(:,1)*sf; %downscale to widefield resolution
sY = xy(:,2)*sf;
[xWF, yWF] = transformPointsForward(t, sX, sY);
clear sX sY

%plot results and save them
hold(wfIm.Parent,'On')
hold(alignIm.Parent,'On')
plot(wfIm.Parent,xWF*2,yWF*2,'w+','MarkerSize',2,'LineWidth',1);
plot(alignIm.Parent,xWF*2,yWF*2,'w+','MarkerSize',2,'LineWidth',1);
hold(wfIm.Parent,'Off')
hold(alignIm.Parent,'Off')

