function [Mask_ROI] = ROI_standardizationApp(Parameters,order,file_parameters,center_manual)
% Function to align standardized ROIs onto signal intensity-based ROIs
%
%   INPUTS:
%       - Parameters: structure containing all analysis parameters. The relevant
%             ones for this function are :
%             - Parameters.windowMap : xy grayscale image of the window (usually from the green reflectance to see surface vasculature)
%             - Parameters.ROI_IDX: {nMasks x 1} containing the string identifiers of the recording the ROI has been pulled from
%             - Parameters.Masks: xy-nMasks matrix of masks extracted using BarrelMask or the standard ROI alignement
%       - order : {nSeries x 2} containing the string identifiers of the barrels (i.e. E2:B2) in {:1} and the match ROI identifier in {:,2}
%       - file_parameters: path of the data file. Used for saving purposes
%       - center_manual: boolean. 1 to manually select the ROI centroid, 0 for automatic (used in Martineau et al 2024)
%
%   OUTPUTS:
%       - Mask_ROI: structure containing the standard barrel masks and their IDs
%               - Mask_ROI.Tag: string identifier for each standard ROI
%               - Mask_ROI.Mask: xy logical matrix for each standard mask
%               - Mask_ROI.Centroid: centroids of the signal based-ROI used to align the standard ROI
%
% Written by Antoine Malescot and Éric Martineau - Universite de Montreal
%
% Notice:
%   This file is provided as a reference of the scripts used to process data
%   aquired with a labeoTech OiS2000 system in Martineau et al 2024. This
%   script is dependent on the UMIT toolbox v1.5.8. We make no claims as to 
%   the functionality or intended application under a different set of conditions
%   and the user assumes all responsibility for its use. 

%%
fprintf('============================================================== \n')
fprintf('Standardization \n')
fprintf('Loading Mask... \n')
load('MASK_stock.mat')
load('Whisk_barrel.mat')

%% Preparation

overview = Parameters.windowMap;
original_label = MASK_stock(6).Edge;

size_BCR = 83; %pixel per mm
size_exp = 59.3120; %pixel per mm 

coord_label = [718,451;...
    736,487;...
    756,514;...
    778,532]; %don't touch

% Resize image
facteur_conversion = size_BCR/size_exp;
coord_label = coord_label./facteur_conversion;

for i =1:length(MASK_stock)
    MASK_stock(i).Mask = imresize(MASK_stock(i).Mask,1./facteur_conversion);
    if ~isempty(MASK_stock(i).Edge)
        MASK_stock(i).Edge = imresize(MASK_stock(i).Edge,1./facteur_conversion);
    end
end

for i =1:length(Whisk_barrel)
    Whisk_barrel(i).Mask = imresize(Whisk_barrel(i).Mask,1./facteur_conversion);
end

%% Resize image and stuff
facteur_conversion = size_BCR/size_exp;
label = mat2gray(imresize(original_label,1/facteur_conversion));

%% Verify the order done during experiment
% Struct to save Whiskers barrel masks
Mask_ROI = struct('Tag',[],'Mask',[],'Centroid',[]);
centers_fixed = zeros(4,2);


for i = 1:size(order,2)
   Index = find(contains(Parameters.ROI_IDX, order{2,i}));
   order{3,i} = Index;
    
   if isempty(Index) == 0
       % Take the coordinates of barrels on stim
       reposition = zeros(size(label));
       reposition(1:size(Parameters.Masks(:,:,Index),1),1:size(Parameters.Masks(:,:,Index),2)) = Parameters.Masks(:,:,Index);
       
       % If the option to choose the center manually is selected
       if center_manual == 1
            figure()
            imshow(Parameters.Masks(:,:,Index))
            hold on
            title('Select the center')
            center_tmp = regionprops(Parameters.Masks(:,:,Index),'Centroid');
            plot(center_tmp.Centroid(1),center_tmp.Centroid(2),'r+')
            [col,row] = ginput(1);
            center.Centroid = [col,row];
            close gcf
       else
           center = regionprops(reposition,'centroid');
       end

       Parameters_right_order(i).center = center.Centroid(1);
       new_roi = zeros(size(Parameters.Masks(:,:,Index)));
       new_roi(round(center.Centroid(2)),round(center.Centroid(1))) = 1;
       new_roi = imdilate(new_roi,strel('disk',5)); % Note: radius on image is 8pixels
       Parameters_right_order(i).Whisk = new_roi;
       centers_fixed(i,:) = [center.Centroid(1),center.Centroid(2)];
   end
end

%% Reorganize and extract centroids of selected ROIs
idx = [1:length(order)];
idx = idx(~cellfun(@isempty,order(3,:)));
coord_source = coord_label(idx,:);
centers_target = centers_fixed(idx,:);
clear idx
%% Looking for the right transformation
[regParams,~,~]=absor(coord_source.',centers_target.');
tform = rigid2d(regParams.R.',regParams.t.');

%% Applying transformation on the brain map

EDGE_gene = zeros(size(MASK_stock(1).Mask));

MASK_stock(6).Edge = imtranslate(MASK_stock(6).Edge,[2,0]);

for i=1:length(MASK_stock) % FOR loop on the whole brain Mask

    if (i == 6) % If BC_R 
        MASK_stock(i).Edge = imerode(MASK_stock(i).Edge,strel('square',2)); % Erode slighlty the image
        EDGE_gene = EDGE_gene + double(MASK_stock(i).Edge); % Add to general edge
        MASK_stock(i).Edge = imwarp(MASK_stock(i).Edge,tform,'OutputView',imref2d(size(overview))); % Transform Edge
        MASK_stock(i).Mask = imwarp(MASK_stock(i).Mask,tform,'OutputView',imref2d(size(overview))); % Transform Mask
        MASK_stock(i).center = center_label(mat2gray(MASK_stock(i).Mask)); % center
    else
        EDGE_gene = EDGE_gene + edge(MASK_stock(i).Mask); % Add to general edge
        MASK_stock(i).Mask = imwarp(MASK_stock(i).Mask,tform,'OutputView',imref2d(size(overview))); % Trasnform Mask
        MASK_stock(i).Edge = edge(MASK_stock(i).Mask); % Find the edge 
        MASK_stock(i).center = center_label(MASK_stock(i).Mask); % center
    end
end

EDGE_gene = imwarp(EDGE_gene,tform,'OutputView',imref2d(size(overview))); % Transform general edge
EDGE_gene(EDGE_gene~=0) = 1; % Binarize Edge
C = EDGE_gene+mat2gray(overview); % overlap edge and img window


whisk_done = sum(reshape([Parameters_right_order.Whisk],size(overview,1),size(overview,2),[]),3);

whisk_done(whisk_done>0) = 256;
Mask_aligned(:,:,2) = C.*(1-logical(whisk_done));
Mask_aligned(:,:,3) = C.*(1-logical(whisk_done));
Mask_aligned(:,:,1) = 256.*(logical(whisk_done))+C.*(1-logical(whisk_done));

Mask_aligned(5:10,5:size_exp,1) = 205;
Mask_aligned(5:10,5:size_exp,2) = 256;
Mask_aligned(5:10,5:size_exp,3) = 256;

imwrite(Mask_aligned,fullfile(file_parameters,'Analysis','Mask_aligned.png'))

%% Applying transformation on the whisker barrel mask

roi_overview = zeros(size(overview));
for i = 1:size(order,2) % FOR loop on the 4 whiskers barrel
    Mask_ROI(i).Tag = Whisk_barrel(i).Tag; % Name
    Mask_ROI(i).Mask =  imwarp(Whisk_barrel(i).Mask,tform,'OutputView',imref2d(size(overview))); % Transform
    Mask_ROI(i).Mask = Mask_ROI(i).Mask & Parameters.windowMask; % excludes pixels outside window boundary
    centroid = regionprops(Mask_ROI(i).Mask,'Centroid'); % Find center
    if isempty(centroid)==0 %in case ROI is fully outside window
        Mask_ROI(i).Centroid = [centroid(1).Centroid(1),centroid(1).Centroid(2)]; % New centroid
    end
    roi_overview = roi_overview + Mask_ROI(i).Mask;
end

roi_overview(roi_overview>0) = 256;
Mask_new_roi(:,:,2) = 128.*(logical(roi_overview))+C.*(1-logical(roi_overview));
Mask_new_roi(:,:,3) = 256.*(logical(roi_overview))+C.*(1-logical(roi_overview));
Mask_new_roi(:,:,1) = C.*(1-logical(roi_overview));

imwrite(Mask_new_roi,fullfile(file_parameters,'Analysis','StandardROI_Mask.png'))

fprintf('Finished!\n')
fprintf('============================================================== \n')
end

