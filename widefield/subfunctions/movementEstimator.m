function [tform, moveFlag] = movementEstimator(chan_name,info_name,zthreshold,mthreshold)
%% DESCRIPTION %%
% Function that estimates rigid motion from a channel and flags potential
% focus changes for elimination. 

% INPUTS:
%   chan_name : name of .dat file to estimate movement from
%   info_name : info file associate with data file
%   zthreshold : threshold for z-score over which a pixel value is considered significantly changed
%   mthreshold : threshold for the percent of pixels changing in a frame. When more than this number of pixels are significantly changing, a frame is considered as moving.

% OUTPUTS:
%   tform : cell-array of rigid transformation to align all frames to the first frame
%   moveFlag : flag for frames with potential focus changes

% Written by Éric Martineau - Universite de Montreal
% Uses Waitbar for Parfor from: Yun Pu (2023). Waitbar for Parfor (https://www.mathworks.com/matlabcentral/fileexchange/71083-waitbar-for-parfor), MATLAB Central File Exchange. Retrieved July 26, 2023.
%% Open images and infos
fid = fopen(chan_name,'r+');
infos = matfile(info_name,'Writable',false);
data = fread(fid, Inf, 'single=>double');
data = reshape(data,infos.datSize(1,1),infos.datSize(1,2),infos.datLength);

%% Rigid registration %%
% Parameters
[optimizer, metric] = imregconfig("monomodal");
optimizer.GradientMagnitudeTolerance = 1.0e-04;
optimizer.MinimumStepLength = 1.0e-05;
optimizer.MaximumStepLength = 0.0625;
optimizer.MaximumIterations = 100;
optimizer.RelaxationFactor = 0.05;
WaitMessage = parfor_wait(infos.datLength, 'Waitbar', true);

%Calculate tform and align the current channel
fixedRefObj = imref2d(size(data(:,:,1)));
tform = cell(1,infos.datLength);
fixed = data(:,:,1);
aligned = zeros(1,infos.datLength); %ma

disp(strcat("Estimating within-trial movement from ",chan_name));
tic()
parfor i = 2:infos.datLength
    tform{i} = imregtform(data(:,:,i),fixed,"rigid",optimizer,metric);
    data(:,:,i) = imwarp(data(:,:,i),tform{i},'Outputview',fixedRefObj);
    WaitMessage.Send;
    pause(0.002);
end
toc()
WaitMessage.Destroy;

%% Remove noise coming from focus changes
%z score pixel map
data = reshape(data,[],infos.datLength); %for simplicity of manipulation
base = mean(data,2);
SD = std(data,1,2);
z_score = (data-base)./SD;

%identify changing pixels
move_thresh = size(z_score,1)*mthreshold; 
thresh_z = (z_score>zthreshold) + (z_score<-zthreshold); %find pixels with a positive or negative change over the z score threshold

%identifiy moving frames
movingPixels = sum(thresh_z,1); %count changing pixels
moveFlag = movingPixels > move_thresh; %if to many changing pixels, frame is considered as moving

%% plot results
fig = figure();
plot(movingPixels);
ylabel('Number of pixels');
xlabel('Frame number');
title(chan_name);
hold on
yyaxis right
trace = mean(data,1);
plot(trace);
ylabel('Grayscale value');
plot(find(moveFlag == 1),trace(moveFlag),'*r')
yyaxis left
plot(repmat(move_thresh,1,infos.datLength));
ylim([0 5000]);
text(infos.datLength,move_thresh(end),strcat("= ",num2str(ceil(move_thresh))," pixels"));
legend({'Pixels exceeding z-score threshold','Pixel count threshold','Mean pixel value','Discarded frames',});

savefig(fig,fullfile(pwd,'discarded_frames.fig'));
