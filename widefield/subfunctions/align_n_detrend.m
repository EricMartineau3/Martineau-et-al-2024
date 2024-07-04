function [data] = align_n_detrend(chan_name,infos,tform,moveflag,bounds)

%% DESCRIPTION %%
%Subfuction to simplify and streamline the preprocessing loop

% Written by Eric Martineau and Antoine Malescot - Universite de Montreal

%% Open a channel %%
fid = fopen(chan_name,'r+');
data = fread(fid, Inf, 'single=>double');
data = reshape(data,infos.datSize(1,1),infos.datSize(1,2),infos.datLength);  

%% Register frames and discard focus changes %%
fixedRefObj = imref2d(size(data(:,:,1)));
%Registration
if isempty(tform) == 0
    parfor i = 2:infos.datLength
        data(:,:,i) = imwarp(data(:,:,i),tform{i},'linear','Outputview',fixedRefObj);
    end

    % Discard focus changes
    data = reshape(data,[],infos.datLength);
    data(:,moveflag) = NaN; %give a NaN value to focus changes
else
    data = reshape(data,[],infos.datLength);
end

%% Detrend channel %%
% Exclude stim and NaN frame and concatenate the rest
x = linspace(1/infos.Freq,infos.datLength/infos.Freq,infos.datLength); %timepoints

[T, ~,gof] = nl_detrend_JacobianRevolt_2(x, data, bounds,chan_name); %exponential fit on average video to extract trend in LED power changes. Can be done in pixelwise fashion but process time is much slower
plotDetrend(x,data,T,chan_name,gof);
data = data./T;
data = reshape(data,256,256,infos.datLength);
clear T gof

%Method 2
edgeflag = data==0; % identify pixels with a 0 value due to imwarp edging
edgeflag = repmat(max(edgeflag,[],3),1,1,infos.datLength); %exclude any pixel with a 0 value at any timepoint.Otherwise causes an issue with normalization.
data(edgeflag) = NaN;

% Save data
data = reshape(data,[],1);
frewind(fid);
fwrite(fid, data ,'single');
fclose(fid);