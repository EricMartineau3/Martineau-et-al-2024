function [dataFilt, data] = normalize_signal(data,infos,Trial,Parameters,fSize,deltaonly)

%% DESCRIPTION %%
%Subfuction to simplify and streamline the preprocessing loop

%Function to open channel and make dF/F
% INPUTS : 
%   - data : frames in X-Y-T format
%   - infos: metadata for the 

%deltaonly : if 1, perform delta from baseline only, if 0 perform dI/I0

%Written by Eric Martineau and Antoine Malescot - Universite de Montreal

%%

%%% Average baseline %%%
preStim = fix((infos.Freq)*((Trial.StimStart - Trial.Start)/1000));
baseline = data(:,:,(1 + preStim - round(Parameters.base_perc*preStim)) :preStim);
S0 = mean(baseline,3,"omitnan");
clear baseline

%%%% dI/I %%%
if deltaonly == 0
    data = (data - S0)./S0;
elseif deltaonly == 1
    data = data-S0;
end

%%% gaussian filter %%%
nanFlag = isnan(data);
data(nanFlag) = 0; %remove nan before gaussian filt to prevent enlargement of NaN value areas.
dataFilt = imgaussfilt(data,Parameters.Sigma);
data(nanFlag) = NaN; %re-exclude pixels
dataFilt(nanFlag) = NaN;
clear S0

%%% Adjust number of frames %%%
if isempty(fSize) == 0
    if fSize < size(dataFilt,3)%remove last frame if too many frames (happens sometimes due to clock jitter).
        dataFilt = dataFilt(:,:,1:fSize);
        data = data(:,:,1:fSize);
    elseif fSize > size(dataFilt,3)
        dataFilt = cat(3,dataFilt,repmat(zeros(size(dataFilt,[1,2])),1,1,fSize-size(dataFilt,3)));
        data = cat(3,data,repmat(zeros(size(data,[1,2])),1,1,fSize-size(data,3)));
    end
end


