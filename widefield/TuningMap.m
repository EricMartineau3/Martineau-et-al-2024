function [TunSig] = TuningMap(Sig,Start, End, ROIs) 
% Function to calculate the average spatial tuning (i.e. normalizes the response
% in each pixel to the max response across all stims) in each ROI over a
% specific time window
%
%   INPUTS:
%       - sig: xy-nSeries-t matrix, where the signal from multiple series (i.e whisker stim)
%               can be extracted from the each ROI (barrel).
%       - Start: frame number corresponding to the start of the time window
%       - End: frame number corresponding to the end of the time window
%       - ROIs: xy-nMasks matrix of masks extracted using BarrelMask or the standard ROI alignement
% WARNING - Series and ROIs must be ordered so that adjacent whisker stims/barrel ROIs are in adjacent dimensions (i.e. sig(:,:,1:4,:) = E2:B2 stim and ROI(:,:,1:4) = E2:B2 ROIs). 
%
%   OUTPUTS:
%       - TunSig: nSeries-barrelDistance matrix of tuning indexes. Indexes in TunSig(:,1) correspond to the associated barrel, TunSig(:,2) to the 1st neighbors, etc. 
%
% Written by Eric Martineau and Antoine Malescot - Universite de Montreal

%% Parameters %%
dim = size(Sig);
stimNum = dim(3);
clear A

%% Extract AUC of signal %%
AUC = zeros(dim(1), dim(2), stimNum);

for i = 1:stimNum
    A = Sig(:,:,i,:);
    A = A(:,:,Start:End);
    flag = isnan(A);
    if sum(flag,"all") > 0
         disp("Warning! NaN values founds in AUC area. Consider changing interval or discarding series")
         TunSig = [];
         return
    else
        Q = trapz(A,3);
        Q(Q<0) = 0; %to avoid negative values
        AUC(:,:,i) = Q;
    end
end

%% Identify responding pixels
maxRep = max(AUC,[],3);

SD = std(AUC,1,'all');
respPix = AUC > 2*SD;
respPix = max(respPix,[],3);

%% Compute tuning %%
TunMap = AUC ./ maxRep;
TunMap(isnan(TunMap)) = 0;
TunMap = TunMap .* respPix;

%% Extract tuning signals 
TunSig = zeros(stimNum,stimNum); %% Tuning values stim x barrel(E2 -> B2)

for i = 1 : stimNum
    map = TunMap(:,:,i);
    for j = 1: size(ROIs,3)
        mask = logical(ROIs(:,:,j));
        TunSig(i,j) = mean(map(mask),'all','omitnan');
    end
end

%% Unmix tuning signals
unmixMask = zeros(stimNum);
for i = 1:stimNum
    v = ones(1,stimNum-i+1)*i;
    unmixMask = unmixMask + diag(v,(i-1));   
    if i > 1
         unmixMask = unmixMask + diag(v,-1*(i-1));
    end   
end

clear i v

for i = 1:stimNum
    mask = unmixMask(i,:);    
    A = TunSig(i,:);
    for j = 1:stimNum
        B = mask == j;
        C = A(B);
        C = mean(C,'omitnan');
        TunSig(i,j) = C;
    end
end