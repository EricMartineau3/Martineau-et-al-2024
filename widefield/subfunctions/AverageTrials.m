function [av_S1]  = AverageTrials(S,fps,final_tBin, chanTag, TrialType)
%% Mean of all trials by type
av_S1 = mean(S(:,:,:,TrialType),4,"omitnan");

if sum(TrialType)<size(S,4)
    av_S1(:,:,:,2) = mean(S(:,:,:,~TrialType),4,"omitnan");
end

%% Temporal averaging
binSize = round(fps / final_tBin);
A = zeros(size(av_S1,1),size(av_S1,2),floor(size(av_S1,3)/binSize),size(av_S1,4));

for i = 1:size(A,3)
        C = av_S1(:,:,(binSize*(i-1)+1):i*binSize,:);
        A(:,:,i,:) = mean(C,3, "omitnan");
end
clear C D

%Save images for Fiji
binTag = int2str(final_tBin);

im1 = mat2gray(A(:,:,:,1));
saveTifStack(im1,strcat(chanTag, "_Real_", binTag, "hzBin"));
if size(A,4)>1
    im1 = mat2gray(A(:,:,:,2));
    saveTifStack(im1,strcat(chanTag, "_Sham_", binTag, "hzBin"));
end