%Function to extract vessel diameter using vasometrics on a series of 2P
%time series (t-stacks).
%
% AllData : structure containing all the aligned videos and metadata, obtained using Thor_MultitrialExtract2
% VasChan : channel number of vessel labeling (i.e. Alexa680)
% VesselID : text string to identify vessel
% Run batch process on one vessel segment at the time using different
% vesselID
%
% Written by Eric Martineau, adapted by Antoine Malescot - Universite de Montreal

%% Inputs
%runTBin = 12;
VasChan = 2;
VesselID = 'Vessel2';
stimOnset = 4;
penVessel = 0; %boolean variable, assign a value of 1 if the vessel you want to measure a penetrating vessel.

%if real=C2stim and shams=D2stim ares stims
RealID = "C2"; 
ShamID = "D2";

%% Create Vasometrics output fodler
root = pwd;
mkdir Vasometrics;

out = strcat(root,'\Vasometrics\');

%% Batch vasometrics cmd and average table %%
close all
CLspace = 1; %space between crosslines, in micrometer
minDist = 0; %0 = take the maximal fwhm, 1= take minimal fwhm % 0 was used in Martineau et al 2024
cLength = 40; % For standard method: input crossline length here if automatic not working. Empty array ([]) if you want auto.
cLengthR = 20; % For rosas method: length (in pixels) to add on each side of the diameter to determine crossline length.
alpha = 10; % Rotation angle between crossline for penetrating vessels, min value of 10o. 

% Preallocation
diam = zeros(size(AllData,1), size(AllData(1).Frames,3));

%Loop
for k = 1:length(AllData)
    % Run through Vasometrics
    disp(strcat("Analyzing trial",num2str(k))); 
    if k == 1
        if penVessel == 0
            [diam(k,:), fwhms, CLs, fig] = Vasometrics(AllData(k).Frames(:,:,:,VasChan),AllData(k).metadata, {}, CLspace, minDist, cLength,[]);
        elseif penVessel == 1
            [diam(k,:), fwhms, CLs, fig] = Vasometrics_penetrating(AllData(k).Frames(:,:,:,VasChan),AllData(k).metadata, {}, minDist, cLengthR,alpha,[]);
        end
        savefig(fig,fullfile(out,strcat(VesselID,'Crosslines')));
    else
        if penVessel == 0
            [diam(k,:), fwhms, ~, ~] = Vasometrics(AllData(k).Frames(:,:,:,VasChan),AllData(k).metadata, CLs, CLspace, minDist, cLength,[]);
        elseif penVessel == 1
            [diam(k,:), fwhms, ~, ~] = Vasometrics_penetrating(AllData(k).Frames(:,:,:,VasChan),AllData(k).metadata, CLs, minDist, cLengthR,alpha,[]);
        end
    end
    FWHMS{k} = fwhms;
end

%Store av results
Parameters = struct('FinalFramerate', AllData(1).metadata.FrameRate,'CLspace',CLspace,'VesselChan',VasChan);
save(fullfile(out,strcat('diameters',VesselID,'.mat')),'diam','CLs','FWHMS','Parameters','-mat');

%% Plot diameter and save the figure
trialFlag = [AllData.Type];

fBin = round(0.2/ (1/AllData(1).metadata.FrameRate)); %number of frames in 200 ms
nBin = round(length(diam)/fBin);
bDiam = zeros(size(diam,1),nBin);
t = zeros(1,nBin);
for i = 1:nBin %bin Diameters
    bDiam(:,i) = mean(diam(:,((i-1)*fBin+1):i*fBin),2); 
    t(1,i) = (i*fBin)/AllData(1).metadata.FrameRate; %timepoint = end of bin
end

fig = figure('Position', [256,128,1280,768]);
avDiam = zeros(length(bDiam),2);
relAvDiam = zeros(length(bDiam),2);

figure(fig)
Y1 = bDiam(trialFlag == "Real",:);
% Y2 = Y1;
Y2 = bDiam(trialFlag ~= "Real",:);
subplot(2,2,1);
plot(permute(t,[2,1]), permute(Y1,[2,1]));
title(strcat('Indiv Trials ',RealID));
xlabel('Time (s)');
ylabel('Diameter (um)');        

subplot(2,2,2);
plot(permute(t,[2,1]), permute(Y2,[2,1]));
title(strcat('Indiv Trials ',ShamID));
xlabel('Time (s)');
ylabel('Diameter (um)');       

subplot(2,2,3);
avDiam(:,1) = mean(Y1,1,'omitnan');
avDiam(:,2) = mean(Y2,1,'omitnan');
plot(permute(t,[2,1]), avDiam(:,1), permute(t, [2,1]), avDiam(:,2));
title('Average');
xlabel('Time (s)');
ylabel('Diameter (um)');
legend(RealID, ShamID);

baselineIdx = find(t(1,:) < stimOnset);
relAvDiam = avDiam./ mean(avDiam(1:max(baselineIdx),:),1,"omitnan")-1;
subplot(2,2,4);
plot(permute(t,[2,1]), relAvDiam(:,1)*100, permute(t, [2,1]), relAvDiam(:,2)*100);
title('Average - Relative');
xlabel('Time (s)');
ylabel('Diameter (%)');
legend(RealID, ShamID);


sgtitle(fig,VesselID);
savefig(fig,fullfile(out,strcat("diameterPlot_",VesselID)));

save(fullfile(out,strcat('diameters',VesselID,'.mat')),'bDiam','avDiam','relAvDiam','t','-mat','-append');