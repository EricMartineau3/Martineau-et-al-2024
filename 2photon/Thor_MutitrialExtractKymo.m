%% Temporary ana routine for 2PI files %%
% Routine to extract kymographs form multi trial acquisition.
% Data must be saved as a series of T x X frames.
%
% The experiment folder must contain individual image folders for each trial,
% each containing a single multi-dimension image file named 'Image_scan_1_region_0_0.tif'
% and an 'Experiment.xml' file containing the acquisition metadata. 
% Image Files must should be saved in OME TIFF format as best practice, but the script is compatible with raw .tiff. 
%
% The script also requires  a matching set of trial file containing the trial infos.
% Trial: structure containing the trial infos. The relevant parameters for this function are:
%             - Trial.Type: string identifier of the trial type ("real" or .sham")
%             - Trial.Start: time of the trial start in us
%             - Trial.StimStart: time of the stimulation start, in us.
%             - Trial.StimEnd: time of the stimulation end, in us.
%             - Trial.End: time of the end of the trial, in us. 
%
%
% Written by Eric Martineau and Antoine Malescot - Universite de Montreal
%
% Notice:
%   This file is provided as a reference of the scripts used to process data
%   aquired with a Thorlab Bergamo II system, using ThorImage versions 4.1.2021.5031 to 4.3.2022.4062)
%   in Martineau et al 2024. This script is made to process files saved in the 
%   OME TIFF format. We make no claims as to the functionality or intended application
%   under a different set of conditions and the user assumes all responsibility for its use. 
%
% Written by Eric Martineau
%% Extract 
discardTrials = []; %number of the trials that need to be excluded
rbcChan = 2; %channel with vessel labeling
stimOnset = 4; %time of stim onset in seconds

%if real stim are stims and shams are shams
% RealID = "Real"; 
% ShamID = "Sham";

%if real=C2stim and shams=D2stim ares stims
RealID = "C2"; 
ShamID = "D2";

oemTiff = 1; % Value of 0 if images saved in raw tiff (series of images), value of 1 if saved in oem tiff format (one image contains stack)

%% Import image files and trial data
%Get file path
path = uigetdir(pwd,'Select image folder');
[trialFile, trialPath] = uigetfile('.mat','Select Trial file');

%List image files
cd(path);
D = dir;
D = D([D.isdir]);
flag = contains({D.name},'.');
D = D(~flag); 
flag = contains({D.name},'AQuA');
D = D(~flag);
flag = contains({D.name},'Vasometrics');
D = D(~flag); 
flag = contains({D.name},'Fourrier');
D = D(~flag); 

%load trial file
load(fullfile(trialPath, trialFile),'Trial');

%remove discarded trials
for i = 1:length(discardTrials)
    D(discardTrials(i,1)-i+1) = [];
    Trial(discardTrials(i,1)-i+1) = [];
end
clear i

%Create Data file structure
AllData = struct('ID',cell(length(D),1),'Frames', cell(length(D),1),'Type',cell(length(D),1),'Rewarded',cell(length(D),1),'metadata',cell(length(D),1),'xmlMeta',cell(length(D),1));

%Import and alignement loop
for k = 1:length(D)
    %Import
    if oemTiff == 1
        [stack, metadata, xmlMeta] = ThorImport(fullfile(path,D(k).name), 'Image_scan_1_region_0_0.tif');
    elseif oemTiff == 0
        [stack, metadata, xmlMeta] = ThorImport_noOEM(fullfile(path,D(k).name));
    end
    disp(strcat('Imported stack', num2str(k)));
    
    %Stiching frames together
    dims = size(stack);
    nStack = zeros(dims(1)*dims(3), dims(2), size(stack,4),'uint16');
    for i = 1:dims(3)
        tIdx1 = (i-1)*dims(1)+1;
        tIdx2 = i*dims(1);
        nStack(tIdx1:tIdx2,:,:) = squeeze(stack(:,:,i,:));
    end
    clear stack
    disp(strcat('Stiched stack', num2str(k)));
    
    %Store data
    AllData(k).ID = D(k).name;
    AllData(k).Frames = nStack;
    AllData(k).metadata = metadata;
    AllData(k).xmlMeta = xmlMeta;
    
    %Import relevant trial data
    if Trial(k).Type == 'Real'
        AllData(k).Type = RealID;
    elseif Trial(k).Type == 'Sham'
        AllData(k).Type = ShamID;
    else
        disp('Trial type error');
    end
    AllData(k).Rewarded = Trial(k).Rewarded;
    %Add any data you wish to filter trials with here (early licks, movement
    %etc.)
end 
clear k

% Save everything
fname = strcat(extractBefore(trialFile,'stack'),'AllData');
save(fullfile(path,fname),'AllData','discardTrials','-v7.3');

%% Fourrier transform processing
vesselID = 'Vessel1'; %modifiy if multiple vessel segments in same line and run multiple times
nStdev = 3; % Input a number to exclude x nb of stdev above or below mean, input NaN if you want to keep all points 
stall = "none";% Method to exclude outliers or detect stalls:
            % "auto" will use the matlab findchangepts to detect stalls with user input
            % "manual" will allow the user to select stalls
            % "none" will skip %this option was used in Martineau et al 2024
linebin = 1; %number of lines to bin to increase SNR. Use 1 for GG scans, 10 for GR scans
dT_win = 0.2; %window over which speed and flux are calculated in seconds. Default = 0.2, increase if vessel very slow.
dT_step = 0.2; %step for window in seconds. Default = 0.2, increase if vessel very slow.

mkdir Fourrier
out = strcat(pwd,'/Fourrier');
mkdir Fourrier FluxOutput
linerate = AllData(1).metadata.FrameRate*AllData(1).metadata.SizeY/linebin;
clear AllSpeed AllFlux

for k = 1:length(AllData)    
    dat = AllData(k).Frames(:,:,rbcChan);
    [dat] = lineBinning(dat,linebin);  
    if k == 1
        [Time, Speed, Parameters, X, maxFreqT, maxFreqX,xStall] = speed_fourier(dat, linerate, AllData(k).metadata.PixelSize, dT_win, dT_step, [],[], [], nStdev,stall);
        [flux, x, kmean, thresh,fig] = rbcCount(dat,[],[], linerate, dT_win, dT_step, X,[]);
    else
       [~, Speed,~ ,~ ,~,xStall] = speed_fourier(dat, linerate, AllData(k).metadata.PixelSize, dT_win, dT_step, X, maxFreqT, maxFreqX, nStdev,stall);
       [flux, ~, ~,fig] = rbcCount(dat,kmean, thresh, linerate, dT_win, dT_step, X,[]);
    end

    %Save flux measurments figure
    savefig(fig,fullfile(out,'FluxOutput',strcat(vesselID,"_Flux",num2str(k))));
    close all

    %Calculate stall duration and exclude stalls from flux
    stallDur = 0;
    for i = 1:size(xStall,1)
       flux(xStall(i,1):xStall(i,2)) = NaN;       
       stallDur = stallDur + Time(xStall(i,2))-Time(xStall(i,1));%total duration of stalls       
    end

    %Store data
    AllSpeed(:,k) = Speed;
    AllFlux(:,k) = flux;
    AllStalls{1,k} = xStall;
    AllStalls{2,k} = stallDur;
    clear dat
end
Parameters.stimOnset = stimOnset;
Parameters.nStdev = nStdev;
Parameters.linebin = linebin;
Parameters.kmean = kmean;
Parameters.thresh = thresh;
Parameters.stall = stall;

%% Density, Normalization and figure creation
out = strcat(pwd,'/Fourrier'); 

% Density calculation
maxD = min(length(AllSpeed), length(AllFlux)-1);
AllDensity = zeros(maxD, size(AllSpeed,2));
AllDensity = AllFlux(2:(1+maxD),:)./AllSpeed(1:maxD,:);
AllDensity(AllDensity < 0) = NaN;
AllDensity(AllDensity > 400) = NaN; %~5um cell, theoretical max in a mm = 200 so 400 to exclude absurd data points

% Normalization and figure creation
baselineIdx1 = find(Time(1,:) < stimOnset);
baselineIdx2 = find(x(:,1) < stimOnset);
RelativSpeed = AllSpeed./ mean(AllSpeed(1:max(baselineIdx1),:),'omitnan');
RelativFlux = AllFlux./ mean(AllFlux(1:max(baselineIdx2),:),'omitnan');
RelativDensity = AllDensity./mean(AllDensity(1:max(baselineIdx1),:),'omitnan');

flag = contains([AllData.Type],RealID);

%Split trials
Y1 = AllSpeed(:,flag);
Y2 = AllSpeed(:,~flag);
Z1 = AllFlux(:,flag);
Z2 = AllFlux(:,~flag);
D1 = AllDensity(:,flag);
D2 = AllDensity(:,~flag);

%Plot speed
fig = figure();
subplot(3,4,1)
plot(permute(Time,[2 1]), permute(Y1,[2 1]));
title(strcat('IndivTraces - ',RealID));
xlabel('Time (s)');
ylabel('RBC velocity (mm/s');

subplot(3,4,2)
plot(permute(Time,[2 1]), permute(Y2,[2 1]));
title(strcat('IndivTraces - ',ShamID))
xlabel('Time (s)');
ylabel('RBC velocity (mm/s');

subplot(3,4,3)
Y1 = mean(Y1,2,'omitnan');
Y2 = mean(Y2,2,'omitnan');
plot(permute(Time,[2 1]), permute(Y1,[2 1]), permute(Time,[2 1]), permute(Y2,[2 1]));
title('Average Traces');
xlabel('Time (s)');
ylabel('RBC velocity (mm/s');
AvSpeed = [Y1 Y2];
legend(RealID, ShamID);

subplot(3,4,4)
Y1 = RelativSpeed(:,flag);
Y2 = RelativSpeed(:,~flag);
Y1 = mean(Y1,2,'omitnan');
Y2 = mean(Y2,2,'omitnan');
plot(permute(Time,[2 1]), permute(Y1,[2 1]), permute(Time,[2 1]), permute(Y2,[2 1]));
title('Average Traces - Relative');
xlabel('Time (s)');
ylabel('RBC velocity (mm/s');
legend(RealID, ShamID);
RelAvSpeed = [Y1 Y2];

%Plot flux
subplot(3,4,5)
plot(permute(x,[2 1]), permute(Z1,[2 1]));
title(strcat('IndivTraces - ',RealID));
xlabel('Time (s)');
ylabel('RBC flux (RBCs/s');

subplot(3,4,6)
plot(permute(x,[2 1]), permute(Z2,[2 1]));
title(strcat('IndivTraces - ',ShamID))
xlabel('Time (s)');
ylabel('RBC flux (RBCs/s');

subplot(3,4,7)
Z1 = mean(Z1,2,'omitnan');
Z2 = mean(Z2,2,'omitnan');
plot(permute(x,[2 1]), permute(Z1,[2 1]), permute(x,[2 1]), permute(Z2,[2 1]));
title('Average Traces');
xlabel('Time (s)');
ylabel('RBC flux (RBCs/s');
AvFlux = [Z1 Z2];
legend(RealID, ShamID);

subplot(3,4,8)
Z1 = RelativFlux(:,flag);
Z2 = RelativFlux(:,~flag);
Z1 = mean(Z1,2,'omitnan');
Z2 = mean(Z2,2,'omitnan');
plot(permute(x,[2 1]), permute(Z1,[2 1]), permute(x,[2 1]), permute(Z2,[2 1]));
title('Average Traces - Relative');
xlabel('Time (s)');
ylabel('RBC flux (RBCs/s');
legend(RealID, ShamID);
RelAvFlux = [Z1 Z2];

%Plot Density
subplot(3,4,9)
plot(permute(Time(:,1:maxD),[2 1]), permute(D1,[2 1]));
title(strcat('IndivTraces - ',RealID));
xlabel('Time (s)');
ylabel('RBC density (RBCs/mm)');

subplot(3,4,10)
plot(permute(Time(:,1:maxD),[2 1]), permute(D2,[2 1]));
title(strcat('IndivTraces - ',ShamID))
xlabel('Time (s)');
ylabel('RBC density (RBCs/mm)');

subplot(3,4,11)
D1 = mean(D1,2,'omitnan');
D2 = mean(D2,2,'omitnan');
plot(permute(Time(:,1:maxD),[2 1]), permute(D1,[2 1]), permute(Time(:,1:maxD),[2 1]), permute(D2,[2 1]));
title('Average Traces');
xlabel('Time (s)');
ylabel('RBC density (RBCs/mm)');
AvDensity = [D1 D2];
legend(RealID, ShamID);

subplot(3,4,12)
D1 = RelativDensity(:,flag);
D2 = RelativDensity(:,~flag);
D1 = mean(D1,2,'omitnan');
D2 = mean(D2,2,'omitnan');
plot(permute(Time(:,1:maxD),[2 1]), permute(D1,[2 1]), permute(Time(:,1:maxD),[2 1]), permute(D2,[2 1]));
title('Average Traces - Relative');
xlabel('Time (s)');
ylabel('RBC density (RBCs/mm)');
legend(RealID, ShamID);
RelAvDensity = [D1 D2];


savefig(fig,fullfile(out,strcat("rbcPlot_",vesselID)));
save(fullfile(out,strcat('Fourier_',vesselID)),'Time', 'AllSpeed','RelativSpeed','AvSpeed','RelAvSpeed', 'x', 'AllFlux', 'RelativFlux', 'AvFlux', 'RelAvFlux', 'AllDensity', 'RelativDensity', 'AvDensity', 'RelAvDensity','Parameters', 'linerate');
