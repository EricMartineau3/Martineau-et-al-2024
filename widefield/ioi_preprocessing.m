function [SignalsFilt, FPS] = ioi_preprocessing(Trial, Parameters, hemoCorr, saveRaw, alignment_inter_trial, Discard_Movement, indivTrial)
% Trial: structure containing the trial infos. The relevant parameters for this function are:
%             - Trial.Type: string identifier of the trial type ("real" or .sham")
%             - Trial.Start: time of the trial start in us
%             - Trial.StimStart: time of the stimulation start, in us.
%             - Trial.StimEnd: time of the stimulation end, in us.
%             - Trial.End: time of the end of the trial, in us. 
%
% Parameters: structure containing all analysis parameters. The relevant
%             ones for this function are :
%             - Parameters.HemoCorrChan = string identifier of the channels used for hemodynamic corrections, ex {'green', 'red'};
%             - Parameters.Sigma: value for gaussian filter
%             - Parameters.base_perc: baseline proportion to keep, value between 0 and 1 
%             - Parameters.tBin: temporal binning value
%             - Parameters.sBin: spatial binning value
%             - Parameters.final_tBin: final framerate for video preview, in Hz
%             - Parameters.ComputeHbO: 0 to skip, 1 to compute hbo and hbr
%             - Parameters.FilterSet: char vector defining filterset, ex 'jRGECO'. Consult the HemoCompute_HomeMade function documentation for valid inputs
%             - Parameters.DetrendBounds: time interval to exclude from detrending

% hemoCorr: cell array containing the name of the IOS channels to use for
%            hemodynamic correction. e.g : {'Green' 'Red'} will use both
%            the green and red channel to estimate the hemodynamic
%            correction. Input empty cell array to skip hemodynamic
%            correction
%
% saveRaw: 1 or 0 to determine if .dat and data.mat files should be saved
%
% alignment_inter_trial: Was not used in Martineau et al 2024. Input empty
%                       vector.
%
% Discard_Movement: Was not used in Martineau et al 2024. Input empty
%                       vector.
%
% indivTrial: set value to 1 to save individual trial recordings in the final
%             mat file and set it to 0 to only save average real and sham trials. Set to
%             0 by default on the batch process script to save storage space. Run this
%             script to save the individual trials if you whish to access them.
%
%
% Written by Eric Martineau and Antoine Malescot - Universite de Montreal
%
% Notice:
%   This file is provided as a reference of the scripts used to process data
%   aquired with a labeoTech OiS2000 system in Martineau et al 2024. This
%   script is dependent on the UMIT toolbox v1.5.8. We make no claims as to 
%   the functionality or intended application under a different set of conditions
%   and the user assumes all responsibility for its use. 

%% Preparation
fclose('all');
flag_Rx567 = 0;
InterTrialMovement ={};
count = 0;

p = gcp('nocreate');
if isempty(p)
    parpool(8);
end

%% Gather list of trials
root = pwd;
D = dir;
D = D(3:end);
D = D([D.isdir] == 1);
D = D(~contains({D.name},'Analysis'));

IndivSignalsFilt = struct('Green', cell(1,1),...
    'Amber',cell(1,1),...
    'Red', cell(1,1),...
    'HbO', cell(1,1),'HbR', cell(1,1),...
    'F1', cell(1,1),'F1corr', cell(1,1),...
    'F2', cell(1,1), 'F2corr', cell(1,1),...
    'TrialType',cell(1,length(D)));

% Sort trials in ascending order
if length(D)>1
    id = str2double({D.name});
    [~,I] = sort(id);
    D = D(I);
elseif isempty(D)
    D(1).name = pwd;
end
clear id

% Exclude trials without stim
stimFlag = [Trial.StimStart] ~= 0;
Trial = Trial(stimFlag);
D = D(stimFlag);
clear stimFlag

%% Pre-processing loop %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pre-processing loop for every trial  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(D) %COMPILATION OF TRIALS
    cd(D(i).name);
    disp(strcat('Processing folder #',D(i).name));

    %Clear interp_1.bin thingy
    if exist(fullfile(D(i).folder,D(i).name,'img_interp_1.bin'),'file')==2
        disp('Deleting img_interp1_bin...')
        delete(fullfile(D(i).folder,D(i).name,'img_interp_1.bin'));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         Image classification         %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ImagesClassification(pwd,pwd,Parameters.sBin,Parameters.tBin,0); %new version
    AcqInfoStream = ReadInfoFile(pwd);
    nbColors = sum(cellfun(@(X) contains(X,'Illumination'),fieldnames(AcqInfoStream)));
    chanList = cell(2,nbColors); %%Row 1 = chanID, Row 2 = FPS

    % Rename speckle for Rx_567
    if isfile('Rx_567.dat')
        delete('Rx_567.dat')
        delete('Rx_567.mat')
        movefile speckle.dat Rx_567.dat
    elseif isfile('speckle.dat')
        movefile speckle.dat Rx_567.dat
        movefile speckle.mat Rx_567.mat
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%    Multi-CAM check and alignement    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if AcqInfoStream.MultiCam
        cams2_flag = 1;
        camIllumFlag = cell(nbColors,3);
        for q = 1:nbColors
            chanID = eval(['AcqInfoStream.Illumination' int2str(q) ';']);
            camIllumFlag{q,3} = chanID.FrameIdx;
            camIllumFlag{q,2} = chanID.CamIdx;
            camIllumFlag{q,1} = chanID.Color;
        end
        nbIllumPerCam(1) = sum([camIllumFlag{:,2}]==1);
        nbIllumPerCam(2) = sum([camIllumFlag{:,2}]==2);
        nBCamFrame = length(unique([camIllumFlag{:,3}]));

        %Align cam2 to cam1
        idx_cam1 = find([camIllumFlag{:,2}] == 1);
        idx_cam2 = find([camIllumFlag{:,2}] == 2);
        camreg = load(fullfile(D(1).folder,'tform.mat'));

        chanID = camIllumFlag{idx_cam1(1),1};
        [chan_name, info_name,~] = selectChannelFiles(chanID);
        infos = matfile(info_name,'Writable',false);   
        fid = fopen(chan_name,'r+');
        dat1 = fread(fid, inf, '*single');
        dat1 = reshape(dat1, infos.datSize(1,1), infos.datSize(1,2),[]);
        dat1 = dat1(:,:,1);

        for k = 1:length(idx_cam2)
            chanID = camIllumFlag{idx_cam2(k),1};
            [chan_name, info_name,~] = selectChannelFiles(chanID);
            infos = matfile(info_name,'Writable',false);
            fid = fopen(chan_name,'r+');
            dat = fread(fid, inf, '*single');
            dat = reshape(dat, infos.datSize(1,1), infos.datSize(1,2),[]);
            dat = imwarp(dat, camreg.tform, 'OutputView', imref2d(infos.datSize));

            fclose(fid);
            fid = fopen(chan_name,'w');
            fwrite(fid,dat,'single');
            fclose(fid);

            figure()
            C = imfuse(dat(:,:,1),dat1,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
            imshow(C)
            savefig(['Coregistration' camIllumFlag{idx_cam2(k),1}(1:end-4) '.fig'])
            close gcf
        end
    else
        cams2_flag = 0;
    end   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Image registration and frame discard %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Align and detrend every channel
    for j = 1:nbColors
        chanID = eval(['AcqInfoStream.Illumination' int2str(j) ';']);
        chanList{1,j} = chanID.Color;
        [chan_name, info_name,structField] = selectChannelFiles(chanID.Color);        

        if contains(chan_name,'Rx_567.dat')==1 %flag lime reflectance
            flag_Rx567 = 1;
        end          
        
        %%% Open metadata %%%        
        infos = matfile(info_name,'Writable',true);        
        if cams2_flag == 1
            chanList{3,j} = chanID.CamIdx; %associate camIDX with color channel        
            infos.Freq = round(AcqInfoStream.FrameRateHz/nBCamFrame);  % temporary correction of framerate metadata for 2 cam
        else
            chanList{3,j} = 1;%if 1 cam, all CAM1
        end
        chanList{2,j} = infos.Freq;
        
        %%% Align and detrend channel %%%
        [data] = align_n_detrend(chan_name,infos,{},[],Parameters.DetrendBounds);            
        data = reshape(data,infos.datSize(1,1),infos.datSize(1,2),infos.datLength);          
    end %END of alignement and detrending loop
    clear data    

    if keep_trial
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%        Hemodynamic correction        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(hemoCorr)
            fList = dir([strcat(pwd, filesep) 'fluo*.dat']);

            % Case where multiple fluo channels exist
            if length(fList) > 1
                [NewFluo1, NewFluo2] = HemoCorrection_HomeMade(pwd, hemoCorr);
            elseif length(fList) == 1
                NewFluo1 = HemoCorrection_HomeMade(pwd, hemoCorr);
            else
                disp('No fluo channel detected, hemodynamic correction skipped');
            end
            for w = 1:size(fList)
                nbColors = nbColors +1;
                fname = strcat('corr_', fList(w).name);
                if contains(fname, 'fluo_475')
                    mockIllum = struct('ID',nbColors,'Color','Fluo #1 475 nmCorrected');
                elseif contains(fname, 'fluo_567')
                    mockIllum = struct('ID',nbColors,'Color','Fluo #2 567 nmCorrected');
                end
                fid = fopen(fname,'w');
                fwrite(fid, eval(strcat('NewFluo', int2str(w))),'single');
                eval(['AcqInfoStream.Illumination' int2str(nbColors) '= mockIllum ;']);
            end
            clear fList fname fid mockIllum NewFluo1 NewFluo2
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Normalization for  each channel    %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        count = count + 1;
        IndivSignalsFilt(count).TrialType = Trial(i).Type;
        
        for j = 1:nbColors %new loop because new channels created
            %%% Select first channel %%%
            chanID = eval(['AcqInfoStream.Illumination' int2str(j) ';']);
            chanList{1,j} = chanID.Color;
            [chan_name, info_name,structField] = selectChannelFiles(chanID.Color);         

            %%% Open metadata %%%
            infos = matfile(info_name);
            chanList{2,j} = infos.Freq;  

            %%% Open channel %%%
            fid = fopen(chan_name,'r+');
            data = fread(fid, Inf, 'single=>double');
            data = reshape(data,infos.datSize(1,1),infos.datSize(1,2),infos.datLength);

            %%% dI/I0 %%%
            if i == 1 && j == 1
                [dataFilt, data] = normalize_signal(data,infos,Trial(i),Parameters,[],0);
                nFrames = size(dataFilt,3); %number of frames on first trial, to pass to subsequent channels/trials
            else
                [dataFilt, data] = normalize_signal(data,infos,Trial(i),Parameters,nFrames,0);
            end

            %%% Store in structure %%%
            IndivSignalsFilt(count).(structField) = dataFilt;
            clear dataFilt

            %%% Save normalized data for Hemocompute
            if Parameters.ComputeHbO == 1
                if contains(chan_name,{'green.dat','red.dat','yellow.dat','Rx_567.dat'})
                    metaData = load(info_name);
                    metaData.datLength = nFrames;
                    mkdir normData;
                    save2Dat(fullfile(pwd, 'normData', chan_name), single(data+1), metaData);
                    copyfile(fullfile(pwd,'AcqInfos.mat'),fullfile(pwd, 'normData','AcqInfos.mat'));
                end
            end
            clear data
        end
        clear info_name chan_name

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Compute HbO and HbR if necessary   %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Parameters.ComputeHbO == 1
            folder = fullfile(pwd, 'normData');
            datFile = dir(fullfile(folder,'*.dat'));
            Illumination = lower(extractBefore({datFile.name},'.dat'));
            idx = or(or(or(contains(Illumination,'green'),contains(Illumination,'red')),contains(Illumination,'yellow')),contains(Illumination,'rx_567'));
            disp("Starting HbO and HbR computation...");
            [HbO, HbR] = HemoCompute_HomeMade(folder, pwd, Parameters.FilterSet, Illumination(idx),0);
            disp("HbO and HbR computation done!")

            IndivSignalsFilt(count).HbO = normalize_signal(HbO,infos,Trial(i),Parameters,nFrames,1);
            IndivSignalsFilt(count).HbR = normalize_signal(HbR,infos,Trial(i),Parameters,nFrames,1);
            clear datFile HbO HbR
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Delete info and dat files - save space%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fclose('all');
        fList = [dir([strcat(pwd, filesep) '*.mat']); dir([strcat(pwd, filesep) '*.dat'])];
        for w = 1:length(fList)
            if saveRaw == 0
                delete(fList(w).name);
            end
        end
        if exist('normData','dir')==7
            rmdir normData s
        end
        clear fList
        cd(root);
    else %discard trial if movement and user wants to exclude them
        fprintf('Discard Trial %s \n',D(i).name)
        fclose('all');
        fList = [dir([strcat(pwd, filesep) '*.mat']); dir([strcat(pwd, filesep) '*.dat'])];
        for w = 1:length(fList)
            delete(fList(w).name);
        end
        cd(root)
        movefile(fullfile(root,D(i).name),fullfile(root,'Discarded',D(i).name))
    end
end %%END OF COMPILATION OF TRIALS

%% Split, average, bin and save channels %%
SignalsFilt = struct('Green', cell(2,1), 'Amber',cell(2,1), 'Red', cell(2,1), 'HbO', cell(2,1),'HbR', cell(2,1), 'F1', cell(2,1), 'F1corr', cell(2,1), 'F2', cell(2,1), 'F2corr', cell(2,1));
chanIDX = ~structfun(@isempty, IndivSignalsFilt(1));
newChanList = fieldnames(IndivSignalsFilt);
newChanList = newChanList(and(chanIDX,~contains(newChanList,'TrialType')));
newChanList(:,2) = chanList(2,1);
nbColors = length(newChanList); %check non-empty fields in structure, -1 to count for TrialType

count = 1;
for i = 1:nbColors % loop on each color
    chanID = newChanList{i,1};
    gFiltSig = cat(4,IndivSignalsFilt(:).(chanID));
    [Av_gfiltSig] = AverageTrials(gFiltSig,newChanList{i,2},Parameters.final_tBin, strcat(chanID,'_gFilter'), contains([IndivSignalsFilt.TrialType],'Real'));
    SignalsFilt(1).(chanID) = Av_gfiltSig(:,:,:,1);
    if size(Av_gfiltSig,4)== 2
        SignalsFilt(2).(chanID) = Av_gfiltSig(:,:,:,2);
    end
    
    if i == 1      
        fig = figure("Position",[1 41 1920 963]);
        %Count number of frames not discarded per timepoint
        real_nan = ~isnan(squeeze(mean(gFiltSig(:,:,:,contains([IndivSignalsFilt.TrialType],'Real')),[1,2],"omitnan"))); %1 if frame for each trial is ok, 0 if all pixels are NaN on that frame
        sham_nan = ~isnan(squeeze(mean(gFiltSig(:,:,:,~contains([IndivSignalsFilt.TrialType],'Real')),[1,2],"omitnan")));
        nRealFrames = sum(real_nan,2,"omitnan"); %Count number of non-NaN frames
        nShamFrames = sum(sham_nan,2,"omitnan");
       
        %Plot graph
        subplot(2,3,1);
        plot(nRealFrames);
        hold on
        plot(nShamFrames);
        hold off
        title("Frames per timepoint")
        ylim([0, max(nRealFrames)+5]);
        legend("Real","Sham");
        ylabel("Number of averaged frames");
        xlabel("Frame number")        
        clear real_nan sham_nan
    end
    if ~contains(chanID,'Hb') 
        count = count+1;
        %Plot mean+sem value for each channel except HbO and HbR
        mReal = squeeze(mean(gFiltSig(:,:,:,contains([IndivSignalsFilt.TrialType],'Real')),[1,2,4],"omitnan")); %Mean pixel value for real
        mSham = squeeze(mean(gFiltSig(:,:,:,~contains([IndivSignalsFilt.TrialType],'Real')),[1,2,4],"omitnan"));%Mean pixel value for sham
        SEM_Real = squeeze(std(gFiltSig(:,:,:,contains([IndivSignalsFilt.TrialType],'Real')),1,[1,2,4],"omitnan"))./nRealFrames; %stdev for real
        SEM_Sham = squeeze(std(gFiltSig(:,:,:,~contains([IndivSignalsFilt.TrialType],'Real')),1,[1,2,4],"omitnan"))./nShamFrames; %Stdev for sham

        subplot(2,3,count)
        x1 = 1:length(mReal);
        plot(x1,mReal,"Color",[0 0.4470 0.7410]) %mean real
        hold on
        plot(x1,mSham,"Color",[0.8500 0.3250 0.0980]) %mean sham        

        %plot error bars
        x2 = [x1, fliplr(x1)]';
        inBetween = [(mReal - SEM_Real); flipud(mReal + SEM_Real)];
        x2(isnan(inBetween)) = [];
        inBetween(isnan(inBetween)) = [];
        fill(x2, inBetween, [0 0.4470 0.7410],'FaceAlpha',0.25,'EdgeColor','none'); %real SEM        
        inBetween = [(mSham + SEM_Sham); flipud(mSham - SEM_Sham)];
        x2 = [x1, fliplr(x1)]';
        x2(isnan(inBetween)) = [];
        inBetween(isnan(inBetween)) = [];
        fill(x2, inBetween, [0.8500 0.3250 0.0980],'FaceAlpha',0.25,'EdgeColor','none'); %sham SEM
        legend("Real","Sham","Real - SEM", "Sham - SEM");

        title(chanID)
        ylabel("Average pixel value (% dI/I0)");
        xlabel("Frame number")         
        hold off 
    end
    savefig(fig,fullfile(root,"TrialAverageFrameNumbers.fig"));
    clear gFiltSig Av_gfiltSig inBetween x2 mReal mSham SEM_Real SEM_Sham
end

% Calculate HbT
if sum(contains(newChanList(:,1),'HbO'))
    HbO = cat(4,IndivSignalsFilt(:).HbO);
    HbR = cat(4,IndivSignalsFilt(:).HbR);
    HbT = HbO + HbR;
    [~] = AverageTrials(HbT,newChanList{find(contains(newChanList(:,1),'HbO')),2},Parameters.final_tBin, 'HbT_gFilter', contains([IndivSignalsFilt.TrialType],'Real'));
    clear HbO HbR HbT
end
% Extract and average whisker video
if isempty(Trial(1).vid) == 0
    vid = extractVid(Trial, contains([IndivSignalsFilt.TrialType],'Real'));
    saveTifStack(vid, "WhiskVid");
end

%% Save individual trials if necessary
if indivTrial == 1
    save("IndivProcessedSignals",'IndivSignalsFilt','-v7.3');
elseif indivTrial == 0 && exist(fullfile(pwd,"IndivProcessedSignals.mat"),'file')==2
    delete(fullfile(pwd,"IndivProcessedSignals.mat"));
end
if flag_Rx567
    writelines('Warning: Green is Lime REFLECTANCE!','Warning.txt')
end
FPS = min([chanList{2,:}]);
if isequal(FPS,max([chanList{2,:}])) == 0
    writelines('Warning: FPS is not equal between channels!','Error.txt') %% shouldn't be the case, but just to be sure
end
end