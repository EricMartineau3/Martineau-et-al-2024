function [traces, fig] = pullROIsignal(sig, Masks, avBin, fps, SeriesID, ROI_ID)
% Function to extract the signal from a ROI 
%
%   INPUTS:
%       - sig: xy-nSeries-t matrix, where the signal from multiple series (i.e whisker stim)
%               can be extracted from the each ROI (barrel). 
%       - Masks: xy-nMasks matrix of masks extracted using BarrelMask or the standard ROI alignement
%       - avBin: runing average window size (in number of frames)
%       - fps: framerate (in Hz) to convert frame numbers to time (in seconds)
%       - SeriesID: {nSeries x 1} containing the string identifiers of the recording being analyzed
%       - ROI_ID: {nMasks x 1} containing the string identifiers of the recording the ROI has been pulled from
%
%   OUTPUTS:
%       - traces: nSeries-nMasks-t traces for each series and ROI combination
%
% Written by Eric Martineau and Antoine Malescot - Universite de Montreal

%% Extract ROI Signals

NbROI = size(Masks,3);
NbSeries = size(sig,3);
time = size(sig,4);

traces = zeros(NbSeries, NbROI, time);

for i = 1:NbSeries
    serie = squeeze(sig(:,:,i,:));
    for j = 1:NbROI
        Mask = logical(Masks(:,:,j));
        for t = 1:time            
            frame = serie(:,:,t);
            traces(i,j,t) = mean(frame(Mask),'omitnan'); %series-roi-t matrix
        end
    end    
end
clear serie Mask frame i j t

%% Sliding window average
traces = permute(traces, [3 2 1]); %t-roi-series
traces = movmean(traces, [(avBin-1) 0],'omitnan');
traces = permute(traces, [3 2 1]); %series-roi-t

%% Create figures
x = linspace(0,round(time/fps,1),time);
ymax = max(traces, [], 'all');
ymin = min(traces, [], 'all');

tiledlayout(NbSeries,1);

newcolors = [0 0 1
             0 1 0
             1 0 1
             1 0 0
             1 1 0
             0 1 1
             1 0.5 0.3
             0.3 1 1
             0.3 0.5 1];
colororder(newcolors);

for i = 1:NbSeries
    nexttile      
    for j = 1:NbROI
        y = squeeze(traces(i,j,:));
        plot(x,y,'-');
        if j == 1
            hold on
            title(SeriesID{i,1});  
        elseif j == NbROI
            ylim([ymin ymax]);
            legend(strcat('ROI ', ROI_ID));
            hold off
        end
    end
end

fig = gcf;
    
