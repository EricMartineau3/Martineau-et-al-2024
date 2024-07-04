function [Time_speed, Speed, par, X,max_freq_T, max_freq_X,xStall] = speed_fourier(im, freq, dx, dT_win, dT_step, X, max_freq_T, max_freq_X, nStdev,stall)

% Script written by Bruno-Félix Osmanski, Charpak lab, INSERM for Ravi Rungta
% Modified by Éric Martineau - Universite de Montreal to convert to function for integration in
% other analysis routines
%
% INPUTS
%
% im : stiched_linescan, must be T x X format (only 1 channel)
% path : self explanatory
% TrialID : identifier of trial for saving purposes
% freq : linerate of acquisition (in Hz)
% dx : pixel size in ummeters (um)
% dT_win : delta-timewindow, in seconds. Determines the timewindow
% overwhich the velocity is calculated
% dT_step : step by which the window is moved to calculate the velocity
% over time. Must be equal or inferior to dT_win. Increases runtime if
% too low.
% X, max_freq_X and max_freq_T : user-inputed variables feeded to function
% after first run. 


%% parameters
norm_on = 0; % zero = no normalisation

% median filter to be applied on image
medfilt_on = [1 1]; %[2 2]; % put zero to remove medfilt
Time_plot_select = [2 2.2]; %(s)
dx = dx * 1e-3; % convert to mm
time_1_line = 1/freq;

%% loading image  
im = single(im);
 
if norm_on == 1
    im = bsxfun(@minus,im,mean(im,1));
    im = bsxfun(@rdivide,im,sqrt(sum(abs(im).^2,1)));
end

%% Apply median filter on image
if medfilt_on(2) > 0 
    im = medfilt2(im,medfilt_on);
end
%%    
X_vect = [0:size(im,2)-1]*dx;
T_vect = [0:size(im,1)-1]/freq;

figure; 
imagesc([1 length(X_vect)],T_vect,  im)
ylabel('time (s)')
ylim(Time_plot_select)


%%
if isempty(X) == 1
    [X Y] = ginput(2);
    X = round(X);
else
    close gcf
end

%%
N_X_ex = [X(1):X(2)];
im = im(:,N_X_ex);
X_vect = [0:size(im,2)-1]*dx;
%%

win_N_T_anal = round(dT_win * freq );
step_N_T_anal = round(dT_step * freq );

N_T_anal = floor(  (size(im,1) - win_N_T_anal) / step_N_T_anal + 1 );
T_frame = [0:N_T_anal-1]/freq*step_N_T_anal; %Offset by dT_win

im_resh = zeros(win_N_T_anal,size(im,2),N_T_anal,'single');

for i = 1:N_T_anal
    im_resh(:,:,i) = im([1:win_N_T_anal] + (i-1)*step_N_T_anal ,: );
    
end

im_resh = bsxfun(@minus,im_resh,mean(mean(im_resh,1),2));


N_fft_T = 1024; %2 ^ ceil( log(size(im_resh,1)) / log(2) );
N_fft_X = 1024; %2 ^ ceil( log(size(im_resh,2)) / log(2) );
freq_X = [-N_fft_X/2:N_fft_X/2-1]/N_fft_X*(1/dx);
freq_T = [-N_fft_T/2:N_fft_T/2-1]/N_fft_T*freq;
%%
fft_im_resh = abs(fft2(im_resh,N_fft_T,N_fft_X));
    
tmp = zeros(size(fft_im_resh),'single');
tmp(1:N_fft_T/2,1:N_fft_X/2,:) = fft_im_resh(N_fft_T/2+1:end,N_fft_X/2+1:end,:);
tmp(N_fft_T/2+1:end,N_fft_X/2+1:end,:) = fft_im_resh(1:N_fft_T/2,1:N_fft_X/2,:);
tmp(1:N_fft_T/2,N_fft_X/2+1:end,:) = fft_im_resh(N_fft_T/2+1:end,1:N_fft_X/2,:);
tmp(N_fft_T/2+1:end,1:N_fft_X/2,:) = fft_im_resh(1:N_fft_T/2,N_fft_X/2+1:end,:);

fft_im_resh = tmp;

%% set fft frequency

N_fft_T_all = 2 ^ ceil( log(size(im,1)) / log(2) );
freq_T_all = [-N_fft_T_all/2:N_fft_T_all/2-1]/N_fft_T_all*freq;
im_fft_T_all = fftshift(mean(abs(fft(im,N_fft_T_all,1)),2));

N_fft_X_all = 512;
freq_X_all = [-N_fft_X_all/2:N_fft_X_all/2-1]/N_fft_X_all/dx;
im_fft_X_all = fftshift(mean(abs(fft(im,N_fft_X_all,2)),1));

figure;
imagesc(freq_X,freq_T,sum(abs(fft_im_resh),3))
xlabel('mm-1')
ylabel('Hz')

%%
if isempty(max_freq_T) == 1
    uiwait();
    prompt = {'max_freqT', 'max_freq_X'};
    dlgtitle = "Processing parameters";
    dims = [1 50];
    definputs = {'120', '100'};

    answer = inputdlg(prompt,dlgtitle, dims, definputs);

    max_freq_T = str2double(answer{1}); %(Hz) % zero for all frequency
    max_freq_X = str2double(answer{2}); %(mm-1) % zero for all frequency
end

%%
if max_freq_T>0
    [tmp,N_max_freq_T] = min(abs(max_freq_T-freq_T)); %(Hz)
    N_max_freq_T = N_max_freq_T - N_fft_T/2;
    freq_T_ex = freq_T(N_fft_T/2 + [-N_max_freq_T+1:N_max_freq_T]);
else
    N_max_freq_T = N_fft_T/2;
    freq_T_ex = freq_T;
end

if max_freq_X>0
    [tmp,N_max_freq_X] = min(abs(max_freq_X-freq_X)); %(mm-1)
    N_max_freq_X = N_max_freq_X - N_fft_X/2;
    freq_X_ex = freq_X(N_fft_X/2 + [-N_max_freq_X+1:N_max_freq_X]);
else
    N_max_freq_X = N_fft_X/2;
    freq_X_ex = freq_X;
end

fft_im_resh_ex = fft_im_resh( N_fft_T/2 + [-N_max_freq_T+1:N_max_freq_T] , N_fft_X/2 + [-N_max_freq_X+1:N_max_freq_X] ,:);
%%
clear ind_Max Max ft_Max_ex fx_Max fx_Max_ex

for i = 1:size(fft_im_resh_ex,3)
    [Max(:,i),ind_Max(:,i)] = max(fft_im_resh_ex(:,:,i),[],2);
    fx_Max_ex(:,i) = freq_X_ex(ind_Max(:,i));
end

ft_Max_ex = freq_T_ex;


%% moindre carré tamisé

tic
clear coef
for i = 1:size(fx_Max_ex,2)
    %disp(i/size(fx_Max_ex,2))
    tmp_y = -ft_Max_ex';
    tmp_x = fx_Max_ex(:,i);
    
    for j = 1:round(size(fx_Max_ex,1)/2)
               
        [coef(:,i)] = polyfit(tmp_x,tmp_y,1);
        fit = polyval(coef(:,i),tmp_x);
        [tmp,ind_tmp] = max(abs(fit - tmp_y));
        
        tmp_y(ind_tmp) = [];
        tmp_x(ind_tmp) = [];
                
    end
    
end
toc

%% display

figure;
plot(T_frame(1:size(coef,2)) + T_frame(1),sign(median(squeeze(coef(1,:))))*coef(1,:),'')
hold on
xlabel('time (s)') 
ylabel('speed (mm/s)')
title(['column ' num2str(N_X_ex(1)) ' to ' num2str(N_X_ex(end))])

figure;
plot(sign(median(squeeze(coef(1,:))))*coef(1,:),'')
%% save data %%
par.dx = dx; 
par.freq = freq; 
par.max_freq_T = max_freq_T; 
par.max_freq_X = max_freq_X; 
par.dT_win = dT_win; 
par.dT_step = dT_step; 
par.N_X_ex = N_X_ex;
par.norm_on = norm_on;
par.medfilt_on = medfilt_on;

% Time_speed = T_frame(1:size(coef,2)) + T_frame(1); %original time calc
Time_speed = T_frame(1:size(coef,2)) + win_N_T_anal/freq; %correction by eric

Speed = sign(median(squeeze(coef(1,:))))*coef(1,:);

%% stall detection
if contains(stall,"manual") == 1 %Manual method to select regions where there is a stall
    %Plot trace
    figure;
    plot(Time_speed,Speed);
    xlabel('Time (s)');
    ylabel('RBC velocity (mm/s');

    %Draw stall areas
    ok = "Yes";
    idx = 0;
    while contains(ok,"Yes")
        idx = idx + 1;
        roi = drawrectangle('Color','r');
        tmp = [roi.Position(1) roi.Position(1)+roi.Position(3)];
        xStall(idx,1) = find(Time_speed >= tmp(1),1);
        xStall(idx,2) = find(Time_speed >= tmp(2),1);
        Speed(xStall(idx,1):xStall(idx,2)) = NaN;

        ok = questdlg("Select another area?","Stall selection","Yes","No","No");
    end
elseif contains(stall,"auto") == 1 % Automatic method which uses the matlab findchangepts to detect stalls with user input
    ok = "No";
    while contains(ok,"No")
        %Plot trace
        close all
        figure;
        plot(Time_speed,Speed);
        xlabel('Time (s)');
        ylabel('RBC velocity (mm/s');

        %find change points
        answer = inputdlg("Number of stalls:","Input",[1,35],"1");
        numc = str2double(answer);
        findchangepts(Speed,'MaxNumChanges',numc*2,'Statistic','std');

        %Validate
        ok = questdlg("Stall detection ok?","Validation","Yes","No","Yes");
    end
    ipt = findchangepts(Speed,'MaxNumChanges',numc*2,'Statistic','std');
    
    if rem(length(ipt),2) > 0 %if not even number
        ipt = [ipt,length(Speed)]; %add endpoint ... could be a problem if stall is present at start but I assume the user would have stopped the recording...
    end
    xStall(:,1) = ipt(1:2:end);
    xStall(:,2) = ipt(2:2:end);

    for i = 1:size(xStall,1)
        Speed(xStall(i,1):xStall(i,2)) = NaN;
    end
elseif contains(stall,"none")
    xStall = [];
elseif contains(stall,"none") == 0
    disp("Input error for nStdev");
end

%Filter out outliers
if isnan(nStdev) == 0    
    Speed(Speed < 0) = NaN;
    m = mean(Speed,2);
    stdev = std(Speed,0,2);

    Speed(Speed > m + nStdev*stdev) = NaN;
    Speed(Speed < m - nStdev*stdev) = NaN;
end
























