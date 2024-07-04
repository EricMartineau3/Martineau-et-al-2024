function [chan_name, info_name, structField] = selectChannelFiles(Color)

switch Color
    case 'Green'
        chan_name = 'green.dat';
        info_name = 'green.mat';
        structField = Color;
    case 'Amber'
        chan_name = 'yellow.dat';
        info_name = 'yellow.mat';
        structField = Color;
    case 'Red'
        chan_name = 'red.dat';
        info_name = 'red.mat';
        structField = Color;
    case 'Fluo #1 475 nm'
        chan_name = 'fluo_475.dat';
        info_name = 'fluo_475.mat';
        structField = 'F1';
    case 'Fluo #2 567 nm'
        chan_name = 'fluo_567.dat';
        info_name = 'fluo_567.mat';
        structField = 'F2';
    case 'Fluo #1 475 nmCorrected'
        chan_name = 'corr_fluo_475.dat';
        info_name = 'fluo_475.mat';
        structField = 'F1corr';
    case 'Fluo #2 567 nmCorrected'
        chan_name = 'corr_fluo_567.dat';
        info_name = 'fluo_567.mat';
        structField = 'F2corr';
    case 'Rx 567 nm'
        chan_name = 'Rx_567.dat';
        info_name = 'Rx_567.mat';
        structField = 'Green';
end