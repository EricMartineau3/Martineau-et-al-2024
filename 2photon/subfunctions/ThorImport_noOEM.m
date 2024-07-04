function [stack, metadata, xmlMeta] = ThorImport_noOEM(path)

% Imports image file and metadata from a ThorImage experiment and create a
% multidimentional stack out of it.
%
% References scripts :
% - Wouter Falkena (2021). xml2struct (https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct), MATLAB Central File Exchange. Retrieved May 7, 2021.
% - https://docs.openmicroscopy.org/bio-formats/5.7.1/developers/matlab-dev.html
%
% Used when data is saved in raw tiff and not OEM-TIFF
%
% Written by Eric Martineau - Universite de Montreal 

%% Import tiff image
cd(path);

% Filter file names
D = dir;
dirflag = [D.isdir];
D = D(~dirflag);

flag = contains({D.name},"Preview");
D = D(~flag);

flag = contains({D.name},"Chan");
D = D(flag);

k = length(D);
fID = fullfile(path,D(1).name);
im = imread(fID);
x = size(im,1);
y = size(im,2);

clear flag dirflag fID im

% Loop for import
data = zeros(x,y,k);

for i = 1:length(D)        
    fID = fullfile(path,D(i).name);
    data(:,:,i) = imread(fID);
end

clear dirflag flag fID D
%%  Extract metadata
xmlMeta = xml2struct(fullfile(path,'Experiment.xml'));
xmlMeta = xmlMeta.ThorImageExperiment;

stackSizeX = str2double(xmlMeta.LSM.Attributes.pixelX); % image width, pixels
stackSizeY = str2double(xmlMeta.LSM.Attributes.pixelY); % image height, pixels
stackSizeZ = str2double(xmlMeta.ZStage.Attributes.steps); % number of Z slices
stackSizeT = str2double(xmlMeta.Timelapse.Attributes.timepoints); % number of T slices

PMTs = {'A', str2double(xmlMeta.PMT.Attributes.enableA);'B', str2double(xmlMeta.PMT.Attributes.enableB); 'C', str2double(xmlMeta.PMT.Attributes.enableC)}; 
%NbChan = sum(cell2mat(PMTs(:,2))); %bug in xmlMeta making this method
%unreliable
NbChan = size(xmlMeta.Wavelengths.Wavelength,2);

PixelSize = str2double(xmlMeta.LSM.Attributes.pixelSizeUM); %pixel size in micro-meters
FrameRate = str2double(xmlMeta.LSM.Attributes.frameRate); %Acquisition framerate
dwellTime = str2double(xmlMeta.LSM.Attributes.dwellTime); %Pixel dwell time
if str2double(xmlMeta.LightPath.Attributes.GalvoResonance) == 1  
    ScanMode = 'GR';
elseif str2double(xmlMeta.LightPath.Attrbutes.GalvoGalvo) == 1
    ScanMode = 'GG';
else
    ScanMode = 'NaN';
end
parts = strsplit(path, '\');
stackID = parts{end-1};

ExpNotes = xmlMeta.ExperimentNotes.Attributes.text;

%Store main metadata values needed for analysis etc. Orginial ThorImage xml metadata will be save with it. 
metadata = struct("StackID", stackID, "SizeX",stackSizeX, "SizeY",stackSizeY, "SizeZ", stackSizeZ, "SizeT",stackSizeT, "NbChannels", NbChan,"PixelSize", PixelSize, "FrameRate", FrameRate, "DwellTime", dwellTime, "EnabledPMTs", {PMTs}, "ScanMode", ScanMode, "ExpNotes", ExpNotes);
clear stackID stackSizeX stackSizeY stackSizeZ stackSizeT PMTs NbChan FrameRate PixelSize dwellTime ScanMode ExpNotes parts

%% Reorganize and cleanup data
dims = [metadata.SizeY metadata.SizeX metadata.SizeZ metadata.SizeT metadata.NbChannels];

logic = dims==1;
stackOrder = {'Y' 'X' 'Z' 'T' 'C'};
stackOrder = stackOrder(~logic);
clear logic

rawIDX = 1;

stack = zeros(dims,'uint16');
for C = 1:dims(1,5)
   for T = 1:dims(1,4)
       for Z = 1:dims(1,3)
           stack(:,:,Z,T,C)= data(:,:,rawIDX);
           rawIDX = rawIDX+1;
       end
   end
end

stack = squeeze(stack);
metadata.stackOrder = rot90(char(stackOrder));
