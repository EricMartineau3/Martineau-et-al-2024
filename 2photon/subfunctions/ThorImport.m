function [stack, metadata, xmlMeta] = ThorImport(path, file)

% Imports image file and metadata from a ThorImage experiment and create a
% multidimentional stack out of it.
%
% References scripts :
% - Wouter Falkena (2021). xml2struct (https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct), MATLAB Central File Exchange. Retrieved May 7, 2021.
% - https://docs.openmicroscopy.org/bio-formats/5.7.1/developers/matlab-dev.html
%
% Written by Eric Martineau - Universite de Montreal

%% Import using bioformat toolbox
fID = fullfile(path,file);
data = bfopen(fID);

%%  Extract metadata
omeMeta = data{1,4};
xmlMeta = xml2struct(fullfile(path,'Experiment.xml'));
xmlMeta = xmlMeta.ThorImageExperiment;

stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
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
           stack(:,:,Z,T,C)= data{1,1}{rawIDX,1};
           rawIDX = rawIDX+1;
       end
   end
end

stack = squeeze(stack);
metadata.stackOrder = rot90(char(stackOrder));
