function [maxPeakDist] = getCLineLength(x,y,position, mip)
% Vasometric function get cross line length
% Inputs :
%   - x : coordinates of line in x
%   - y : coordinates of line in y
%   - mip : maximum intensity projection of stack
%
% Outputs : 
%   - cLength : crossline length

excl = 10; %pixels to exclude from the edge of the line

%% Script
mip = medfilt2(mip, [10, 10]);
dims = size(mip);

maxPeakDist = 0;
maxPeakDist2 = []; % for diagnotic purposes only
imshow(mat2gray(mip));
hold on
plot(x,y,'-g');
idx = 1;
C = {}; %for diagnostic purposes only
peakDists = []; %for diagnostic purposes only

for i = 1:size(position,1)-1
    invM = -1/((position(i+1,2) - position(i,2)) / (position(i+1,1) - position(i,1)));    
    maxX = position(i+1,1);
    maxY = position(i+1,2);       

    while x(idx,1) ~= maxX || y(idx,1) ~= maxY
        invB = ((y(idx+1,1) + y(idx,1)) / 2) + (-1 * invM * ((x(idx+1,1) + x(idx,1)) / 2));     
        if invM == -Inf
           [cx,cy,c] = improfile(mip,[((x(idx+1,1) + x(idx,1)) / 2),((x(idx+1,1) + x(idx,1)) / 2)],[0,dims(1)]);
        else            
           [cx,cy,c] = improfile(mip,[round((0-invB)/invM),round((dims(1)-invB)/invM)],[0,dims(1)]);
        end
        imLim = find(cx > 0 & cx < dims(2));
        c = c(imLim);
        cx = cx(imLim);
        cy = cy(imLim);
        
        c = c(excl:length(c)-excl);
        cx = cx(excl:length(cx)-excl);
        cy = cy(excl:length(cy)-excl);
        plot(cx,cy,'-r');
        c = (c-min(c)) / (max(c) - min (c));
        
        halfMax = max(c)/2;
        dist = find(c == max(c),1);  %approximate vessel center
        fwhm = fwhmFromProfile(c, halfMax, dist,0);
        idx1 = find(c >= halfMax, 1, 'first');
        idx2 = find(c >= halfMax, 1, 'last');
        
        fwhm_x = cx(idx2)-cx(idx1);
        fwhm_y = cy(idx2)-cy(idx1);
        fwhm = sqrt(fwhm_x^2 + fwhm_y^2);
        
        peakDist = fwhm + 0.65*fwhm;
        
        maxPeakDist = max([maxPeakDist, peakDist]);
        maxPeakDist2 = [maxPeakDist2, peakDist];
        idx = idx + 1;
    end 
end


maxPeakDist = maxPeakDist/2;
    
