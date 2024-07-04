function fwhms = getFWHMFromLineProfiles(C, PixelSize, minDist)
% Vasometric function get cross line length
% Inputs :
%   - C : array of crossline profiles. First dimension is the line profile, second
%   dimension is the various lines
%   - PixelSize : selft explanatory
%   - minDist: 0 = take the maximal fwhm, 1= take minimal fwhm
%
% Outputs : 
%   - fwhms : calculated fwhms a 1x crossline vector

%% Script %%
fwhms = zeros(size(C,2),1);
for i = 1 : size(C,2)
    profile = C{i};
    
    % Obtain the FWHM value for this profile
    halfMax = max(profile)/2;
    intersects = getYIntersects(halfMax,profile);        
    fwhm = (max(intersects) - min(intersects))*PixelSize;
    
    % Determine the derivative of this profile and use it to expand the bounds of FWHM height
    derivative = diff(profile);   
    medianV = median(profile);
    intersectsNew = getYIntersects(0,profile);
    leftIntChange = 1;
    rightIntChange = length(profile);
    
    for j = 1:length(intersectsNew)
        if intersectsNew(j) > leftIntChange && min(intersects) > intersectsNew(j)
            leftIntChange = intersectsNew(j);
        end
        if intersectsNew(j) < rightIntChange && intersectsNew(j) > max(intersects)
            rightIntChange = intersectsNew(j);
        end
    end
    
    % Using the adjusted intersects, recalculate the half max and find the x distance
    halfMax = (max(profile) - min(profile(round(leftIntChange):round(rightIntChange))))/2;
    fwhms(i) = fwhmFromProfile(profile, halfMax, round((rightIntChange + leftIntChange)/2),minDist) * PixelSize;
end

    
    
       
    

