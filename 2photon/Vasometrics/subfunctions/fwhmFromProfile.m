function fwhm = fwhmFromProfile(profile, halfMax, vesselCenterX, minDist)
%Vasometric function get real fwhm value (in pixel units) and filters value
%to keep only those closest to the vessel center. 

% Inputs :
%   - profile : crossline profile
%   - halfMax : halfMax value
%   - vesselCenterX : x-axis value corresponding to center of vessel
%   - minDist : 0 = take the maximal fwhm, 1= take minimal fwhm
%
% Outputs : 
%   - fwhm :width between values where Y = halfmax;
%% Script
intersects = getYIntersects(halfMax,profile);
if halfMax < 0.2
    fwhm = NaN;
elseif    length(intersects) < 2
    fwhm = NaN;
else
    leftX = intersects(1);
    rightX = intersects(length(intersects));

    for x = 1:length(intersects)
        if minDist == 0
            if intersects(x) < vesselCenterX && (vesselCenterX - intersects(x)) > (vesselCenterX - leftX)
                leftX = intersects(x);
            elseif intersects(x)> vesselCenterX && (intersects(x) - vesselCenterX) > (rightX - vesselCenterX) 
                rightX = intersects(x);
            end
        elseif minDist == 1
            if intersects(x) < vesselCenterX && (vesselCenterX - intersects(x)) < (vesselCenterX - leftX)
                leftX = intersects(x);
            elseif intersects(x)> vesselCenterX && (intersects(x) - vesselCenterX) < (rightX - vesselCenterX) 
                rightX = intersects(x);
            end
        end
    end
    fwhm = rightX - leftX;
end