function intersects = getYIntersects(halfMax,profile)
%Vasometric function get x-axis value of halfmax, interpolating between pixels
% Inputs :
%   - halfMax : halfMax value
%   - profile : crossline profile
%
% Outputs : 
%   - intersects : x-values where Y = halfmax

%% Script 
intersects = [];
for i = 1:length(profile)-1
    profileSlope = profile(i+1) - profile(i);
    profileYInt = profile(i) + (-1 * profileSlope * i);
    xInt = ((halfMax - profileYInt) / profileSlope);
    
    if xInt >= i && xInt <= i+1
        intersects = [intersects xInt];
    end
end