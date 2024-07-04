function [diam, fwhms, CLs, fig] = Vasometrics(stack,metadata, CLs, CLspace, minDist, cLength,handle)
% Function to analyse vessel segment diameter change overtime
%
% Inputs : 
%   - stack : time series with vessel-labeling channel in XYT format stack, must be motion corrected
%   - metadata : stackmetadata
%   - CLs : cell containing x and y coordinates of crossline endpoints. Use
%   {} if this is the first run on this segment
%   - CLspace : space between crosslines in um
%   -  minDist: 0 = take the maximal fwhm, 1= take minimal fwhm
%   - cLength : crossline length. Use [] if you want to calculate
%   automatically. 
%   - handle : axes handle to trace on.
%
% Ouput : 
%   - diam : average diameter measurments over time
%   - fwhms : fwhm for each line, for each frame
%   - CLs : cell containing x and y coordinates of crossline endpoints.
%   Used to draw lines in subsequent runs
%   - fig : figure display of crosslines
%
% Based on Andy Shih's macro Vasometerics in ImageJ 
% McDowell K, et al (2021)VasoMetrics: unbiased spatiotemporal analysis of
% microvascular diameter in multi-photon imaging applications.
% https://github.com/mcdowellkonnor/ResearchMacros.git
%
% Translated to Matlab by Éric Martineau

%% Get stack infos %%
dims = size(stack);
PixelSize = metadata.PixelSize;

if isempty(handle)==0 %calculate resize factor if figure already exists
    rFactor = size(handle.Children.CData)/dims(1:2);
end

%% UI to draw polyline across segment of interest %%
if isempty(CLs)  
    mip = max(stack,[],3);
    if isempty(handle)
        imshow(mat2gray(imresize(mip,3)));
        L = drawpolyline('Color', 'g');
        Position = round(L.Position/3);
    else
        L = drawpolyline('Color', 'g','Parent',handle);
        Position = round(L.Position/rFactor);
    end

    %Obtain coordinates for each line segment
    x = [];
    y = [];
    for z = 1:length(Position)-1
        [cx, cy, ~] = improfile(mip,[Position(z,1) Position(z+1,1)],[Position(z,2) Position(z+1,2)]);
        x = [x; cx];
        y = [y; cy];
    end 
    delete(L);
end
clear cx cy

%% Determine crosslines length and placement
if isempty(CLs)
    if isempty(cLength)
        cLength = getCLineLength(x, y,Position, mip);
    end
    cSpace = CLspace / PixelSize;

    totalLength = 0;
    for i = 1:length(x)-1
        totalLength = totalLength + sqrt(power(x(i+1,1) - x(i,1),2) + power(y(i+1,1) - y(i,1),2));
    end
    excessLength = rem(totalLength, cSpace);
end

%% Draw crosslines and measure FWHM at each slice
for t = 1:dims(3)    
    if t == 1 && isempty(CLs)%Draw line on first slice
        lastX = 0;
        lastY = 0;
        idx = 1;
        C = [];
        imshow(mat2gray(stack(:,:,t)));
        hold on
        
        for i = 1:size(Position,1)-1
            slope = (Position(i+1,2) - Position(i,2)) / (Position(i+1,1) - Position(i,1));            
            invSlope = -1/slope;
            
            if slope ~= Inf || slope ~= -Inf                
                yInt = Position(i,2) + (-1 * slope * Position(i,1));
            
                % Move the starting point based on the excess so that the cross-lines are centered
                if i == 1
                    val = getIntersection(x(i,1), y(i,1), excessLength/2, slope, yInt);
                    if Position(i+1,1) > Position(i,1)
                        movingX = max(val);
                    else
                        movingX = min(val);
                    end
                % Move the starting points on line-segments following the first so that equal cross-line distance is maintained
                else
                    val = getIntersection(lastX, lastY, cSpace, slope, yInt);
                    if Position(i+1,1) > Position(i,1)
                        movingX = max(val);
                    else
                        movingX = min(val);
                    end
                end

                while movingX <= max([Position(i,1), Position(i+1,1)]) && movingX >= min([Position(i,1), Position(i+1,1)])
                    movingY = slope * movingX + yInt;
                    if invSlope ~= -Inf               
                        invYInt = movingY + (-1 * invSlope * movingX);
                        val = getIntersection(movingX, movingY, cLength, invSlope, invYInt); %get x-values of endpoints of crossline
                    end

                    % Draw crossline
                    if invSlope == -Inf
                        endptX = [movingX, movingX];
                        endptY = [movingY + cLength, movingY - cLength];                                                
                    else
                        endptX = [val(1), val(2)];
                        endptY = [val(1)*invSlope+invYInt, val(2)*invSlope+invYInt];                        
                    end
                    [cx,cy,c] = improfile(stack(:,:,t),endptX,endptY);
                    plot(cx, cy, '-r');
                    hold on
                    CLs(idx, :) = {endptX, endptY};
                    c = (c-min(c)) / (max(c) - min (c));
                    C{idx} = c;

                    % Move the current x calue
                    lastX = movingX;
                    lastY = movingY;
                    idx = idx + 1;

                    val = getIntersection(lastX, lastY, cSpace, slope, yInt);
                    if Position(i+1,1) > Position(i,1)
                        movingX = max(val);
                    else
                        movingX = min(val);
                    end
                end
            else
                disp("Line segment is completely vertical ... redraw line segment with a slope different than Infinity");
            end           
        end
    else
        C = [];
        for cl = 1:length(CLs)
            endptX = CLs{cl,1};
            endptY = CLs{cl,2};     
            c = improfile(stack(:,:,t),endptX,endptY);
            c = (c-min(c)) / (max(c) - min (c));
            C{cl}= c;
        end    
    end
    % Parse through the cross-lines and obtain the FWHM values
    fwhms(:,t) = getFWHMFromLineProfiles(C, PixelSize, minDist); 
end

diam = mean(fwhms,1,'omitnan'); % filter out NaN    

fig = gcf;
