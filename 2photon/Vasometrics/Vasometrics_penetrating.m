function [diam, fwhms, CLs, fig] = Vasometrics_penetrating(stack,metadata, CLs, minDist, cLength,alpha,handle)
% Function to analyse penetrating vessel segment diameter change overtime
%
% Inputs : 
%   - stack : time series with vessel-labeling channel in XYT format stack must be motion corrected
%   - metadata : stackmetadata
%   - CLs : cell containing x and y coordinates of crossline endpoints. Use
%   {} if this is the first run on this segment
%   -  minDist: 0 = take the maximal fwhm, 1= take minimal fwhm
%   - cLength : crossline length to add on each side of the circle
%   diameter.
%   - alpha : angle between each cross-lines in degrees
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
% Adapted to penetrating arteriole by Antoine Malescot


%% Get stack infos %%
dims = size(stack);
PixelSize = metadata.PixelSize;

if isempty(handle)==0 %calculate resize factor if figure already exists
    rFactor = size(handle.Children.CData)/dims(1:2);
end

%% Extract Rosas position
if isempty(CLs)
    mip = max(stack,[],3);
    if isempty(handle)
        imshow(mat2gray(imresize(mip,3)));
        L = drawcircle('Color', 'g'); % Draw circle
        Radius = round(L.Radius/3);
        Center = round(L.Center/3);

        % Obtain coordinates for circle
        mask = createMask(L);
        boundaries = bwboundaries(mask);
        x = round(boundaries{1}(:,2)./3); % x coordinates
        y = round(boundaries{1}(:,1)./3); % y coordinates
    else
        L = drawcircle('Color', 'g','Parent',handle); % Draw circle
        Radius = round(L.Radius/rFactor);
        Center = round(L.Center/rFactor);

        % Obtain coordinates for circle
        mask = createMask(L);
        boundaries = bwboundaries(mask);
        x = round(boundaries{1}(:,2)./rFactor); % x coordinates
        y = round(boundaries{1}(:,1)./rFactor); % y coordinates
    end   

    % Create several small slopes (sampling)
    slope_qty = round(180/alpha); % Number of tangantes

    Position(1,:) = [x(1)-cLength y(1)];
    Position(2,:) = [x(round(length(x)/2))+cLength y(round(length(x)/2))];
    delete(L)

end

%% Obtain coordinates for each line segment
if isempty(CLs)

    x = [];
    y = [];
    for z = 1:length(Position)-1
        [cx, cy, ~] = improfile(mip,[Position(z,1) Position(z+1,1)],[Position(z,2) Position(z+1,2)]);
        x = [x; cx];
        y = [y; cy];
    end

    clear cx cy
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

        xy_start = [x(1),y(1)]; % Coordinates of  vector's start 
        xy_end = [x(end),y(end)]; % Coordinates of vector's end
        
        for i = 0:slope_qty-1 % For loop on the number of angles to do

            R = [cosd(i*alpha) -sind(i*alpha); sind(i*alpha) cosd(i*alpha)]; % Matrix for rotation

            % Slide the vector to the origin 
            xy_start_origin = xy_start - Center; % Vector's start
            xy_end_origin = xy_end - Center; % Vector's end

            % Apply Rotation
            XY_0_rot = R*xy_start_origin'; % Start
            XY_end_rot = R*xy_end_origin'; % End

            % Slide the vector to the circle's center
            XY_0_rot = XY_0_rot  + Center; % Start
            XY_end_rot = XY_end_rot  + Center; % End
    
            % Extract coordinates
            movingX = [XY_0_rot(1) XY_end_rot(1)]; % x axis
            movingY = [XY_0_rot(end) XY_end_rot(end)]; %y axis

            % If the vector rotates above 90deg, then the coordinates are
            % reversed
            if (i*alpha > 90)
                movingX = [XY_end_rot(1) XY_0_rot(1)];
                movingY = [XY_end_rot(end) XY_0_rot(end)];
            end

            % Slope evalue the line position
            slope = (movingY(1) - movingY(end)) / (movingX(1) - movingX(end));
            invSlope = -1/slope;

            % Draw crossline
            endptX = [movingX(1), movingX(end)];
            endptY = [movingY(1), movingY(end)];

            % Evaluate the cx, cy values 
            [cx,cy,c] = improfile(stack(:,:,t),endptX,endptY);
            plot(cx, cy, '-r');
            hold on
            CLs(idx, :) = {endptX, endptY};
            c = (c-min(c)) / (max(c) - min (c));
            C{idx} = c;

            % Move the current x value
            lastX = movingX;
            lastY = movingY;
            idx = idx + 1;

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
% stdev = 1

fig = gcf;
end

