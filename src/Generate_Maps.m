% Maps Mars' magnetic field from MAVEN data
% @author Brenden Lech

clearvars

% Reads in MAVEN data from a csv file
data = readtable("MAVEN 2017-2019.csv");
"Data read from file"
rMars = 3389.5; % Mars radius (km)
% Filters out erroneous B field measurements and data where MAVEN is >0.3
% Mars radii from the surface
filter = sqrt(data.SPACECRAFT_X_GEO____1___km.^2 + ...
    data.SPACECRAFT_Y_GEO____2___km.^2 + ...
    data.SPACECRAFT_Z_GEO____3___km.^2) < rMars * 1.3 & ...
    abs(data.MAG_FIELD_GEOX____1___nT) < 1e30 & ...
    abs(data.MAG_FIELD_GEOY____2___nT) < 1e30 & ...
    abs(data.MAG_FIELD_GEOZ____3___nT) < 1e30;
data = data(filter,:);
pos = data{:,7:9}; % MAVEN position in GEO coordinates (km)
b = data{:,4:6}; % B field measurements in GEO cooridnates (nT)
longitudes = data{:,2}; % MAVEN longitude
latitudes = data{:,3}; % MAVEN latitude
numData = size(longitudes,1); % Number of data points in the dataset
resolution = 400; % The resolution of the map produced
"Data filtered"

% Finds the maximum measured magnetic field strength
bMagnitudes = zeros(numData, 1);
for i = 1:numData
    bMagnitudes(i) = norm(b(i,:));
end
bMax = max(bMagnitudes);

"Found maximum B field strength"

% Counts the number of data points in each pixel of the map
numDataInPixels = zeros(resolution);
dataMapCoordinates = zeros(numData, 2); % The x,y coordinate of each data point in map pixels
for i = 1:numData
    
    longitude = longitudes(i);
    x = 1 + floor(longitude / 360 * resolution);
    x = min(x, resolution);
    
    % surfacem seems to put matrix pixel x=0 at longitude -180 rather than
    % longitude 0. Correcting for that
    x = floor(x - resolution / 2);
    if x < 1
        x = x + resolution;
    end
    
    latitude = 90 + latitudes(i);
    y = 1 + floor(latitude / 180 * resolution);
    y = min(y, resolution);
    
    dataMapCoordinates(i, 1) = x;
    dataMapCoordinates(i, 2) = y;
    numDataInPixels(y, x) = numDataInPixels(y, x) + 1;
    
end

"Counted data points per map pixel"

% Prepares the map matrix containing magnetic field magnitude for each map
% pixel, averaging all measurements in each pixel
magnitudeMapMatrix = zeros(resolution);
for i = 1:numData
    
    x = dataMapCoordinates(i, 1);
    y = dataMapCoordinates(i, 2);
    magnitudeMapMatrix(y, x) = magnitudeMapMatrix(y, x) + ...
        bMagnitudes(i) / numDataInPixels(y, x);
    
end

"Prepared magnitude map matrix"

% Prepares the map matrix containing magnetic field direction for each map
% pixel, averaging all measurements in each pixel
directionMapMatrix = zeros(resolution);
for i = 1:numData
    
    sphereNormalUnit = pos(i,:) / norm(pos(i,:));
    bUnit = b(i,:) / bMagnitudes(i);
    anglecos = dot(sphereNormalUnit, bUnit);
    
    x = dataMapCoordinates(i, 1);
    y = dataMapCoordinates(i, 2);
    directionMapMatrix(y, x) = directionMapMatrix(y, x) + ...
        anglecos / numDataInPixels(y, x);
    
end

"Prepared direction map matrix"

% Prepares the map matrix containing altitude for each map pixel,
% averaging all measurements in each pixel
altitudeMapMatrix = zeros(resolution);
for i = 1:numData
    
    altitude = sqrt(data.SPACECRAFT_X_GEO____1___km(i).^2 + ...
    data.SPACECRAFT_Y_GEO____2___km(i).^2 + ...
    data.SPACECRAFT_Z_GEO____3___km(i).^2) - rMars;
    
    x = dataMapCoordinates(i, 1);
    y = dataMapCoordinates(i, 2);
    altitudeMapMatrix(y, x) = altitudeMapMatrix(y, x) + ...
        altitude / numDataInPixels(y, x);
    
end

"Prepared altitude map matrix"

% Plots MAVEN's magnetic field direction measurements using a robinson
% projection
figure;
worldmap("World");
mlabel("south");
title("Mars Magnetic Field Polarity", ...
    "1=North, -1=South (Robinson Projection)")
surfacem([-90, 90], [-180, 180], directionMapMatrix);
colorResolution = 2065; % Must be odd
mapColors = turbo(colorResolution);
mapColors(ceil(colorResolution / 2), :) = [1, 1, 1];
colormap(mapColors);
clim([-1, 1]);
colorbar;

% Plots MAVEN's magnetic field direction measurements using an equidistant
% cylindrical projection
figure;
eqdcylinDirectionMap = worldmap("World");
setm(eqdcylinDirectionMap, "mapprojection", "eqdcylin");
setm(eqdcylinDirectionMap, "Origin", [0, 180]);
mlabel("south");
title("Mars Magnetic Field Polarity", ...
    "1=North, -1=South (Equidistant Cylindrical Projection)")
surfacem([-90, 90], [-180, 180], directionMapMatrix);
colorResolution = 2065; % Must be odd
mapColors = turbo(colorResolution);
mapColors(ceil(colorResolution / 2), :) = [1, 1, 1];
colormap(mapColors);
clim([-1, 1]);
colorbar;

% Plots MAVEN's magnetic field magnitude measurements
figure;
worldmap("World");
mlabel("south");
title("Mars Magnetic Field Magnitude (nT)", ...
    "(Robinson Projection)")
surfacem([-90, 90], [-180, 180], magnitudeMapMatrix);
colorResolution = 2065;
mapColors = jet(colorResolution);
mapColors(1, :) = [0, 0, 0];
colormap(mapColors);
colorbar;

% Plots MAVEN's average altitudes
figure;
worldmap("World");
mlabel("south");
title("MAVEN Measurement Altitude (km)", ...
    "With respect to Mars radius of 3389.5 km (Robinson Projection)")
surfacem([-90, 90], [-180, 180], altitudeMapMatrix);
colorResolution = 2065;
mapColors = gray(colorResolution);
mapColors(1, :) = [1, 1, 1];
colormap(mapColors);
colorbar;