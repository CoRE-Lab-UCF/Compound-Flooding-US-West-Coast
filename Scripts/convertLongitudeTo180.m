function lon180 = convertLongitudeTo180(lon)
    % This function converts longitude from the 0 to 360 degree system
    % to the -180 to +180 degree system.

    % Input:
    % lon - longitude in the 0 to 360 degree system

    % Output:
    % lon180 - longitude in the -180 to +180 degree system

    % Wrap longitudes greater than 180 degrees
    lon180 = lon;
    lon180(lon > 180) = lon(lon > 180) - 360;
end
