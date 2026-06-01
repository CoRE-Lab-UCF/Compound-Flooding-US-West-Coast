function lon360 = convertLongitudeTo360(lon)
    % This function converts longitude from the -180 to +180 degree system 
    % to the 0 to 360 degree system.

    % Input:
    % lon - longitude in the -180 to +180 degree system

    % Output:
    % lon360 - longitude in the 0 to 360 degree system
    
    % Convert negative longitudes to positive equivalents
    lon360 = mod(lon, 360);
    
    % Ensure that longitudes that were originally negative are adjusted properly
    lon360(lon360 < 0) = lon360(lon360 < 0) + 360;
end
