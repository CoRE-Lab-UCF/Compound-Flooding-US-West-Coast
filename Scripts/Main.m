
clear
clc

% Stratification
load('H:\OneDrive - University of Central Florida\Analysis_AR_two_way_8_events_splitted_marg_comb_SS_Hmo\Creating_WLONC\Time_series_data.mat')
load('H:\OneDrive - University of Central Florida\Analysis_AR\AR_data\AR_data_California.mat');
lon=ncread('H:\OneDrive - University of Central Florida\Analysis_AR\AR_data\globalARcatalog_ERA5_1940-2024.nc','lon');
lat=ncread('H:\OneDrive - University of Central Florida\Analysis_AR\AR_data\globalARcatalog_ERA5_1940-2024.nc','lat');
load('H:\OneDrive - University of Central Florida\Analysis_AR_two_way_8_events_splitted_marg_comb_SS_Hmo\Two_way_sampling\POT_both_extreme.mat')

% location of catchment centroid
latt_gauge = Centroids_catchment(:,1);% Coresponds to the sites listed in RF accumulation data file
long_gauge = Centroids_catchment(:,2);


% % Creating the mesh grid
% [lonGrid,latGrid] = meshgrid(convertLongitudeTo180(lon),lat);  % Create meshgrid for lat/lon
% lat_flat = latGrid(:);  % Flatten latitude grid
% lon_flat = lonGrid(:);  % Flatten longitude grid


%% reading the time arra of the AR data
T_num=ncread("H:\OneDrive - University of Central Florida\Analysis_AR\AR_data\globalARcatalog_ERA5_1940-2024.nc","time");
T1=datenum(1900,1,1,0,0,0);
Time_AR=T_num./24+T1;



AR_data={};
non_AR_data={};



%%
for jj=1:length(Centroids_catchment)
    data = POT_both_extreme{jj,1};
    [AR, non_AR] = stratify_AR_RF(data, ...
        'Column_of_Time_RF', 1, ...
        'AR_Time',           Time_AR, ...
        'AR_lng',            lon, ...
        'AR_ltt',            lat, ...
        'time_bef_RF',       0.5, ...
        'time_aft_RF',       0.5, ...
        'latt_gauge',        latt_gauge(jj,1), ...
        'long_gauge',        long_gauge(jj,1), ...
        'dist_thres_RF',     200);

    AR_data{jj,1}=AR;
    non_AR_data{jj,1}=non_AR;
end

save('AR_data.mat','AR_data','-v6');
save('non_AR_data.mat','non_AR_data','-v6');

%% For AR
for j = 12
[lonGrid_plot, latGrid_plot] = meshgrid(convertLongitudeTo180(lon), lat);

% Find the NC file index matching this event's time
nc_idx = find(Time_AR>AR(j,1)-0.124 & Time_AR<AR(j,1)+0.124);

    
    for steps=1:4
        figure
        shape_plot = ncread('H:\OneDrive - University of Central Florida\Analysis_AR\AR_data\globalARcatalog_ERA5_1940-2024.nc', ...
            'shapemap', [1,1,1,nc_idx-3+steps,1], [inf,inf,1,1,1]);
        shape_plot = shape_plot';
    
        geoscatter(latGrid_plot(:), lonGrid_plot(:), 4, shape_plot(:), "filled", ...
            "MarkerEdgeColor", "none", "MarkerFaceColor", 'flat'); hold on
        geoscatter(latt_gauge(jj,1), long_gauge(jj,1), "filled", ...
            "MarkerEdgeColor", "none", "MarkerFaceColor", 'flat')
    end
end


%% For non-AR
for j = 2
[lonGrid_plot, latGrid_plot] = meshgrid(convertLongitudeTo180(lon), lat);

% Find the NC file index matching this event's time
nc_idx = find(Time_AR>non_AR(j,1)-0.124 & Time_AR<non_AR(j,1)+0.124); %to account for 6 hr inetrval

    
    for steps=1:4
        figure
        shape_plot = ncread('H:\OneDrive - University of Central Florida\Analysis_AR\AR_data\globalARcatalog_ERA5_1940-2024.nc', ...
            'shapemap', [1,1,1,nc_idx-3+steps,1], [inf,inf,1,1,1]);
        shape_plot = shape_plot';
    
        geoscatter(latGrid_plot(:), lonGrid_plot(:), 4, shape_plot(:), "filled", ...
            "MarkerEdgeColor", "none", "MarkerFaceColor", 'flat'); hold on
        geoscatter(latt_gauge(jj,1), long_gauge(jj,1), "filled", ...
            "MarkerEdgeColor", "none", "MarkerFaceColor", 'flat')
    end
end


