% Author: Pravin
% Description: Stratifies compound POT events into AR and non-AR events using the AR catalog.

clear
clc

% -------------------------------------------------------------------------
% Load input data
% -------------------------------------------------------------------------
base_path = fullfile('data');

load(fullfile(base_path,'Creating_WLONC','Time_series_data.mat'));
load(fullfile(base_path,'AR_data','AR_data_California.mat'));
load(fullfile(base_path,'Two_way_sampling','POT_both_extreme.mat'));

ar_catalog = fullfile(base_path,'AR_data','globalARcatalog_ERA5_1940-2024.nc');

lon = ncread(ar_catalog,'lon');
lat = ncread(ar_catalog,'lat');

% Location of catchment centroids
latt_gauge = Centroids_catchment(:,1);
long_gauge = Centroids_catchment(:,2);

% -------------------------------------------------------------------------
% Read AR catalog time
% -------------------------------------------------------------------------
T_num = ncread(ar_catalog,'time');
T1 = datenum(1900,1,1,0,0,0);
Time_AR = T_num./24 + T1;

AR_data = {};
non_AR_data = {};

% -------------------------------------------------------------------------
% Stratify POT events into AR and non-AR events
% -------------------------------------------------------------------------
for jj = 1:length(Centroids_catchment)

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

    AR_data{jj,1} = AR;
    non_AR_data{jj,1} = non_AR;

end

% -------------------------------------------------------------------------
% Save stratified event data
% -------------------------------------------------------------------------
save(fullfile(base_path,'Stratification','AR_data.mat'),'AR_data','-v6');
save(fullfile(base_path,'Stratification','non_AR_data.mat'),'non_AR_data','-v6');
