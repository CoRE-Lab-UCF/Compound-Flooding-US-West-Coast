% Author: Pravin
% Description: Identifies compound and one-way POT events and estimates SED parameters for rainfall and water level extremes.

clear
clc

% Analysis settings
Rf_Acc = 24;
compound_window = 1.5;
dec_tim = 1.5;

base_path = fullfile('data');
load(fullfile(base_path,'Creating_WLONC','Time_series_data.mat'));

% Percentile thresholds used to define extreme events
RF_percentile = 90;
SS_percentile = 90;

POT_both_extreme = {};
POT_RF_only_extreme = {};
POT_NTTWL_only_extreme = {};

SED_RF = nan(numel(shp),2);
SED_WLAM = nan(numel(shp),2);
SED_SS = nan(numel(shp),2);

for loc = 1:numel(shp)

    Data = Time_series_data{loc,1};

    % Compute daily rainfall totals for threshold estimation
    dateDates = floor(Data(:,1));
    [uniqueDays,~,idx] = unique(dateDates);
    dailyRainfall = [uniqueDays, accumarray(idx, Data(:,2))];

    % Determine rainfall and water-level thresholds
    threshold_RF = prctile(dailyRainfall(:,2), RF_percentile);
    threshold_WL = prctile(Data(:,5) + Data(:,6), SS_percentile);
    threshold_SS = prctile(Data(:,5), SS_percentile);

    % Estimate SED parameters for rainfall and water-level series
    [SED_RF(loc,1),SED_RF(loc,2)] = find_parameters( ...
        [Data(:,1),Data(:,2)], threshold_RF);

    [SED_WLAM(loc,1),SED_WLAM(loc,2)] = find_parameters( ...
        [Data(:,1),Data(:,5)+Data(:,6)], threshold_WL);

    [SED_SS(loc,1),SED_SS(loc,2)] = find_parameters( ...
        [Data(:,1),Data(:,5)], threshold_SS);

    % Identify compound rainfall–water-level extremes
    out = find_POT_compound_RF_WLNC( ...
        Data, threshold_RF, threshold_WL, ...
        'RfAccHours', Rf_Acc, ...
        'compound_window', compound_window, ...
        'dec_tim', dec_tim, ...
        'doPlots', false);

    POT_both_extreme{loc,1} = out.POT;
    POT_both_extreme{loc,2} = threshold_RF;
    POT_both_extreme{loc,3} = threshold_WL;

    % Identify rainfall-only extreme events
    out2 = find_POT_oneway_RF( ...
        Data, threshold_RF, ...
        'RfAccHours', Rf_Acc, ...
        'compound_window', compound_window, ...
        'dec_tim', dec_tim, ...
        'doPlots', false);

    % Identify water-level-only extreme events
    out3 = find_POT_oneway_NTTWL( ...
        Data, threshold_WL, ...
        'compound_window', compound_window, ...
        'dec_tim', dec_tim, ...
        'doPlots', false);

    POT_RF_only_extreme{loc,1} = out2.POT;
    POT_RF_only_extreme{loc,2} = threshold_RF;

    POT_NTTWL_only_extreme{loc,1} = out3.POT;
    POT_NTTWL_only_extreme{loc,2} = threshold_WL;

end

% Save extracted events and SED estimates
save(fullfile(base_path,'POT_both_extreme.mat'), ...
    'POT_both_extreme', ...
    'POT_RF_only_extreme', ...
    'POT_NTTWL_only_extreme', ...
    'SED_RF', ...
    'SED_WLAM', ...
    'SED_SS');

close all
