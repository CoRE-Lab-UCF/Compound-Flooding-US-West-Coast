% Author: Pravin
% Description: Stratifies compound POT events into AR and non-AR events using the AR catalog.

clear
clc


%% Analysis settings
Rf_Acc = 24;              % Rainfall accumulation period (hours)
compound_window = 1.5;    % Half-width of compound event window (days)
dec_tim = 1.5;            % Half-width of declustering window (days)

RF_percentile = 90;
nTWL_percentile = 90;

%% Load time-series data
load('Time_series_data.mat');

%% Preallocate outputs
nSites = numel(shp);

POT_both_extreme = cell(nSites, 3);
POT_RF_only_extreme = cell(nSites, 2);
POT_NTTWL_only_extreme = cell(nSites, 2);

%% Identify POT events for each site
for loc = 1:nSites

    Data = Time_series_data{loc,1};

    % Calculate daily rainfall totals
    dayIndex = floor(Data(:,1));
    [uniqueDays, ~, idx] = unique(dayIndex);

    dailyRainfall = [ ...
        uniqueDays, ...
        accumarray(idx, Data(:,2))];

    % Define rainfall and nTWL thresholds
    threshold_RF = prctile( ...
        dailyRainfall(:,2), ...
        RF_percentile);

    nTWL = Data(:,5) + Data(:,6);

    threshold_nTWL = prctile( ...
        nTWL, ...
        nTWL_percentile);

    % Compound RF-nTWL POT events
    out = find_POT_compound_RF_WLNC( ...
        Data, ...
        threshold_RF, ...
        threshold_nTWL, ...
        'RfAccHours', Rf_Acc, ...
        'compound_window', compound_window, ...
        'dec_tim', dec_tim, ...
        'doPlots', false);

    POT_both_extreme{loc,1} = out.POT;
    POT_both_extreme{loc,2} = threshold_RF;
    POT_both_extreme{loc,3} = threshold_nTWL;

    % RF-only POT events
    out_RF = find_POT_oneway_RF( ...
        Data, ...
        threshold_RF, ...
        'RfAccHours', Rf_Acc, ...
        'compound_window', compound_window, ...
        'dec_tim', dec_tim, ...
        'doPlots', false);

    POT_RF_only_extreme{loc,1} = out_RF.POT;
    POT_RF_only_extreme{loc,2} = threshold_RF;

    % nTWL-only POT events
    out_nTWL = find_POT_oneway_NTTWL( ...
        Data, ...
        threshold_nTWL, ...
        'compound_window', compound_window, ...
        'dec_tim', dec_tim, ...
        'doPlots', false);

    POT_NTTWL_only_extreme{loc,1} = out_nTWL.POT;
    POT_NTTWL_only_extreme{loc,2} = threshold_nTWL;

end

%% Save outputs
save( ...
    'POT_both_extreme.mat', ...
    'POT_both_extreme', ...
    'POT_RF_only_extreme', ...
    'POT_NTTWL_only_extreme');

disp('POT event identification complete.');

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

%% Load data (or use the ones generated above)
% Required files:
%   AEP_Grid_all.mat
%   AR_data.mat
%   non_AR_data.mat

load('AEP_Grid_all.mat');
load('AR_data.mat');
load('non_AR_data.mat');

%% Analysis settings
doneSites = 1:77;

aepLevels = [0.5 0.2 0.1 0.02 0.01];
aepLabels = {'0.5','0.2','0.1','0.02','0.01'};

grid_all = Grid(:,:,1);

% Folder where figures and calculated values will be saved
outRoot = fullfile(pwd, 'Figures_AEP_AllSites');

if ~exist(outRoot, 'dir')
    mkdir(outRoot);
end

%% Calculate AR contribution for each site
nSites = numel(doneSites);

RC_all = cell(nSites,1);
marginalx_probs = cell(nSites,1);
marginaly_probs = cell(nSites,1);

for k = 1:nSites

    site = doneSites(k);

    fprintf('Processing site %d of %d (site ID %d)\n', ...
        k, nSites, site);

    out = calc_and_plot_AEP_AR_nonAR_site( ...
        'Site', site, ...
        'Grid', grid_all, ...
        'AEP_P1', AEP_P1, ...
        'AEP_P2', AEP_P2, ...
        'UseFullAEP', true, ...
        'MLE_full', MLE, ...
        'UseFullMLE', true, ...
        'Ensemble', Ensemble, ...
        'AR_data', AR_data, ...
        'non_AR_data', non_AR_data, ...
        'aep', aepLevels, ...
        'AAEP', aepLabels, ...
        'OutRoot', outRoot, ...
        'x_U_lim', 20);

    RC_all{k} = out.RC_AR_at_MLE_pct;
    marginalx_probs{k} = out.probs_xmrginal;
    marginaly_probs{k} = out.probs_ymrginal;

end

%% Save results
resultsFile = fullfile(outRoot, 'RC_at_MLE_ALL_SITES.mat');

save(resultsFile, ...
    'RC_all', ...
    'marginalx_probs', ...
    'marginaly_probs', ...
    'doneSites', ...
    'aepLevels');

%% Convert cell outputs to matrices
average_MLE = NaN(nSites, numel(aepLevels));
X_marginal = NaN(nSites, numel(aepLevels));
Y_marginal = NaN(nSites, numel(aepLevels));

for k = 1:nSites

    average_MLE(k,:) = RC_all{k}(:)';
    X_marginal(k,:) = marginalx_probs{k}(:)';
    Y_marginal(k,:) = marginaly_probs{k}(:)';

end

%% Store AR contributions by return period
all_MLE_2yr   = NaN(max(doneSites),1);
all_MLE_5yr   = NaN(max(doneSites),1);
all_MLE_10yr  = NaN(max(doneSites),1);
all_MLE_50yr  = NaN(max(doneSites),1);
all_MLE_100yr = NaN(max(doneSites),1);

all_MLE_2yr(doneSites)   = average_MLE(:,1);
all_MLE_5yr(doneSites)   = average_MLE(:,2);
all_MLE_10yr(doneSites)  = average_MLE(:,3);
all_MLE_50yr(doneSites)  = average_MLE(:,4);
all_MLE_100yr(doneSites) = average_MLE(:,5);

disp('Analysis complete.');
