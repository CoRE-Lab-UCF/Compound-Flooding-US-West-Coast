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

%% Load data (or use the once generated above)
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
