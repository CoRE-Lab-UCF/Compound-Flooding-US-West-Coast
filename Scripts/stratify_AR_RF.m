function [AR, non_AR] = stratify_AR_RF(data, varargin)
% Stratify events by RF criterion only using AR maps
% INPUTS (name-value):
%   data              - numeric array [Time_RF, ...]
%   Column_of_Time_RF - column index of RF peak time
%   AR_Time           - time vector of AR maps
%   AR_lng/AR_ltt     - AR grid lon/lat
%   time_bef_RF       - hours before RF peak to search
%   time_aft_RF       - hours after RF peak to search
%   latt_gauge/long_gauge - gauge location
%   dist_thres_RF     - search radius (km)

p = inputParser;
p.addRequired('data');
p.addParameter('Column_of_Time_RF', 1);
p.addParameter('AR_Time',    []);
p.addParameter('AR_lng',     []);
p.addParameter('AR_ltt',     []);
p.addParameter('time_bef_RF', 0);
p.addParameter('time_aft_RF', 0);
p.addParameter('latt_gauge', NaN);
p.addParameter('long_gauge', NaN);
p.addParameter('dist_thres_RF', NaN);
p.parse(data, varargin{:});
o = p.Results;

nc_file = 'H:\OneDrive - University of Central Florida\Analysis_AR\AR_data\globalARcatalog_ERA5_1940-2024.nc';

% --- Spatial setup ---
[latGrid, lonGrid] = meshgrid(o.AR_ltt, convertLongitudeTo180(o.AR_lng));
dis_km   = deg2km(distance(o.latt_gauge, o.long_gauge, latGrid(:), lonGrid(:)));
ind_sel  = find(dis_km <= o.dist_thres_RF);
if isempty(ind_sel)
    error('No grid points within dist_thres_RF. Please increase it.');
end

% --- Collect unique AR timestep indices across all events ---
all_inds = [];
for i = 1:size(data,1)
    t = data(i, o.Column_of_Time_RF);
    all_inds = [all_inds; find(o.AR_Time >= t - o.time_bef_RF & o.AR_Time <= t + o.time_aft_RF)];
end
unique_inds = unique(all_inds);

% --- Load AR data (block read if contiguous, else individual) ---
AR_cache = containers.Map('KeyType','int32','ValueType','logical');
if ~isempty(unique_inds)
    t_min = min(unique_inds);  t_max = max(unique_inds);
    n_block = t_max - t_min + 1;

    if n_block <= 5 * numel(unique_inds)
        block = ncread(nc_file, 'shapemap', [1,1,1,t_min,1], [inf,inf,1,n_block,1]);
        for k = 1:numel(unique_inds)
            ti = unique_inds(k);
            s  = block(:,:,1,ti-t_min+1,1);
            AR_cache(int32(ti)) = any(~isnan(s(ind_sel)));
        end
        clear block
    else
        for k = 1:numel(unique_inds)
            ti = unique_inds(k);
            s  = ncread(nc_file, 'shapemap', [1,1,1,ti,1], [inf,inf,1,1,1]);
            AR_cache(int32(ti)) = any(~isnan(s(ind_sel)));
        end
    end
end

% --- Classify events ---
AR = []; non_AR = [];
for i = 1:size(data,1)
    t      = data(i, o.Column_of_Time_RF);
    inds   = find(o.AR_Time >= t - o.time_bef_RF & o.AR_Time <= t + o.time_aft_RF);
    is_AR  = any(cellfun(@(k) AR_cache(int32(k)), num2cell(inds)));

    if is_AR
        AR     = [AR;  data(i,:)];
    else
        non_AR = [non_AR; data(i,:)];
    end
end


end