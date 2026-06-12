function out = find_POT_oneway_RF(Data, threshold_RF, varargin)
% Author: Pravin
% Description: Identifies declustered one-way POT events based on accumulated rainfall threshold exceedances.

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.addRequired('Data',         @(x) isnumeric(x) && size(x,2) >= 3);
p.addRequired('threshold_RF', @(x) isnumeric(x) && isscalar(x));
p.addParameter('RfAccHours',      24,    @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('compound_window', 1.5,   @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('dec_tim',         1.5,   @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('doPlots',         false, @(x) islogical(x) && isscalar(x));
p.addParameter('plotDownsample',  50,    @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(Data, threshold_RF, varargin{:});
opt = p.Results;

% -------------------------------------------------------------------------
% Unpack data
% -------------------------------------------------------------------------
t_dn = Data(:,1);
RF   = Data(:,2);

% -------------------------------------------------------------------------
% 1. Accumulate rainfall
% -------------------------------------------------------------------------
RfAcc = movsum(RF, [opt.RfAccHours - 1, 0], 'omitnan');

% -------------------------------------------------------------------------
% 2. Decluster accumulated rainfall peaks
% -------------------------------------------------------------------------
RfAcc_work = RfAcc;
RfAcc_work(RfAcc_work < opt.threshold_RF) = NaN;

POT = zeros(0, 2);

while true

    [mx_rf, imx] = max(RfAcc_work);
    if ~isfinite(mx_rf), break; end

    t_rf = t_dn(imx);

    inWin = (t_dn >= t_rf - opt.dec_tim) & ...
            (t_dn <= t_rf + opt.dec_tim);

    RfAcc_work(inWin) = NaN;

    POT(end+1, :) = [t_rf, mx_rf]; %#ok<AGROW>

end

% Sort events by rainfall peak time
POT  = sortrows(POT, 1);
nPOT = size(POT, 1);

% -------------------------------------------------------------------------
% 3. Optional diagnostic plot
% -------------------------------------------------------------------------
if opt.doPlots

    figure('Color','w','Name','One-way POT - Accumulated RF');
    plot(t_dn, RfAcc, 'b', 'LineWidth', 0.8); hold on;

    if nPOT > 0
        plot(POT(:,1), POT(:,2), 'or', 'MarkerSize', 7, ...
             'LineWidth', 1.5, 'DisplayName', 'Rainfall peaks');
    end

    yline(opt.threshold_RF, '--k', sprintf('RF threshold = %.3g', opt.threshold_RF), ...
          'LabelHorizontalAlignment','left');

    ylabel(sprintf('Accumulated RF (%d-hr)', opt.RfAccHours));
    title(sprintf('Accumulated RF - %d rainfall events found', nPOT));
    legend({'AccRF series','Rainfall peaks'}, 'Location','best');
    grid on; box on;

end

% -------------------------------------------------------------------------
% 4. Pack outputs
% -------------------------------------------------------------------------
out.POT    = POT;
out.RfAcc  = RfAcc;
out.th_RF  = opt.threshold_RF;
out.nPOT   = nPOT;

end
