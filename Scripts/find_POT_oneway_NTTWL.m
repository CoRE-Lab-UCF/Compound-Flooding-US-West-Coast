function out = find_POT_oneway_NTTWL(Data, threshold_NTTWL, varargin)
% Author: Pravin
% Description: Identifies declustered one-way POT events based on total water level threshold exceedances.

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.addRequired('Data', @(x) isnumeric(x) && size(x,2) >= 3);
p.addRequired('threshold_NTTWL', @(x) isnumeric(x) && isscalar(x));
p.addParameter('compound_window', 1.5, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('dec_tim', 1.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('doPlots', false, @(x) islogical(x) && isscalar(x));
p.addParameter('plotDownsample', 50, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(Data, threshold_NTTWL, varargin{:});
opt = p.Results;

% -------------------------------------------------------------------------
% Unpack data
% -------------------------------------------------------------------------
t_dn = Data(:,1);
nttwl = Data(:,5) + Data(:,6);

% -------------------------------------------------------------------------
% 1. Prepare nTWL water level series
% -------------------------------------------------------------------------
NTTwl = nttwl;

% -------------------------------------------------------------------------
% 2. Decluster total water level peaks
% -------------------------------------------------------------------------
RfAcc_work = NTTwl; % RfAcc_work represents the ntwl here
RfAcc_work(RfAcc_work < opt.threshold_NTTWL) = NaN;

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

% Sort events by peak time
POT = sortrows(POT, 1);
nPOT = size(POT, 1);

% -------------------------------------------------------------------------
% 3. Optional diagnostic plot
% -------------------------------------------------------------------------
if opt.doPlots

    figure('Color','w','Name','One-way POT - Total Water Level');
    plot(t_dn, NTTwl, 'b', 'LineWidth', 0.8); hold on;

    if nPOT > 0
        plot(POT(:,1), POT(:,2), 'or', ...
             'MarkerSize', 7, ...
             'LineWidth', 1.5, ...
             'DisplayName', 'Water-level peaks');
    end

    yline(opt.threshold_NTTWL, '--k', ...
          sprintf('NTTWL threshold = %.3g', opt.threshold_NTTWL), ...
          'LabelHorizontalAlignment','left');

    ylabel('NTTWL (m)');
    title(sprintf('Total Water Level - %d events found', nPOT));
    legend({'NTTWL series','Water-level peaks'}, 'Location','best');
    grid on;
    box on;

end

% -------------------------------------------------------------------------
% 4. Pack outputs
% -------------------------------------------------------------------------
out.POT = POT;
out.RfAcc = NTTwl;
out.th_RF = opt.threshold_NTTWL;
out.nPOT = nPOT;

end
