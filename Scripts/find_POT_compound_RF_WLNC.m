function out = find_POT_compound_RF_WLNC(Data, threshold_RF, threshold_Hm0, varargin)
% Author: Pravin
% Description: Identifies declustered compound POT events based on accumulated rainfall and water level thresholds.

% Parse inputs
p = inputParser;
p.addRequired('Data', @(x) isnumeric(x) && size(x,2) >= 3);
p.addRequired('threshold_RF', @(x) isnumeric(x) && isscalar(x));
p.addRequired('threshold_Hm0', @(x) isnumeric(x) && isscalar(x));
p.addParameter('RfAccHours', 24, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('compound_window', 1.5, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('dec_tim', 2.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('doPlots', false, @(x) islogical(x) && isscalar(x));
p.addParameter('plotDownsample', 50, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(Data, threshold_RF, threshold_Hm0, varargin{:});
opt = p.Results;

% Unpack data
t_dn = Data(:,1);
RF = Data(:,2);
SS_plus_R2 = Data(:,5) + Data(:,6); % Storm surge + Runnup

% Accumulate rainfall over the selected moving window
RfAcc = movsum(RF, [opt.RfAccHours - 1, 0], 'omitnan');

% Decluster accumulated rainfall peaks and identify compound events
RfAcc_work = RfAcc;
RfAcc_work(RfAcc_work < opt.threshold_RF) = NaN;

POT = zeros(0,4);

while true

    [mx_rf, imx] = max(RfAcc_work);
    if ~isfinite(mx_rf)
        break
    end

    t_rf = t_dn(imx);

    % Search for the water-level peak within the compound event window
    hWin = (t_dn >= t_rf - opt.compound_window) & ...
           (t_dn <= t_rf + opt.compound_window);

    SS_sub = SS_plus_R2(hWin);
    t_sub = t_dn(hWin);

    if isempty(SS_sub) || all(isnan(SS_sub))
        inWin = (t_dn >= t_rf - opt.dec_tim) & ...
                (t_dn <= t_rf + opt.dec_tim);
        RfAcc_work(inWin) = NaN;
        continue
    end

    [mx_hs, ihs] = max(SS_sub);
    t_hs = t_sub(ihs);

    if mx_hs < opt.threshold_Hm0
        inWin = (t_dn >= t_rf - opt.dec_tim) & ...
                (t_dn <= t_rf + opt.dec_tim);
        RfAcc_work(inWin) = NaN;
        continue
    end

    % Remove the selected event window from further POT selection
    t_first = min(t_rf, t_hs);
    t_last = max(t_rf, t_hs);
    inWin = (t_dn >= t_first - opt.dec_tim) & ...
            (t_dn <= t_last + opt.dec_tim);
    RfAcc_work(inWin) = NaN;

    POT(end+1,:) = [t_rf, mx_rf, t_hs, mx_hs]; %#ok<AGROW>

end

POT = sortrows(POT,1);
nPOT = size(POT,1);

% Optional diagnostic plots
if opt.doPlots

    figure('Color','w','Name','Compound POT - Accumulated RF');
    plot(t_dn, RfAcc, 'b', 'LineWidth', 0.8);
    hold on

    if nPOT > 0
        plot(POT(:,1), POT(:,2), 'or', ...
            'MarkerSize', 7, ...
            'LineWidth', 1.5, ...
            'DisplayName', 'Compound RF peaks');
    end

    yline(opt.threshold_RF, '--k', ...
        sprintf('RF threshold = %.3g', opt.threshold_RF), ...
        'LabelHorizontalAlignment','left');

    ylabel(sprintf('Accumulated RF (%d-hr)', opt.RfAccHours));
    title(sprintf('Accumulated RF - %d compound events found', nPOT));
    legend({'AccRF series','Compound RF peaks'}, 'Location','best');
    grid on
    box on

    figure('Color','w','Name','Compound POT - Water Level');
    plot(t_dn, SS_plus_R2, 'b', 'LineWidth', 0.8);
    hold on

    if nPOT > 0
        plot(POT(:,3), POT(:,4), 'or', ...
            'MarkerSize', 7, ...
            'LineWidth', 1.5, ...
            'DisplayName', 'Compound water-level peaks');
    end

    yline(opt.threshold_Hm0, '--k', ...
        sprintf('Water-level threshold = %.3g', opt.threshold_Hm0), ...
        'LabelHorizontalAlignment','left');

    ylabel('Water level');
    title(sprintf('Water level - %d compound events found', nPOT));
    legend({'Water-level series','Compound water-level peaks'}, 'Location','best');
    grid on
    box on

end

% Pack outputs
out.POT = POT;
out.RfAcc = RfAcc;
out.th_RF = opt.threshold_RF;
out.th_Hm0 = opt.threshold_Hm0;
out.nPOT = nPOT;

end
