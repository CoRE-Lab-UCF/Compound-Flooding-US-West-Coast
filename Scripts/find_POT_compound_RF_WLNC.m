function out = find_POT_compound_RF_WLNC(Data, threshold_RF, threshold_Hm0, varargin)
%FIND_POT_COMPOUND  Find declustered compound POT events (AccRF + Hm0).
%
% Usage
%   out = find_POT_compound(Data, threshold_RF, threshold_Hm0)
%   out = find_POT_compound(Data, threshold_RF, threshold_Hm0, ...
%             'RfAccHours', 24, 'compound_window', 1.5, 'dec_tim', 2.5, ...
%             'doPlots', true)
%
% Inputs
%   Data(:,1)      = time (datenum)
%   Data(:,2)      = hourly rainfall (RF)
%   Data(:,3)      = significant wave height Hm0 (m)
%   threshold_RF   = minimum accumulated RF for a compound event
%   threshold_Hm0  = minimum Hm0 for a compound event
%
% Options (name-value)
%   RfAccHours       moving-sum window in hours for RF accumulation (default 24)
%   compound_window  half-window in days to search for Hm0 peak around RF peak (default 1.5)
%   dec_tim          decluster half-window in days centred on RF peak (default 2.5)
%   doPlots          true/false diagnostic plots (default false)
%   plotDownsample   plot every Nth point for long series (default 50)
%
% Output (struct out)
%   out.POT    = [t_RF_peak, AccRF_peak, t_Hm0_peak, Hm0_peak]  (nEvents x 4)
%   out.RfAcc  = accumulated RF time series (same length as Data)
%   out.th_RF  = threshold_RF  (echo)
%   out.th_Hm0 = threshold_Hm0 (echo)
%   out.nPOT   = number of compound events found

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
p.addRequired('Data',         @(x) isnumeric(x) && size(x,2) >= 3);
p.addRequired('threshold_RF', @(x) isnumeric(x) && isscalar(x));
p.addRequired('threshold_Hm0',@(x) isnumeric(x) && isscalar(x));
p.addParameter('RfAccHours',      24,    @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('compound_window', 1.5,   @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('dec_tim',         2.5,   @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('doPlots',         false, @(x) islogical(x) && isscalar(x));
p.addParameter('plotDownsample',  50,    @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(Data, threshold_RF, threshold_Hm0, varargin{:});
opt = p.Results;

% -------------------------------------------------------------------------
% Unpack data
% -------------------------------------------------------------------------
t_dn = Data(:,1);
RF   = Data(:,2);
SS_plus_R2   = Data(:,5)+Data(:,6); % storm surge + Runnup

% -------------------------------------------------------------------------
% 1. Accumulate RF (causal: current hour + past RfAccHours-1 samples)
% -------------------------------------------------------------------------
RfAcc = movsum(RF, [opt.RfAccHours - 1, 0], 'omitnan');

% -------------------------------------------------------------------------
% 2. Decluster AccRF series and find compound events
%    - Iteratively pick max AccRF above threshold_RF
%    - Blank ±dec_tim around each picked peak (whether compound or not)
%    - Check if max Hm0 within ±compound_window of that RF peak >= threshold_Hm0
% -------------------------------------------------------------------------
RfAcc_work = RfAcc;                        % working copy — gets blanked
RfAcc_work(RfAcc_work < opt.threshold_RF) = NaN;  % ignore sub-threshold from start

POT = zeros(0, 4);   % [t_RF, AccRF, t_Hm0, Hm0]

while true
    % -- Pick dominant AccRF peak ------------------------------------------
    [mx_rf, imx] = max(RfAcc_work);
    if ~isfinite(mx_rf), break; end

    t_rf = t_dn(imx);

    % -- Search for Hm0 peak within ±compound_window of RF peak -----------
    hWin   = (t_dn >= t_rf - opt.compound_window) & (t_dn <= t_rf + opt.compound_window);
    SS_sub = SS_plus_R2(hWin);
    t_sub  = t_dn(hWin);

    % -- Blank-out window is always dec_tim before first peak to
    %    dec_tim after last peak (dynamic, uses both peak times) ----------
    if isempty(SS_sub) || all(isnan(SS_sub))
        % No Hm0 data in window — blank around RF peak only and skip
        inWin = (t_dn >= t_rf - opt.dec_tim) & (t_dn <= t_rf + opt.dec_tim);
        RfAcc_work(inWin) = NaN;
        continue
    end

    [mx_hs, ihs] = max(SS_sub);
    t_hs = t_sub(ihs);

    if mx_hs < opt.threshold_Hm0
        % Hm0 does not exceed threshold → not compound, blank ±dec_tim around RF peak only
        inWin = (t_dn >= t_rf - opt.dec_tim) & (t_dn <= t_rf + opt.dec_tim);
        RfAcc_work(inWin) = NaN;
        continue
    end

    % Confirmed compound event — dynamic blank:
    % dec_tim before the earlier peak to dec_tim after the later peak
    t_first = min(t_rf, t_hs);
    t_last  = max(t_rf, t_hs);
    inWin   = (t_dn >= t_first - opt.dec_tim) & (t_dn <= t_last + opt.dec_tim);
    RfAcc_work(inWin) = NaN;

    % -- Store compound event ---------------------------------------------
    POT(end+1, :) = [t_rf, mx_rf, t_hs, mx_hs]; %#ok<AGROW>
end

% Sort by RF peak time
POT  = sortrows(POT, 1);
nPOT = size(POT, 1);

% -------------------------------------------------------------------------
% 3. Optional diagnostic plots
% -------------------------------------------------------------------------
if opt.doPlots
    
    % -- Plot 1: Accumulated RF + compound RF peaks -----------------------
    figure('Color','w','Name','Compound POT — Accumulated RF');
    plot(t_dn, RfAcc, 'b', 'LineWidth', 0.8); hold on;
    if nPOT > 0
        plot(POT(:,1), POT(:,2), 'or', 'MarkerSize', 7, ...
             'LineWidth', 1.5, 'DisplayName', 'Compound RF peaks');
    end
    yline(opt.threshold_RF, '--k', sprintf('RF threshold = %.3g', opt.threshold_RF), ...
          'LabelHorizontalAlignment','left');
    ylabel(sprintf('Accumulated RF  (%d-hr)', opt.RfAccHours));
    title(sprintf('Accumulated RF — %d compound events found', nPOT));
    legend({'AccRF series','Compound RF peaks'}, 'Location','best');
    grid on; box on;

    % -- Plot 2: Hm0 + compound Hm0 peaks ---------------------------------
    figure('Color','w','Name','Compound POT — Hm0');
    plot(t_dn, SS_plus_R2, 'b', 'LineWidth', 0.8); hold on;
    if nPOT > 0
        plot(POT(:,3), POT(:,4), 'or', 'MarkerSize', 7, ...
             'LineWidth', 1.5, 'DisplayName', 'Compound Hm0 peaks');
    end
    yline(opt.threshold_Hm0, '--k', sprintf('Hm0 threshold = %.3g', opt.threshold_Hm0), ...
          'LabelHorizontalAlignment','left');
    ylabel('Hm0 (m)');
    title(sprintf('Hm0 — %d compound events found', nPOT));
    legend({'Hm0 series','Compound Hm0 peaks'}, 'Location','best');
    grid on; box on;
end

% -------------------------------------------------------------------------
% 4. Pack outputs
% -------------------------------------------------------------------------
out.POT    = POT;          % [t_RF_peak, AccRF_peak, t_Hm0_peak, Hm0_peak]
out.RfAcc  = RfAcc;
out.th_RF  = opt.threshold_RF;
out.th_Hm0 = opt.threshold_Hm0;
out.nPOT   = nPOT;

end