function out = find_POT_oneway_NTTWL(Data, threshold_NTTWL, varargin)
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
p.addRequired('threshold_NTTWL', @(x) isnumeric(x) && isscalar(x));
p.addParameter('compound_window', 1.5,   @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('dec_tim',         1.5,   @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('doPlots',         false, @(x) islogical(x) && isscalar(x));
p.addParameter('plotDownsample',  50,    @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(Data, threshold_NTTWL, varargin{:});
opt = p.Results;

% -------------------------------------------------------------------------
% Unpack data
% -------------------------------------------------------------------------
t_dn = Data(:,1);
nttwl   = [Data(:,5)+Data(:,6)];

% -------------------------------------------------------------------------
% 1. Accumulate RF (causal: current hour + past RfAccHours-1 samples)
% -------------------------------------------------------------------------
NTTwl = nttwl;

% -------------------------------------------------------------------------
% 2. Decluster AccRF series and find compound events
%    - Iteratively pick max AccRF above threshold_RF
%    - Blank ±dec_tim around each picked peak (whether compound or not)

% -------------------------------------------------------------------------
RfAcc_work = NTTwl;                        % working copy — gets blanked
RfAcc_work(RfAcc_work < opt.threshold_NTTWL) = NaN;  % ignore sub-threshold from start

POT = zeros(0, 2);   % [t_RF, AccRF, t_Hm0, Hm0]

while true
    % -- Pick dominant AccRF peak ------------------------------------------
    [mx_rf, imx] = max(RfAcc_work);%I have
    if ~isfinite(mx_rf), break; end

    t_rf = t_dn(imx);

    % -- Blank-out window is always dec_tim before first peak to
    inWin = (t_dn >= t_rf - opt.dec_tim) & (t_dn <= t_rf + opt.dec_tim);
    RfAcc_work(inWin) = NaN;


    % -- Store compound event ---------------------------------------------
    POT(end+1, :) = [t_rf, mx_rf]; %#ok<AGROW>
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
    plot(t_dn, NTTwl, 'b', 'LineWidth', 0.8); hold on;
    if nPOT > 0
        plot(POT(:,1), POT(:,2), 'or', 'MarkerSize', 7, ...
             'LineWidth', 1.5, 'DisplayName', 'Compound RF peaks');
    end
    yline(opt.threshold_NTTWL, '--k', sprintf('RF threshold = %.3g', opt.threshold_NTTWL), ...
          'LabelHorizontalAlignment','left');
    ylabel("NTTWL (m)");
    title(sprintf('Accumulated RF — %d compound events found', nPOT));
    legend({'AccRF series','Compound RF peaks'}, 'Location','best');
    grid on; box on;

end

% -------------------------------------------------------------------------
% 4. Pack outputs
% -------------------------------------------------------------------------
out.POT    = POT;          % [t_RF_peak, AccRF_peak, t_Hm0_peak, Hm0_peak]
out.RfAcc  = NTTwl;
out.th_RF  = opt.threshold_NTTWL;
out.nPOT   = nPOT;

end