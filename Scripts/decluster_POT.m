function [t_out, val_out] = decluster_POT(time, value, threshold, decluster_days, dt_hours)
% DECLUSTER_POT  Peaks-Over-Threshold declustering via iterative censoring.
%
%   Algorithm:
%     1. Find the global maximum of the remaining series.
%     2. Record it as an independent peak.
%     3. NaN-out all values within decluster_days/2 of that peak.
%     4. Repeat until no value exceeds the threshold.
%
%   Inputs:
%     time            - Nx1 sorted MATLAB datenums
%     value           - Nx1 observed values
%     threshold       - scalar; stop when max(value) <= threshold
%     decluster_days  - scalar; full window censored around each peak
%                       (each side = decluster_days/2)
%     dt_hours        - scalar; time step in hours (e.g. 1 for hourly,
%                       3 for 3-hourly, 24 for daily). Default = 1.
%
%   Outputs:
%     t_out           - Mx1 datenums of independent peaks (chronological)
%     val_out         - Mx1 corresponding peak values

    if nargin < 5 || isempty(dt_hours)
        dt_hours = 1;
    end

    v    = value(:);
    t    = time(:);
    N    = numel(t);
    half = decluster_days * 0.5;

    % Number of samples in each half-window (round up to be conservative)
    half_steps = ceil(half * 24 / dt_hours);

    t_out   = zeros(50000, 1);
    val_out = zeros(50000, 1);
    m       = 0;

    while true
        [peak_val, peak_idx] = max(v);

        if peak_val <= threshold
            break
        end

        m = m + 1;
        t_out(m)   = t(peak_idx);
        val_out(m) = peak_val;

        % O(1) index arithmetic — no search needed
        lo = max(1, peak_idx - half_steps);
        hi = min(N, peak_idx + half_steps);

        v(lo:hi) = NaN;
    end

    % Trim pre-allocated output
    t_out   = t_out(1:m);
    val_out = val_out(1:m);

    % Return in chronological order
    [t_out, idx] = sort(t_out);
    val_out      = val_out(idx);
end