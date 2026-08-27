function out = calc_and_plot_AEP_AR_nonAR_site(varargin)
% plot_AEP_AR_nonAR_site
%
% Generates joint AEP isolines for AR, non-AR, and combined populations,
% and calculates the relative AR contribution for a single site.
%
% Main outputs:
%   out.RC_AR_at_MLE_pct
%   out.MLE_xy
%   out.X_marginal
%   out.Y_marginal
%   out.probs_xmrginal
%   out.probs_ymrginal
%   out.files

%% Input arguments
p = inputParser;
p.FunctionName = 'plot_AEP_AR_nonAR_site';

% Site and probability grids
p.addParameter('Site', [], ...
    @(x) isnumeric(x) && isscalar(x) && x == round(x) && x >= 1);

p.addParameter('Grid', [], ...
    @(x) isnumeric(x) && ndims(x) >= 2);

p.addParameter('AEP_AR', [], ...
    @(x) isnumeric(x) && size(x,2) == 2);

p.addParameter('AEP_nonAR', [], ...
    @(x) isnumeric(x) && size(x,2) == 2);

% Optional full AEP arrays
p.addParameter('AEP_P1', [], @isnumeric);
p.addParameter('AEP_P2', [], @isnumeric);
p.addParameter('UseFullAEP', false, ...
    @(x) islogical(x) && isscalar(x));

% MLE and ensemble points
p.addParameter('MLE', [], ...
    @(x) isnumeric(x) && size(x,2) == 2);

p.addParameter('MLE_full', [], @isnumeric);
p.addParameter('Ensemble', [], @isnumeric);

p.addParameter('UseFullMLE', false, ...
    @(x) islogical(x) && isscalar(x));

% Event samples
p.addParameter('AR_data', [], @iscell);
p.addParameter('non_AR_data', [], @iscell);

% Output directory
p.addParameter('OutRoot', '', ...
    @(x) (ischar(x) || isstring(x)) && strlength(string(x)) > 0);

% Plot/grid settings
p.addParameter('Resolution', 100, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 10);

p.addParameter('x_L_lim', 0, @isscalar);
p.addParameter('x_U_lim', 15, @isscalar);
p.addParameter('y_L_lim', 0, @isscalar);
p.addParameter('y_U_lim', 300, @isscalar);

p.addParameter('aep', [0.1 0.02 0.01], ...
    @(x) isnumeric(x) && isvector(x) && all(x > 0));

p.addParameter('AAEP', {'0.1','0.02','0.01'}, ...
    @(x) iscell(x) && isvector(x));

% Figure formatting
p.addParameter('FontName', 'Calibri', ...
    @(x) ischar(x) || isstring(x));

p.addParameter('FontSize', 15, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

p.addParameter('FigWidthIn', 10, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

p.addParameter('FigHeightIn', 8, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

p.addParameter('SavePNG', true, ...
    @(x) islogical(x) && isscalar(x));

p.addParameter('SaveFIG', true, ...
    @(x) islogical(x) && isscalar(x));

p.addParameter('PNGdpi', 300, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

p.addParameter('XLabel', 'Storm surge (m)', ...
    @(x) ischar(x) || isstring(x));

p.addParameter('YLabel', '24-hr rainfall (mm)', ...
    @(x) ischar(x) || isstring(x));

p.addParameter('ScatterARColor', 'red', ...
    @(x) ischar(x) || isstring(x));

p.addParameter('ScatterNonARColor', 'blue', ...
    @(x) ischar(x) || isstring(x));

p.addParameter('ScatterSize', 15, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

p.addParameter('IsoLineWidth', 2, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

p.addParameter('SepLineWidth', 3, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

p.addParameter('Colormap', [], ...
    @(x) isempty(x) || isnumeric(x));

p.parse(varargin{:});
S = p.Results;

%% Basic checks
if numel(S.AAEP) ~= numel(S.aep)
    error('AAEP must have the same number of elements as aep.');
end

site = S.Site;

%% Select AEP input
if S.UseFullAEP
    if isempty(S.AEP_P1) || isempty(S.AEP_P2)
        error('UseFullAEP=true requires AEP_P1 and AEP_P2.');
    end

    AEP_AR    = S.AEP_P1;
    AEP_nonAR = S.AEP_P2;

else
    if isempty(S.AEP_AR) || isempty(S.AEP_nonAR)
        error(['Provide AEP_AR and AEP_nonAR, or use ', ...
               'UseFullAEP=true with AEP_P1 and AEP_P2.']);
    end

    AEP_AR    = S.AEP_AR;
    AEP_nonAR = S.AEP_nonAR;
end

%% Select MLE input
if S.UseFullMLE
    if isempty(S.MLE_full)
        error('UseFullMLE=true requires MLE_full.');
    end

    MLE_xy = squeeze(S.MLE_full(:,:,site));

else
    if isempty(S.MLE)
        error('Provide MLE, or use UseFullMLE=true with MLE_full.');
    end

    MLE_xy = S.MLE;
end

%% Output folders
outRoot = char(S.OutRoot);

dirs.combARNA = fullfile(outRoot, '02_Combined_AR_and_nonAR');
dirs.combPop  = fullfile(outRoot, '03_Combined_Population');
dirs.relContr = fullfile(outRoot, '04_Relative_Contribution');

dirNames = fieldnames(dirs);

for k = 1:numel(dirNames)
    thisDir = dirs.(dirNames{k});

    if ~exist(thisDir, 'dir')
        mkdir(thisDir);
    end
end

%% Interpolation grid
grid_all = S.Grid;

xv = linspace(S.x_L_lim, S.x_U_lim, S.Resolution);
yv = linspace(S.y_L_lim, S.y_U_lim, S.Resolution);

[Xg, Yg] = meshgrid(xv, yv);

%% AR and non-AR AEP surfaces
Z_AR = griddata( ...
    grid_all(:,1), grid_all(:,2), AEP_AR(:), Xg, Yg);

Z_nonAR = griddata( ...
    grid_all(:,1), grid_all(:,2), AEP_nonAR(:), Xg, Yg);

Z_AR(Z_AR < eps) = 0;
Z_nonAR(Z_nonAR < eps) = 0;

Z_AR(isnan(Z_AR)) = 0;
Z_nonAR(isnan(Z_nonAR)) = 0;

AR_contours = contours_from_Z(xv, yv, Z_AR, S.aep);
nonAR_contours = contours_from_Z(xv, yv, Z_nonAR, S.aep);

%% Plot AR and non-AR isolines
hSep = make_fig(S.FigWidthIn, S.FigHeightIn);
hold on;

hAR = gobjects(1);
hNonAR = gobjects(1);

for i = 1:numel(S.aep)

    hAR = plot( ...
        AR_contours(i).X, AR_contours(i).Y, ...
        'r', 'LineWidth', S.SepLineWidth);

    hNonAR = plot( ...
        nonAR_contours(i).X, nonAR_contours(i).Y, ...
        'b', 'LineWidth', S.SepLineWidth);
end

xlabel(char(S.XLabel));
ylabel(char(S.YLabel));

set(gca, ...
    'FontSize', S.FontSize, ...
    'FontName', char(S.FontName));

grid on;
box on;

legend([hAR hNonAR], ...
    {'AR events','non-AR events'}, ...
    'Location','best');

files.sepAR = save_all( ...
    hSep, ...
    dirs.combARNA, ...
    sprintf('AR_site_%03d', site), ...
    S);

%% Combined population
Z_pop = 1 - (1 - Z_AR) .* (1 - Z_nonAR);
Z_pop(Z_pop < eps) = NaN;

pop_contours = contours_from_Z(xv, yv, Z_pop, S.aep);

hPop = make_fig(S.FigWidthIn, S.FigHeightIn);
hold on;

for r = 1:numel(S.aep)

    plot( ...
        pop_contours(r).X, ...
        pop_contours(r).Y, ...
        'k', ...
        'LineWidth', S.IsoLineWidth);

    scatter( ...
        MLE_xy(r,1), ...
        MLE_xy(r,2), ...
        45, ...
        'k', ...
        'filled', ...
        '^');

    label_aep( ...
        pop_contours(r).Y, ...
        S.AAEP{r}, ...
        S.FontSize, ...
        char(S.FontName));
end

xlabel(char(S.XLabel));
ylabel(char(S.YLabel));

set(gca, ...
    'FontSize', S.FontSize, ...
    'FontName', char(S.FontName));

grid on;
box on;

files.combPop = save_all( ...
    hPop, ...
    dirs.combPop, ...
    sprintf('Combined_population_site_%03d', site), ...
    S);

%% Relative contribution of AR events
AEP_contrib = (Z_AR ./ Z_pop) * 100;

%% Relative contribution at ensemble points
F = NaN(numel(S.aep), size(S.Ensemble,2));

for mm = 1:numel(S.aep)

    F(mm,:) = interp2( ...
        xv, ...
        yv, ...
        AEP_contrib, ...
        S.Ensemble(mm,:,1,site), ...
        S.Ensemble(mm,:,2,site));
end

RC_at_MLE = median(F, 2, 'omitnan');

%% Event samples used for plotting
AR_plot = [ ...
    S.AR_data{site,1}(:,4), ...
    S.AR_data{site,1}(:,2)];

nonAR_plot = [ ...
    S.non_AR_data{site,1}(:,4), ...
    S.non_AR_data{site,1}(:,2)];

%% Marginal locations
Hmo_threshold = min(AR_plot(:,1));
RF_threshold  = min(AR_plot(:,2));

nAEP = numel(S.aep);

X_marginal = NaN(nAEP,1);
Y_marginal = NaN(nAEP,1);

probs_xmrginal = NaN(nAEP,1);
probs_ymrginal = NaN(nAEP,1);

for k = 1:nAEP

    corX = pop_contours(k).X;
    corY = pop_contours(k).Y;

    [~, idxX] = min(abs(corX - Hmo_threshold));
    [~, idxY] = min(abs(corY - RF_threshold));

    X_marginal(k) = corX(idxY);
    Y_marginal(k) = corY(idxX);

    probs_xmrginal(k) = interp2( ...
        xv, ...
        yv, ...
        AEP_contrib, ...
        corX(idxY), ...
        RF_threshold);

    probs_ymrginal(k) = interp2( ...
        xv, ...
        yv, ...
        AEP_contrib, ...
        Hmo_threshold, ...
        corY(idxX));
end

%% Plot AR relative contribution
hRC = make_fig(S.FigWidthIn, S.FigHeightIn);
hold on;

contourf( ...
    Xg, ...
    Yg, ...
    AEP_contrib, ...
    100, ...
    'LineColor', 'none');

if isempty(S.Colormap)
    colormap(flipud(turbo));
else
    colormap(S.Colormap);
end

clim([0 100]);

c = colorbar;
c.Label.String = 'Relative contribution to AEP from AR (%)';
c.Label.Rotation = 270;
c.Label.VerticalAlignment = 'bottom';
c.Label.FontSize = S.FontSize;

h1 = scatter( ...
    AR_plot(:,1), ...
    AR_plot(:,2), ...
    S.ScatterSize, ...
    char(S.ScatterARColor), ...
    'filled');

h2 = scatter( ...
    nonAR_plot(:,1), ...
    nonAR_plot(:,2), ...
    S.ScatterSize, ...
    char(S.ScatterNonARColor), ...
    'filled');

for r = 1:numel(S.aep)

    plot( ...
        pop_contours(r).X, ...
        pop_contours(r).Y, ...
        'k', ...
        'LineWidth', S.IsoLineWidth);

    scatter( ...
        MLE_xy(r,1), ...
        MLE_xy(r,2), ...
        45, ...
        'k', ...
        'filled', ...
        '^');

    scatter( ...
        S.Ensemble(r,:,1,site), ...
        S.Ensemble(r,:,2,site), ...
        55, ...
        'yellow', ...
        'filled');

    label_aep( ...
        pop_contours(r).Y, ...
        S.AAEP{r}, ...
        S.FontSize, ...
        char(S.FontName));
end

xlabel(char(S.XLabel));
ylabel(char(S.YLabel));

xlim([S.x_L_lim S.x_U_lim]);
ylim([S.y_L_lim S.y_U_lim]);

legend( ...
    [h1 h2], ...
    {'Sample AR','Sample non-AR'}, ...
    'Location','best');

set(gca, ...
    'FontSize', S.FontSize, ...
    'FontName', char(S.FontName));

grid on;
box on;

files.relContr = save_all( ...
    hRC, ...
    dirs.relContr, ...
    sprintf('Relative_contribution_AR_site_%03d', site), ...
    S);

%% Save site-specific relative contribution results
files.RCmat = fullfile( ...
    dirs.relContr, ...
    sprintf('RC_at_MLE_site_%03d.mat', site));

save(files.RCmat, 'RC_at_MLE', 'MLE_xy');

%% Output
out.site = site;
out.MLE_xy = MLE_xy;
out.RC_AR_at_MLE_pct = RC_at_MLE;

out.X_marginal = X_marginal;
out.Y_marginal = Y_marginal;

out.probs_xmrginal = probs_xmrginal;
out.probs_ymrginal = probs_ymrginal;

out.files = files;

end


%% Local functions

function h = make_fig(widthIn, heightIn)

h = figure('Visible','off');

set(h, ...
    'Units','inches', ...
    'Position',[1 1 widthIn heightIn], ...
    'PaperUnits','inches', ...
    'PaperPosition',[0 0 widthIn heightIn]);

end


function files = save_all(h, outDir, baseName, S)

files = struct();

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

if ~ishghandle(h)
    error('save_all:InvalidHandle', ...
        'Input is not a valid graphics handle.');
end

if strcmp(get(h,'Type'), 'figure')
    hFig = h;
else
    hFig = ancestor(h, 'figure');
end

drawnow;

if S.SavePNG
    files.png = fullfile(outDir, [baseName '.png']);

    exportgraphics( ...
        hFig, ...
        files.png, ...
        'Resolution', S.PNGdpi);
end

if S.SaveFIG
    files.fig = fullfile(outDir, [baseName '.fig']);
    savefig(hFig, files.fig);
end

end


function cord = contours_from_Z(xv, yv, Z, levels)

cord = repmat( ...
    struct('X',[],'Y',[]), ...
    numel(levels), ...
    1);

for i = 1:numel(levels)

    C = contourc( ...
        xv, ...
        yv, ...
        Z, ...
        [levels(i) levels(i)]);

    [x,y] = first_segment(C);

    cord(i).X = x;
    cord(i).Y = y;
end

end


function [x,y] = first_segment(C)

x = [];
y = [];

if isempty(C)
    return
end

k = 1;
nPoints = C(2,k);
k = k + 1;

if nPoints <= 0 || (k + nPoints - 1) > size(C,2)
    return
end

x = C(1, k:(k+nPoints-1));
y = C(2, k:(k+nPoints-1));

end


function label_aep(Yvec, aepStr, fontSize, fontName)

if isempty(Yvec)
    return
end

yMax = max(Yvec);

text( ...
    0, ...
    yMax + 10, ...
    ['AEP = ', aepStr], ...
    'FontSize', fontSize, ...
    'FontName', fontName, ...
    'VerticalAlignment','top', ...
    'HorizontalAlignment','left');

end