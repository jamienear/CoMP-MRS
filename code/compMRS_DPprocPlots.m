% all_DPs_compMRS_improved.m
% Build an improved all-DP summary figure using compMRS_DPplot preprocessing.
%
% This script:
%   1) Builds a tiled all-DP figure with FIXED layout
%   2) Builds and saves a separate 2x2 figure for 4 user-provided DPs
%
% Outputs:
%   - all_DPs_compMRS_improved.png
%   - all_DPs_compMRS_improved.fig
%   - selected_4DPs_compMRS_2x2.png
%   - selected_4DPs_compMRS_2x2.fig
%
% Optional:
%   - Save each DP panel separately by setting saveIndividualPanels = true

clc
clear
close all

%% User settings
dataFile = 'CoMP-MRS-final.mat';   % or 'out_test.mat'
xRange   = [0 4.5];
xTicks   = 0:1:4;
SD_scaling_factor     = 1.96;
uniformLineWidth      = 1.0;
shadeFaceAlpha        = 0.8;
titleFontSize         = 14;
axisFontSize          = 12;
legendFontSize        = 9;
saveIndividualPanels  = false;  % set true if you also want one file per DP

% ---- Four user-selected DPs for separate 2x2 figure ----
selectedDPs = [11 31 10 20];   % <-- change these to your 4 desired DPs

% ---- Fixed layout for all bounded plots ----
nCols = 4;
nRows = 9;

%% Load data just to identify available DPs
struct = load(dataFile);

if isfield(struct,'out') && ~isempty(struct.out)
    nColsOut = size(struct.out,2);
elseif isfield(struct,'out_auto') && ~isempty(struct.out_auto)
    nColsOut = size(struct.out_auto,2);
else
    error('Neither "out" nor "out_auto" was found in %s.', dataFile);
end

%% Find populated DP indices
validIdx = [];

for idx = 1:nColsOut
    hasThisDP = false;

    if isfield(struct,'out') && ~isempty(struct.out) && idx <= size(struct.out,2)
        if ~isempty(struct.out{1,idx})
            hasThisDP = true;
        end
    end

    if ~hasThisDP && isfield(struct,'out_auto') && ~isempty(struct.out_auto) && idx <= size(struct.out_auto,2)
        if ~isempty(struct.out_auto{1,idx})
            hasThisDP = true;
        end
    end

    if hasThisDP
        validIdx(end+1) = idx; %#ok<AGROW>
    end
end

fprintf('Found %d populated DP entries in %s.\n', numel(validIdx), dataFile);
disp(validIdx)

if isempty(validIdx)
    error('No populated DP entries found.');
end

%% Prepare all-DP layout
maxTiles = nRows * nCols;
nPlots = min(numel(validIdx), maxTiles);

if numel(validIdx) > maxTiles
    warning('Found %d populated DPs, but fixed layout only supports %d. Only the first %d DPs will be plotted.', ...
        numel(validIdx), maxTiles, maxTiles);
end

figW = 1400;
figH = 2100;
f = figure('Color','w','Position',[80 80 figW figH]);

tl = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

plotSummary = struct([]);

%% Loop over all populated DPs (up to 5x7 layout capacity)
for k = 1:nPlots
    DPnum = validIdx(k);
    dpLabel = sprintf('DP%02d', DPnum);

    ax = nexttile(tl);
    hold(ax,'on')

    try
        % ---- REQUIRED: use compMRS_DPplot directly ----
        tmpFig = figure('Visible','off');
        [out_as, mu, sigma] = compMRS_DPplot(dataFile, DPnum);
        if isvalid(tmpFig)
            close(tmpFig);
        end

        % Get ppm reference and number of valid spectra from out_as
        [ppmRef, specMat, ~, ~, nSpec] = summarizeAlignedSpectra(out_as);

        if isempty(specMat) || isempty(mu) || isempty(sigma)
            axis(ax,'off')
            text(ax,0.5,0.5,sprintf('%s (no valid spectra)', dpLabel), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontWeight','bold', ...
                'FontSize',titleFontSize);
            continue
        end

        % Force row vectors for plotting
        ppmRef = ppmRef(:).';
        mu = real(mu(:)).';
        sigma = SD_scaling_factor * real(sigma(:)).';

        % Safety check
        if numel(ppmRef) ~= numel(mu) || numel(mu) ~= numel(sigma)
            axis(ax,'off')
            text(ax,0.5,0.5,sprintf('%s (ppm/mu/sigma size mismatch)', dpLabel), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontWeight','bold', ...
                'FontSize',titleFontSize);
            continue
        end

        % shaded mean ± SD using mu and sigma from compMRS_DPplot
        xp = [ppmRef fliplr(ppmRef)];
        yp = [mu - sigma fliplr(mu + sigma)];

        hp = fill(ax, xp, yp, [0.45 0.45 1.00], ...
            'EdgeColor','none', ...
            'FaceAlpha',shadeFaceAlpha, ...
            'DisplayName','\pm SD');

        hl = plot(ax, ppmRef, mu, 'b', ...
            'LineWidth',uniformLineWidth, ...
            'DisplayName','mean');

        set(ax, 'XDir','reverse')
        xlim(ax, xRange)
        xticks(ax, xTicks)
        ylim(ax, [-0.25 2])

        xlabel(ax, 'Frequency (ppm)')
        ylabel(ax, 'Amplitude [a.u.]')
        title(ax, sprintf('%s  (n = %d)', dpLabel, nSpec), ...
            'FontWeight','bold', 'FontSize',titleFontSize)

        ax.FontSize = axisFontSize;
        ax.LineWidth = 0.75;
        box(ax,'off')
        ax.YAxis.Visible = 'off';

        % lgd = legend(ax, [hl hp], {'mean','\pm SD'}, 'Location','northwest');
        % lgd.FontSize = legendFontSize;
        % lgd.Box = 'off';

        plotSummary(k).index    = DPnum;
        plotSummary(k).dpLabel  = dpLabel;
        plotSummary(k).ppm      = ppmRef;
        plotSummary(k).specMat  = specMat;
        plotSummary(k).meanSpec = mu;
        plotSummary(k).stdSpec  = sigma;
        plotSummary(k).nSpec    = nSpec;

        if saveIndividualPanels
            fOne = figure('Color','w','Position',[120 120 700 450]);
            ax1 = axes(fOne); hold(ax1,'on')

            fill(ax1, xp, yp, [0.45 0.45 1.00], ...
                'EdgeColor','none', ...
                'FaceAlpha',shadeFaceAlpha);

            plot(ax1, ppmRef, mu, 'b', 'LineWidth',1.25);

            set(ax1,'XDir','reverse')
            xlim(ax1, xRange)
            xticks(ax1, xTicks)
            ylim(ax1, [-0.25 2])

            xlabel(ax1,'Frequency (ppm)')
            ylabel(ax1,'Amplitude [a.u.]')
            title(ax1, sprintf('%s  (n = %d)', dpLabel, nSpec), 'FontWeight','bold')
            % legend(ax1, {'\pm SD','mean'}, 'Location','northwest', 'Box','off')
            box(ax1,'off')

            exportgraphics(fOne, sprintf('%s_panel.png', dpLabel), 'Resolution', 300)
            savefig(fOne, sprintf('%s_panel.fig', dpLabel))
            close(fOne)
        end

    catch ME
        % Close any hidden temp figures left behind
        figsNow = findall(groot,'Type','figure');
        for ff = 1:numel(figsNow)
            try
                if ~isequal(figsNow(ff), f)
                    close(figsNow(ff));
                end
            catch
            end
        end

        axis(ax,'off')
        text(ax,0.5,0.5,sprintf('%s\nERROR:\n%s', dpLabel, ME.message), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',9)
    end
end

% Blank any unused tiles in the 5x7 layout
for k = (nPlots + 1):maxTiles
    ax = nexttile(tl);
    axis(ax,'off')
end

sgtitle(tl, sprintf('CoMP-MRS mean \\pm SD spectra across populated DPs (%s)', dataFile), ...
    'FontWeight','bold', 'FontSize',14)

exportgraphics(f, 'all_DPs_compMRS_improved.png', 'Resolution', 300)
savefig(f, 'all_DPs_compMRS_improved.fig')

disp('Done: improved all-DP figure created using compMRS_DPplot processing.');

%% ===== Separate 2x2 figure for four selected DPs =====
if numel(selectedDPs) ~= 4
    error('selectedDPs must contain exactly 4 DP indices.');
end

f2 = figure('Color','w','Position',[120 120 1200 900]);
tl2 = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');

for k = 1:4
    DPnum = selectedDPs(k);
    dpLabel = sprintf('DP%02d', DPnum);

    ax = nexttile(tl2);
    hold(ax,'on')

    try
        tmpFig = figure('Visible','off');
        [out_as, mu, sigma] = compMRS_DPplot(dataFile, DPnum);
        if isvalid(tmpFig)
            close(tmpFig);
        end

        [ppmRef, specMat, ~, ~, nSpec] = summarizeAlignedSpectra(out_as);

        if isempty(specMat) || isempty(mu) || isempty(sigma)
            axis(ax,'off')
            text(ax,0.5,0.5,sprintf('%s (no valid spectra)', dpLabel), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontWeight','bold', ...
                'FontSize',titleFontSize);
            continue
        end

        ppmRef = ppmRef(:).';
        mu = real(mu(:)).';
        sigma = real(sigma(:)).';

        if numel(ppmRef) ~= numel(mu) || numel(mu) ~= numel(sigma)
            axis(ax,'off')
            text(ax,0.5,0.5,sprintf('%s (ppm/mu/sigma size mismatch)', dpLabel), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontWeight','bold', ...
                'FontSize',titleFontSize);
            continue
        end

        xp = [ppmRef fliplr(ppmRef)];
        yp = [mu - sigma fliplr(mu + sigma)];

        hp = fill(ax, xp, yp, [0.45 0.45 1.00], ...
            'EdgeColor','none', ...
            'FaceAlpha',shadeFaceAlpha, ...
            'DisplayName','\pm SD');

        hl = plot(ax, ppmRef, mu, 'b', ...
            'LineWidth',1.25, ...
            'DisplayName','mean');

        set(ax, 'XDir','reverse')
        xlim(ax, xRange)
        xticks(ax, xTicks)
        ylim(ax, [-0.25 2])

        xlabel(ax, 'Frequency (ppm)')
        ylabel(ax, 'Amplitude [a.u.]')
        title(ax, sprintf('%s  (n = %d)', dpLabel, nSpec), ...
            'FontWeight','bold', 'FontSize',titleFontSize)

        ax.FontSize = axisFontSize;
        ax.LineWidth = 0.75;
        box(ax,'off')
        ax.YAxis.Visible = 'off';

        % lgd = legend(ax, [hl hp], {'mean','\pm SD'}, 'Location','northwest');
        % lgd.FontSize = legendFontSize;
        % lgd.Box = 'off';

    catch ME
        figsNow = findall(groot,'Type','figure');
        for ff = 1:numel(figsNow)
            try
                if ~isequal(figsNow(ff), f2)
                    close(figsNow(ff));
                end
            catch
            end
        end

        axis(ax,'off')
        text(ax,0.5,0.5,sprintf('%s\nERROR:\n%s', dpLabel, ME.message), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',9)
    end
end

sgtitle(tl2, sprintf('Selected DPs mean \\pm SD spectra (%s)', dataFile), ...
    'FontWeight','bold', 'FontSize',14)

exportgraphics(f2, 'selected_4DPs_compMRS_2x2.png', 'Resolution', 300)
savefig(f2, 'selected_4DPs_compMRS_2x2.fig')

disp('Done: separate 2x2 figure for selected DPs created.');

%% ===== Local function =====
function [ppmRef, specMat, meanSpec, stdSpec, nSpec] = summarizeAlignedSpectra(out_as)
% Convert compMRS_DPplot output (aligned spectra cell array) into matrix stats.

ppmRef = [];
specMat = [];
meanSpec = [];
stdSpec = [];
nSpec = 0;

if isempty(out_as)
    return
end

if ~iscell(out_as)
    error('Expected compMRS_DPplot output to be a cell array.');
end

for ii = 1:numel(out_as)
    if isempty(out_as{ii}) || ~isstruct(out_as{ii})
        continue
    end

    if ~isfield(out_as{ii}, 'ppm') || ~isfield(out_as{ii}, 'specs')
        continue
    end

    ppm = out_as{ii}.ppm(:).';
    spec = real(out_as{ii}.specs(:)).';

    if isempty(ppmRef)
        ppmRef = ppm;
    end

    if numel(ppm) == numel(ppmRef) && numel(spec) == numel(ppmRef)
        specMat = [specMat; spec]; %#ok<AGROW>
    end
end

if isempty(specMat)
    return
end

meanSpec = mean(specMat, 1, 'omitnan');
stdSpec  = std(specMat, 0, 1, 'omitnan');
nSpec    = size(specMat, 1);
end