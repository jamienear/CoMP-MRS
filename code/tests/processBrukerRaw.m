function [out_proc, ref_proc, out_scanner, ref_scanner, out_raw, ref_raw] = processBrukerRaw(inDir)
% This function processes raw (coil-uncombined) Bruker data so that we can
% compare the spectra processed on the scanner with those that we process
% ourselves.

% This should also help us figure out scaling, coil combination, phasing
% etc.
% Georg Oeltzschner, 2025

%% Get the on-scanner processed data first
[out_scanner, ref_scanner, nav_scanner, coilcombos_scanner, isECCed] = io_loadspec_bruk_new(inDir, 'n');

% Some still have averages which we just want to add up here (since the
% scanner does just sum them instead of averaging)
if out_scanner.dims.averages ~= 0
    % Don't have a FID-A function to sum averages, so just average and 
    % scale by number of raw averages
    out_scanner = op_averaging(out_scanner);
    out_scanner = op_ampScale(out_scanner, out_scanner.rawAverages);
    
end

% Apply eddy-current correction
if ~isECCed && ~isempty(ref_scanner)
    % ... to both scans if EDC_OnOff was set to No
    [out_scanner, ref_scanner] = op_eccKlose(out_scanner, ref_scanner);
elseif isECCed && ~isempty(ref_scanner)
    % ... only to the ref scan if EDC_OnOff was set to Yes
    [~, ref_scanner] = op_eccKlose(out_scanner, ref_scanner);
end

%% Now load and process the raw (coil-uncombined) data
[out_raw, ref_raw, nav_raw, coilcombos_raw, isECCed] = io_loadspec_bruk_new(inDir, 'y');

if out_raw.dims.coils ~= 0
    % This coil combination normalizes the coil combination coefficients to
    % 1, but I believe Bruker simply sums up the channels, so we'll need to
    % scale up again by the number of channels
    normalizeWeights = 0;
    out_proc = op_addrcvrs(out_raw, 1, 'h', coilcombos_raw, normalizeWeights);

    % Per Tranh (01/24/25):
    % "In PV 6.0 and 7.0 the channels are summed, while in PV 360 V1.0 and
    % onwards the channels are averaged"
    if contains(out_proc.version, ["PV 6", "PV 7", "PV-7"])
        % do nothing, summation is already done by op_addrcvrs
    elseif contains(out_proc.version, ["PV 360", "PV-360"])
        % divide by number of channels to achieve averaging
        out_proc = op_ampScale(out_proc, 1/length(coilcombos_raw.ph));
    end

    if out_proc.dims.averages ~= 0
        % Don't have a FID-A function to sum averages, so just average and
        % scale by number of raw averages
        out_proc = op_averaging(out_proc);
        if contains(out_proc.version, ["PV 5"])
            %out_proc = op_ampScale(out_proc, out_proc.rawAverages);
        elseif contains(out_proc.version, ["PV-360.3.5"])
            %out_proc = op_ampScale(out_proc, 1/out_proc.rawAverages);
        end
    end

else

% Not all 'raw' datasets actually have multiple coils
    if out_raw.dims.averages ~= 0
        % Average
        out_proc = op_averaging(out_raw);
        % For PV5, these appear to be accumulated
        if contains(out_proc.version, ["PV 5"])
            out_proc = op_ampScale(out_proc, out_proc.rawAverages);
        else
            % do nothing
        end
    else
        out_proc = out_raw;
    end

end

% Eddy-current correction
if ~isempty(ref_raw)
    % Always to the reference raw data *and* the metabolite raw data
    [out_proc, ref_proc] = op_eccKlose(out_proc, ref_raw);
else
    ref_proc = ref_raw;
end

% % Final phase correction? I still have not figured out if we need one.
% if isECCed
%     out_proc = op_addphase(out_proc, coilcombos_raw.ph(end), 0, 4.68, 1);
% end

%% Plots

% Get DP string
pathParts = strsplit(inDir, filesep);
DPString = pathParts{contains(pathParts, 'DP')};

if ~isempty(out_proc) && ~isempty(out_scanner)
figure(1)
subplot(1,2,1);
plot(out_proc.ppm, real(out_proc.specs));
hold on;
plot(out_scanner.ppm, real(out_scanner.specs));
hold off;
set(gca, 'xdir', 'reverse', 'xlim', [0 4.2]);
ylims = get(gca, 'YLim');
legend('Out\_Proc', 'Out\_Scanner');

subplot(1,2,2);
plot(out_proc.ppm, real(out_proc.specs)-real(out_scanner.specs));
set(gca, 'xdir', 'reverse', 'xlim', [0 4.2], 'ylim', ylims);
legend('Difference');

sgtitle(DPString);
end

if ~isempty(ref_proc) && ~isempty(ref_scanner)
figure(2)
subplot(1,2,1);
plot(ref_proc.ppm, real(ref_proc.specs));
hold on;
plot(ref_scanner.ppm, real(ref_scanner.specs));
hold off;
set(gca, 'xdir', 'reverse', 'xlim', [2.7 6.7]);
ylims = get(gca, 'YLim');
legend('Ref\_Proc', 'Ref\_Scanner');

subplot(1,2,2);
plot(ref_proc.ppm, real(ref_proc.specs)-real(ref_scanner.specs));
set(gca, 'xdir', 'reverse', 'xlim', [2.7 6.7], 'ylim', ylims);
legend('Difference');

sgtitle([DPString ' ref']);
end

end
