function [out_proc, ref_proc, out_scanner, ref_scanner, out_raw, ref_raw] = processBrukerRaw(inDir)
% This function processes raw (coil-uncombined) Bruker data so that we can
% compare the spectra processed on the scanner with those that we process
% ourselves.

% This should also help us figure out scaling, coil combination, phasing
% etc.
% Georg Oeltzschner, 2025

%% Get the on-scanner processed data first
[out_scanner, ref_scanner, nav_scanner, coilcombos_scanner] = io_loadspec_bruk_new(inDir, 'n');

% Some still have averages which we just want to add up here (since the
% scanner does just sum them instead of averaging)
if out_scanner.dims.averages ~= 0
    % Don't have a FID-A function to sum averages, so just average and 
    % scale by number of raw averages
    out_scanner = op_averaging(out_scanner);
    out_scanner = op_ampScale(out_scanner, out_scanner.rawAverages);
end

%% Now load and process the raw (coil-uncombined) data
[out_raw, ref_raw, nav_raw, coilcombos_raw] = io_loadspec_bruk_new(inDir, 'y');

if out_raw.dims.coils ~= 0
    % This coil combination normalizes the coil combination coefficients to
    % 1, but I believe Bruker simply sums up the channels, so we'll need to
    % scale up again by the number of channels
    out_proc = op_addrcvrs(out_raw, 1, 'h', coilcombos_raw);
    out_proc = op_ampScale(out_proc, 1/length(coilcombos_raw.ph));

    if out_proc.dims.averages ~= 0
        % Don't have a FID-A function to sum averages, so just average and
        % scale by number of raw averages
        out_proc = op_averaging(out_proc);
        out_proc = op_ampScale(out_proc, out_proc.rawAverages);
    end

else

% Not all 'raw' datasets actually have multiple coils
    if out_raw.dims.averages ~= 0
        % Don't have a FID-A function to sum averages, so just average and
        % scale by number of raw averages
        out_proc = op_averaging(out_raw);
        out_proc = op_ampScale(out_proc, out_proc.rawAverages);
    else
        out_proc = out_raw;
    end

end

%% Ref data (leave as is for now)
ref_proc = ref_raw;

%% Plots

if ~isempty(out_proc) && ~isempty(out_scanner)
figure(1)
subplot(1,2,1);
plot(out_proc.ppm, real(out_proc.specs));
hold on;
plot(out_scanner.ppm, real(out_scanner.specs));
hold off;
set(gca, 'xdir', 'reverse', 'xlim', [0 4.2]);
legend('Out\_Proc', 'Out\_Scanner');

subplot(1,2,2);
plot(out_proc.ppm, real(out_proc.specs)-real(out_scanner.specs));
set(gca, 'xdir', 'reverse', 'xlim', [0 4.2]);
legend('Difference');
end

if ~isempty(ref_proc) && ~isempty(ref_scanner)
figure(2)
subplot(1,2,1);
plot(ref_proc.ppm, real(ref_proc.specs));
hold on;
plot(ref_scanner.ppm, real(ref_scanner.specs));
hold off;
set(gca, 'xdir', 'reverse', 'xlim', [2.7 6.7]);
legend('Ref\_Proc', 'Ref\_Scanner');

subplot(1,2,2);
plot(ref_proc.ppm, real(ref_proc.specs)-real(ref_scanner.specs));
set(gca, 'xdir', 'reverse', 'xlim', [2.7 6.7]);
legend('Difference');
end

end
