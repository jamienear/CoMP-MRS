%compMRS_procBrukerRaw.m
%Jamie Near, Sunnybrook Research Institute, 2025
%Georg Oeltzschner, Johns Hopkins University 2025
%Diana Rotaru, Columbia University, 2025
%
% USAGE:
% [out,outw]=compMRS_procBrukerRaw(inDir);
%
% DESCRIPTION:
% Simple script to processes raw (coil-uncombined) Bruker data
%
% INPUTS:
% inDir     = String variable specifying the path to the scan number
%             directory.  Alternatively, inDir can be an integer specifying
%             number of the scan directory to analyze (assuming it is in the
%             current directory).
%
% OUTPUTS:
% out_proc:      processed water-suppressed spectra
% ref_proc:      processed water-unsuppressed spectra
% Both of these include the following structure fields:
%           - flags
%           - fids, specs, sz, ppm, t, averages, subSpecs
%           - spectralwidth, dwelltime, txfrq, dims, BO, pointsToLeftshift
%           - seq, te, tr
%           - date, version, filepath


function [out_proc, ref_proc] = compMRS_procBrukerRaw(out_raw,outw_raw)


% %% Get the on-scanner processed data first
% [out_scanner, ref_scanner, nav_scanner, coilcombos_scanner, out_raw.isECCed] = compMRS_loadspecBruker(inDir, 'n');
%
% % Some still have averages which we just want to add up here (since the
% % scanner does just sum them instead of averaging)
% if out_scanner.dims.averages ~= 0
%     % Don't have a FID-A function to sum averages, so just average and
%     % scale by number of raw averages
%     out_scanner = op_averaging(out_scanner);
%     out_scanner = op_ampScale(out_scanner, out_scanner.rawAverages);
%
% end
%
% % Apply eddy-current correction
% if ~out_raw.isECCed && ~isempty(ref_scanner)
%     % ... to both scans if EDC_OnOff was set to No
%     [out_scanner, ref_scanner] = op_eccKlose(out_scanner, ref_scanner);
% elseif out_raw.isECCed && ~isempty(ref_scanner)
%     % ... only to the ref scan if EDC_OnOff was set to Yes
%     [~, ref_scanner] = op_eccKlose(out_scanner, ref_scanner);
% end

% %% Now load and process the raw (coil-uncombined) data
% [out_raw, ref_raw, nav_raw, coilcombos_raw, out_raw.isECCed, isRFLed] = compMRS_loadspecBruker(inDir, 'y');
%

% If retro-frequency correction was applied to the scanner-processed data,
% we use the navigator data to align the transients before averaging.
if out_raw.isRFLed == 1
    % If multiple coils, then we first combine the channels of navigator signals
    % Use the scaling and phase factors from the method/reco files.
    if out_raw.dims.coils ~= 0
        out_raw.nav=op_addrcvrs(out_raw.nav, 1, 'h', out_raw.coilcombos, 0);
    end
    [~,freq_drift]=op_freqAlignAverages(out_raw.nav,1,0);
    [out_raw,~]=op_makeFreqDrift(out_raw,freq_drift,0);
else
    % do nothing
end

% There are two cases
% - Case 1: If the reference scan data is available and scanner data was ECCed
%           then we apply ECC before channel combination.
% - Case 2: If the scanner-processed data is not ECCed, then we first
%           combine the coils then perform ECC.

% Just to avoid some crashes later
coilcombos_raw_to_apply = out_raw.coilcombos;
ref_proc = out_raw.ref;

% Eddy-current compensation
if out_raw.isECCed && ~isempty(ref_proc)
    % Always to the reference raw data *and* the metabolite raw data
    [out_raw, ref_proc] = op_eccKlose(out_raw, ref_proc);
    % If ECC is applied to individual channel data, then we don't do any
    % further phase correction.
    coilcombos_raw_to_apply.ph = coilcombos_raw_to_apply.ph*0;
end

% Perform coil combination if multiple receivers
if out_raw.dims.coils ~= 0

    % This coil combination normalizes the coil combination coefficients to
    % 1, but I believe Bruker simply sums up the channels, so we'll need to
    % scale up again by the number of channels
    normalizeWeights = 0;

    out_proc_cc = op_addrcvrs(out_raw, 1, 'h', coilcombos_raw_to_apply, normalizeWeights);

    % Per Thanh (01/24/25):
    % "In PV 6.0 and 7.0 the channels are summed, while in PV 360 V1.0 and
    % onwards the channels are averaged"
    if contains(out_proc_cc.version, ["PV 6", "PV 7", "PV-7"])
        % do nothing, summation is already done by op_addrcvrs
    elseif contains(out_proc_cc.version, ["PV 360", "PV-360"])
        % divide by number of channels to achieve averaging
        out_proc_cc = op_ampScale(out_proc_cc, 1.0/length(coilcombos_raw.ph));
    end

    %We also combine the reference scan if applicable.
    if ~isempty(ref_proc)
        ref_proc = op_addrcvrs(ref_proc, 1, 'h', coilcombos_raw_to_apply, 0);
        if contains(out_proc_cc.version, ["PV 6", "PV 7", "PV-7"])
            % do nothing, summation is already done by op_addrcvrs
        elseif contains(out_proc_cc.version, ["PV 360", "PV-360"])
            % divide by number of channels to achieve averaging
            ref_proc = op_ampScale(ref_proc, 1.0/length(coilcombos_raw.ph));
        end
    end
end

if ~out_raw.isECCed && ~isempty(ref_proc)
    % Perform ECC
    [out_proc_cc, ref_proc] = op_eccKlose(out_proc_cc, ref_proc);
end

% Final phase correction? I still have not figured out if we need one.
if out_raw.isECCed
    out_proc_cc = op_addphase(out_proc_cc, coilcombos_raw.ph(end), 0, 4.68, 1);
end

if ~exist('out_proc_cc')
    out_proc_cc = out_raw;
end

nBadAvgTotal=0;
nbadAverages=1;
rmbadav = 4;
% rmbadav='n';
if rmbadav==0
    out_proc_rm=out_proc_cc;
    nsd='N/A';
else
    sat='n';
    while sat=='n' || sat=='N'
        nsd=rmbadav; %Setting the number of standard deviations;
        iter=1;
        nbadAverages=1;
        nBadAvgTotal=0;
        while nbadAverages>0
            [out_proc_rm,metric{iter},badAverages]=op_rmbadaverages(out_proc_cc,nsd,'t');
            badAverages;
            nbadAverages=length(badAverages)*out_raw.sz(out_raw.dims.averages);
            nBadAvgTotal=nBadAvgTotal+nbadAverages;
            out_proc_cc=out_proc_rm;
            iter=iter+1;
            disp([num2str(nbadAverages) ' bad averages removed on this iteration.']);
            disp([num2str(nBadAvgTotal) ' bad averages removed in total.']);
        end
        %figure('position',[0 50 560 420]);
        %Make figure to show pre-post removal of averages
        h=figure('visible','on');
        subplot(2,1,1);
        plot(out_raw.ppm,real(out_raw.specs(:,:,1)));xlim([1 4]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Before removal of bad averages','FontSize',12);
        box off;
        subplot(2,1,2);
        plot(out_proc_rm.ppm,real(out_proc_rm.specs(:,:,1)));xlim([1 4]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('After removal of bad averages','FontSize',12);
        box off;

        %figure('position',[0 550 560 420]);
        h=figure('visible','on');
        plot([1:length(metric{1})],metric{1},'.r',[1:length(metric{iter-1})],metric{iter-1},'.b','MarkerSize',16);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Deviation Metric','FontSize',10);
        legend('Before rmBadAv','Before rmBadAv','After rmBadAv','After rmBadAv');
        legend boxoff;
        title('Deviation Metric','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 10]);
        close(h);

        %sat1=input('are you satisfied with the removal of bad averages? ','s');
        sat='y';

    end
end

% out_proc_dc=out_proc_rm;

avgAlignDomain = 'f';
if exist('outw_raw') && outw_raw.flags.averaged ~= 1
    outw_proc_aa=op_alignAverages(outw_raw,0.2,'n');
else
    outw_proc_aa=outw_raw;
end
sat='n';
while sat=='n' || sat=='N'
    iter=0;
    iterin=20;
    p=100;
    fscum=zeros(out_proc_rm.sz(2:end));
    phscum=zeros(out_proc_rm.sz(2:end));
    while (abs(p(1))>0.0003 && iter<iterin)
        iter=iter+1
        tmax=0.25+0.03*randn(1);
        ppmmin=1.6+0.1*randn(1);
        ppmmaxarray=[3.5+0.1*randn(1,2),4+0.1*randn(1,3),5.5+0.1*randn(1,1)];
        ppmmax=ppmmaxarray(randi(6,1));
        switch avgAlignDomain
            case 't'
                [out_proc_dc,fs,phs]=op_alignAverages(out_proc_rm,tmax,'y');
            case 'f'
                [out_proc_dc,fs,phs]=op_alignAverages_fd(out_proc_rm,ppmmin,ppmmax,tmax,'y');
            otherwise
                error('ERROR: avgAlignDomain not recognized!');
        end

        x=repmat([1:size(fs,1)]',1,out_proc_dc.sz(out_proc_dc.dims.averages));
        %p=polyfit(x,fs,1);

        fscum=fscum+fs;
        phscum=phscum+phs;

    end
    h=figure('visible','on');
    subplot(2,1,1);
    plot(out_proc_rm.ppm,real(out_proc_rm.specs(:,:,1)));xlim([1 4]);
    set(gca,'FontSize',8);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude(a.u.)','FontSize',10);
    title('Before frequency drift correction','FontSize',12);
    box off;
    subplot(2,1,2);
    plot(out_proc_dc.ppm,real(out_proc_dc.specs(:,:,1)));xlim([1 4]);
    set(gca,'FontSize',8);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude(a.u.)','FontSize',10);
    title('After frequency drift correction','FontSize',12);
    box off;

    h=figure('visible','on');
    plot([1:out_proc_dc.sz(out_proc_dc.dims.averages)],fscum,'.-','LineWidth',2);
    set(gca,'FontSize',8);
    xlabel('Scan Number','FontSize',10);
    ylabel('Frequency Drift [Hz]','FontSize',10);
    box off;
    title('Estimated Frequency Drift','FontSize',12);
    set(h,'PaperUnits','centimeters');
    set(h,'PaperPosition',[0 0 10 10]);
    close(h);

    h=figure('visible','on');
    plot([1:out_proc_dc.sz(out_proc_dc.dims.averages)],phscum,'.-','LineWidth',2);
    set(gca,'FontSize',8);
    xlabel('Scan Number','FontSize',10);
    ylabel('Phase Drift [Deg.]','FontSize',10);
    box off;
    title('Estimated Phase Drift','FontSize',12);
    set(h,'PaperUnits','centimeters');
    set(h,'PaperPosition',[0 0 10 10]);
    close(h);

    sat='y';
    if sat=='n'
        iter=0;
        p=100;
        fscum=zeros(out_proc_rm.sz(2:end));
        phscum=zeros(out_proc_rm.sz(2:end));
        fs2cum=zeros(out_proc_cc.sz(2:end));
        phs2cum=zeros(out_proc_cc.sz(2:end));
        out_proc_cc=out_raw;
    end
    totalFreqDrift=mean(max(fscum)-min(fscum));
    totalPhaseDrift=mean(max(phscum)-min(phscum));
end

% Average the transients
if out_raw.dims.averages ~= 0
    % Don't have a FID-A function to sum averages, so just average and
    % scale by number of raw averages
    out_proc_aa = op_averaging(out_proc_dc);

    % For PV5, these appear to be accumulated
    if contains(out_proc_cc.version, ["PV 5"])
        out_proc_aa = op_ampScale(out_proc_aa, out_proc_aa.rawAverages);
    else
        % PV6 and above: do nothing
    end
end
out_proc = op_autophase(out_proc_aa, 1.7, 2.3);
out_proc.SNR=op_getSNR(out_proc);
out_proc.LW=op_getLW(out_proc,1.8, 2.2,8,'true');

% op_plotspec(out_proc)
% try
%     title_name = extractBetween(out_proc.filepath,'\Data\','_ses-');
% catch
%     disp('filepath to be checked')
% end
% title_name = strrep(title_name,'_','-');
% title(title_name,'FontSize',12);
% xlabel('Frequency (ppm)','FontSize',16);
% ylabel('Amplitude(a.u.)','FontSize',16);


% saveas(gca,[name '-WatSuppCorrfids'],'fig');
% saveas(gca,[name '-WatSuppCorrfids'],'jpg');

%FILLING IN THE FLAGS FOR THE RAW DATA

out.flags.freqcorrected=1;
out.flags.phasecorrected=1;
out.flags.addedrcvrs=1;
out.flags.averaged=1;

end