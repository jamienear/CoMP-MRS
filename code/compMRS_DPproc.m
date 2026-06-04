%compMRS_DPproc.m
%
%Jamie Near, Sunnybrook Research Institute, 2025
%Diana Rotaru, Columbia University, 2025
%Thanh Phong Lê, CIBM Center for Biomedical Imaging and École polytechnique fédérale de Lausanne, 2025
%Eloise Mougel, CIBM, EPFL 2026
%
% USAGE:
% [out,outw]=compMRS_DPproc(DPid);
%
% DESCRIPTION:  
% 
% 
% INPUTS:
% DPid:     Data Packet ID (i.e. DP01, DP02, etc.)
% opt:      Optional parameters:
%
% 
% OUTPUTS:
% out:      m x n cell array where m is the number of subjects in the DP,
%           and n is the number of sessions in the DP.  Each element of the
%           cell array is a water suppressed FID-A data struct. 
% outw:     m x n cell array where m is the number of subjects in the DP,
%           and n is the number of sessions in the DP.  Each element of the
%           cell array is a water unsuppressed FID-A data struct.
%
% Both outputs are spectrally processed. out{m, n} contains the linewidth and SNR.


function [out, outw, out_auto, outw_auto] = compMRS_DPproc(DPid, opt)
disp(['Processing ' num2str(DPid)])

if ~exist('opt','var')
    % Default options for the processing
    opt.doECCbeforeAvg      = 0; % 1 is how it usually done with Bruker data
    opt.predefCoilAmpl      = 0; % Use the predefined coil scaling coefficients, when available
    opt.rmBadAvg            = 1;
    opt.doBlockAveraging    = 0;
    opt.doDriftCorrection   = 1;
    opt.doCompDriftCorrOnOff= 1;
    opt.iterin              = 20;
    opt.tmaxin              = 0.2;
    opt.aaDomain            = 'f';
    opt.autophase           = 1;
    opt.compFracGroupDelay  = 1;
    % *** EM: handle special cases (07.04.2026)
    opt.exception.extWaterOff.ListDPsub     =   [22,4,1; 22,7,1;... % (DPid, idsubject, idses)
                                                 19,1,1; 19,2,1; 19,3,1; 19,4,1; 19,5,1; 19,6,1; 19,7,1; 19,8,1;...% (DPid, idsubject, idses)
                                                 14,1,1; 14,2,1; 14,3,1; 14,4,1; 14,5,1; 14,6,1; 14,7,1; 14,8,1] ;% EM: handle special cases (07.04.2026) (DPid, idsubject, idses): remove separated water
    opt.exception.discardData.ListDPsub     =   [5,2,1; 16,2,1;...
                                                 17,1,1; 17,2,1; 17,3,1; 17,4,1; 17,5,1; 17,6,1; 17,7,1; 17,8,1] ;% EM: handle special cases (07.04.2026) (DPid, idsubject, idses): remove data (corrupted dataset or poor data quality)
    opt.exception.manualRephase.listIdent     =   {['DP30_sub-8_ses-1'  '_sepWater']}; % % EM: handle special cases (07.04.2026)
    opt.exception.manualRephase.listPhase     =   [-70,1.5, 2 ]; % % EM: handle special cases (07.04.2026) (ph0(°), ph1, ppm0) for each DPid, idsubject in opt.exception.manuallyRephase.ListDPsub 
    % *** EM: handle special cases (07.04.2026)
end


 try
    
    [check]=compMRS_DPcheck(DPid);
    if check.allSame

        % Create a folder to save the plots
        if ~isfolder([pwd filesep 'plots' filesep DPid])
               mkdir(['plots' filesep DPid])
        end

        %---EM 27.02.2026-----
        opt.debug_path = [pwd filesep 'plots' filesep 'debug' ];
        mkdir([opt.debug_path])
        date_ana = char(string(datetime('now'),"yyyy-MM-dd_HH-mm-ss"));
        opt.debug_FileName = ['debug_' num2str(DPid) '_' char(date_ana) '.log'];
        debug_id = fopen([opt.debug_path filesep opt.debug_FileName ], 'w');
        fprintf(debug_id,['Processing ' num2str(DPid)]);
        fclose(debug_id);
        %---EM 27.02.2026 (end)-----

        % --- EM: handle special case (EM 07.04.2026)----
        for ii = 1 : size(opt.exception.extWaterOff.ListDPsub,1)
            opt.exception.extWaterOff.listIdent{ii} = (['DP' num2str(opt.exception.extWaterOff.ListDPsub(ii,1)) '_sub-' num2str(opt.exception.extWaterOff.ListDPsub(ii,2)) '_ses-' num2str(opt.exception.extWaterOff.ListDPsub(ii,3))]);
        end

        for ii = 1 : size(opt.exception.discardData.ListDPsub,1)
            opt.exception.discardData.listIdent{ii} = (['DP' num2str(opt.exception.discardData.ListDPsub(ii,1)) '_sub-' num2str(opt.exception.discardData.ListDPsub(ii,2)) '_ses-' num2str(opt.exception.discardData.ListDPsub(ii,3))]);
        end

        % for ii = 1 : size(opt.exception.manualRephase.ListDPsub,1)
        %     opt.exception.manualRephase.listIdent{ii} = (['DP' num2str(opt.exception.manualRephase.ListDPsub(ii,1)) '_sub-' num2str(opt.exception.manualRephase.ListDPsub(ii,2)) '_ses-' num2str(opt.exception.manualRephase.ListDPsub(ii,3))]);
        % end

        %  --- EM: handle special case (end) (EM 07.04.2026)----
    
        
        [in, inw, inw_auto] = compMRS_DPload(DPid);
        autoWaterExists=~isempty(inw_auto);

        %Loop through subjects and sessions.

        out         = cell(check.nSubj, max(check.nSes));
        outw        = cell(check.nSubj, max(check.nSes));
        out_auto    = cell(check.nSubj, max(check.nSes));
        outw_auto   = cell(check.nSubj, max(check.nSes));

        for m = 1:check.nSubj
            for n = 1:check.nSes(m)
                % disp(['Processing ' DPid ' sub-' num2str(m) ' ses-' num2str(n)]) % <--- EM 27.02.2026-----
                disp_writelog_v26(opt.debug_path, opt.debug_FileName,['Processing ' DPid ' sub-' num2str(m) ' ses-' num2str(n)])  % ---> EM 27.02.2026-----

                nsubj = string(extractBetween(in{m,n}.filepath,[filesep 'sub-0'],[filesep 'ses-']));
                ident = ([DPid '_sub-' nsubj '_ses-' num2str(n)]);
                ident = char(join(ident,""));
                %ident = [DPid '_sub-' num2str(m) '_ses-' num2str(n)];
                % If separate water scan is available
                % if ~isempty(inw) && ~isempty(inw{m, n}) 
                if ~isempty(inw) && ~isempty(inw{m, n}) && sum(contains(opt.exception.extWaterOff.listIdent,ident))==0 && sum(contains(opt.exception.discardData.listIdent,ident))==0  % EM: handle special cases (07.04.2026)
                    disp_writelog_v26(opt.debug_path, opt.debug_FileName,['*** Processing water sep'])  % ---> EM 27.02.2026-----
                    [out{m}{n}, outw{m}{n}] = compMRS_DPproc_sub(in{m, n}, inw{m, n}, [ident  '_sepWater'], check, opt);
                end
                % If automatic water scan is available
                if autoWaterExists && ~isempty(inw_auto{m, n}) && sum(contains(opt.exception.discardData.listIdent,ident))==0  % EM: handle special cases (07.04.2026)
                    disp_writelog_v26(opt.debug_path, opt.debug_FileName,['*** Processing water auto'])  % ---> EM 27.02.2026-----
                    [out_auto{m}{n}, outw_auto{m}{n}] = compMRS_DPproc_sub(in{m, n}, inw_auto{m, n}, [ident  '_autoWater'], check, opt);
                end
            end
        end
    end
catch
    % disp([DPid ' error']) % <--- EM 27.02.2026-----
    disp_writelog_v26(opt.debug_path, opt.debug_FileName,[DPid ' error'])  % ---> EM 27.02.2026-----

end
save(['plots/' DPid])
end

function [out, outw] = compMRS_DPproc_sub(in_mn, inw_mn, ident, check, opt)
   
    % Get number of points to be left-shifted
    ls = in_mn.pointsToLeftshift;
    frac_ls = ls-floor(ls); % On Bruker datasets, this number is not an integer. The fractional part should be corrected as a 1st order phase.

    % Thanh 20260305 - Override lsfid for certain Varian data packets as the
    % procpar file does not contain the right lsfid value to achieve proper phasing
    
    if contains(ident, 'DP12') || contains(ident, 'DP13')
        ls = 0;
    elseif contains(ident, 'DP29')
        ls = 1;
    end

    % Left-shift both the MRS and reference scan.
    out_mn = op_leftshift_keepSize(in_mn, floor(ls));
    outw_mn= op_leftshift_keepSize(inw_mn, floor(ls));

    % Average the water
    outw_mn=op_averaging(outw_mn);

    % On some Varian DPs/subjects, the the reference/working frequency is different
    % between the metabolite and reference scan and it does not seem to be intentional.
    % We try to shift the water scan such that it matches the metabolite scan.
    % This correction is only performed on DP32
    if strcmp(check.vendor(1),'VARIAN') && contains(ident, 'DP32')
        diff_freq = out_mn.txfrq - outw_mn.txfrq;
        if abs(diff_freq) > 20 % Hz
            outw_mn=op_freqshift(outw_mn,diff_freq);
        end
    end

    % do ECC before averaging (yes/no)
    if opt.doECCbeforeAvg
        [out_mn, outw_mn]=op_eccKlose(out_mn, outw_mn);
    end

    % Do coil combination WITHOUT averaging (if applicable)
    % Here we use the water scan to compute the coefficients
    if ~out_mn.flags.addedrcvrs
        % The code below is mostly from op_combineRcvrs
        
        %first find the weights using the water unsuppressed data:
        coilcombos_to_apply=op_getcoilcombos(outw_mn, 2,'h');

        % If we choose to use the predefined weights then set the weights read from the method file (when available)
        if opt.predefCoilAmpl && isfield(out_mn, 'coilcombos')
            coilcombos_to_apply.sig=out_mn.coilcombos.sig;
        end

        % If ECC was applied to individual channel data, then we don't do any
        % further phase correction.
        if opt.doECCbeforeAvg
            coilcombos_to_apply.ph = coilcombos_to_apply.ph*0;
        end

        %Now apply the weights to both the water unsuppressed and water suppressed
        %data, combine the coils, but don't combine the averages:
        
        out_mn  = op_addrcvrs(out_mn,2,'h',coilcombos_to_apply);
        outw_mn = op_addrcvrs(outw_mn,2,'h',coilcombos_to_apply);
 
        % Perform some scaling for Bruker, if applicable
        % In PV 360 the channels are averaged, this should be done in case
        % we need to compare to already combined ref scans.
        if contains(in_mn.version, ["PV 360", "PV-360"])
            % divide by number of channels to achieve averaging
            out_mn  = op_ampScale(out_mn, 1.0/length(coilcombos_to_apply.ph));
            outw_mn = op_ampScale(outw_mn, 1.0/length(coilcombos_to_apply.ph));
        end
    end

    

    % Combine subspectra for SPECIAL, if applicable
    
    if strcmp(out_mn.seq, 'SPECIAL') && out_mn.dims.subSpecs>0
        out_mn=op_combinesubspecs(out_mn,"diff");
    end


    numAveragesPerBlock = out_mn.rawAverages/out_mn.averages;
    % do bad averages removal, if applicable (code from Jamie)
    if opt.rmBadAvg
        out_mn=subBadAveragesRemoval(out_mn, ident, opt);
    end

    % Perform partial averaging, if applicable
    if opt.doBlockAveraging
        % Perform partial averaging
    
        % Find the possible blocks sizes for partial averaging
        av_block_sizes = divisors(out_mn.averages);
        
        % Some data are already partially averaged (Varian datasets)
        av_eff_block_sizes = av_block_sizes*numAveragesPerBlock;
    
        % Remove effective block sizes bigger than 32 (does not really make
        % sense to go higher than that)
        av_block_sizes    =av_block_sizes(av_eff_block_sizes<=32);
        av_eff_block_sizes=av_eff_block_sizes(av_eff_block_sizes<=32);
    
    elseif opt.doCompDriftCorrOnOff
        %Compare with (block size 1) and without (block size = number of
        %avg) drift correction
        av_block_sizes    = [1 out_mn.averages];

        % Some data are already partially averaged (Varian datasets)
        av_eff_block_sizes = av_block_sizes*numAveragesPerBlock;
    
    else
        % Do not perform block averaging
        av_block_sizes    = 1;
        av_eff_block_sizes= numAveragesPerBlock;
    end
    
    % Iterate along the list of block sizes, if applicable

    for kk=1:length(av_block_sizes)
        out_part_avg = op_blockAvg(out_mn, av_block_sizes(kk));
        
        % do drift correction (if applicable) (code from Jamie)
        if opt.doDriftCorrection && out_part_avg.averages>1
            out_part_avg = subDriftCorrection(out_part_avg, ident, opt);
        end

        % do averaging (if applicable)
        out_part_avg=op_averaging(out_part_avg);
        
        % do EDC if not already done
        if ~opt.doECCbeforeAvg
            [out_part_avg, ~]=op_eccKlose(out_part_avg, outw_mn);
        end
        
        % 1st order phase correction to compensate for fractional group
        % delay on Bruker datasets.
        if opt.compFracGroupDelay && abs(frac_ls)>0.00001
            ph1 = -frac_ls*in_mn.dwelltime;
            out_part_avg=op_addphase(out_part_avg, 0, ph1, out_part_avg.centerfreq, 1);
        end
        
        % for debug: rephase Varian data with parameters in procpar file
        %out_part_avg=op_addphase(out_part_avg, -in_mn.rp, -(in_mn.lp/360.0*in_mn.dwelltime), max(in_mn.ppm), 1);  

        % Final phasing (0-order phase)
        if opt.autophase && sum(contains(opt.exception.manualRephase.listIdent,ident))==0
            out_part_avg = op_autophase(out_part_avg, 1.8, 2.2);
             % out_part_avg = op_autophase(out_part_avg, 0.4, 4.1); % EM Test 27.02.2026
        elseif sum(contains(opt.exception.manualRephase.listIdent,ident))==1
            idxTmp = find(contains(opt.exception.manualRephase.listIdent,ident));
            out_part_avg=op_addphase(out_part_avg, opt.exception.manualRephase.listPhase(idxTmp,1), (opt.exception.manualRephase.listPhase(idxTmp,2).*in_mn.dwelltime), opt.exception.manualRephase.listPhase(idxTmp,3), 1);  
        end
        
        % Compute the quality metrics
        % Get LW (of NAA) and SNR
        [~,f_NAA]=op_getPeakHeight(out_part_avg,1.8,2.2);
        [FWHM_NAA,FWHM_NAA_hz] = op_getLW(out_part_avg, f_NAA-0.1, f_NAA+0.1, 8, 1);
        [SNR]=op_getSNR(out_part_avg,f_NAA-0.1,f_NAA+0.1,-2, 0, 1);
        
        out_part_avg.SNR = SNR;
        out_part_avg.LW = FWHM_NAA;
        out_part_avg.LW_hz = FWHM_NAA_hz;
        out_part_avg.SNR_LW_ratio = SNR/FWHM_NAA;
        out_part_avg.block_size = av_eff_block_sizes(kk);
        out_all{kk}=out_part_avg;
    end

    % Search for the best output, if we performed block averaging
    % (otherwise return trivial result)
    SNR_LW_ratios = zeros(length(out_all), 0);
    for kk=1:length(out_all)
        SNR_LW_ratios(kk) = out_all{kk}.SNR_LW_ratio;
    end
    
    [~, index] = max(SNR_LW_ratios);

    if opt.doBlockAveraging || opt.doCompDriftCorrOnOff
        % disp(['The best result is with block size ' num2str(av_eff_block_sizes(index)) '.']) % <--- EM 27.02.2026-----
        disp_writelog_v26(opt.debug_path, opt.debug_FileName,['The best result is with block size ' num2str(av_eff_block_sizes(index)) '.'])  % ---> EM 27.02.2026-----
 
    end

    % Output the output with best SNR/LW
    out = out_all{index};
    outw= outw_mn;

    % Plot and save the result to check
    plotlegend = {};
    for ii=1:length(out_all)
        plotlegend{ii} = [num2str(out_all{ii}.block_size) ' avg/block'];
    end
    
    f=figure ('name', ident);
    hold on
    box on
    set(gca, 'XDir','reverse')
    xlabel('Chemical shift (ppm)')
    ylabel('Signal amplitude (a.u.)')

    xlim([0 4.5])

    for ii=1:length(out_all)
        plot(out_all{ii}.ppm, real(out_all{ii}.specs), 'LineWidth', 3)
    end
    legend(plotlegend)
    ax=gca; ax.FontSize=16;
    saveas(f, ['plots' filesep ident(1:4) filesep ident], 'fig')
    saveas(f, ['plots' filesep ident(1:4) filesep ident], 'png')
    %print(f, ['plots' filesep plots filesep ident], '-dpng', '-r300');
 
end


function output = subBadAveragesRemoval(input, ident, opt)
 
    nBadAvgTotal=0;
    nBadAverages=1;

    sat='n';
    output = input;
    while sat=='n' || sat=='N'
        nsd=4; %Setting the number of standard deviations;
        iter=1;
        nBadAverages=1;
        nBadAvgTotal=0;
        while nBadAverages>0
            [output,metric{iter},badAverages]=op_rmbadaverages(output,nsd,'t');
            nBadAverages=length(badAverages);
            nBadAvgTotal=nBadAvgTotal+nBadAverages;
            iter=iter+1;
            disp([num2str(nBadAverages) ' bad averages removed on this iteration.']);
            disp([num2str(nBadAvgTotal) ' bad averages removed in total.']);
            close all;
        end
        %figure('position',[0 50 560 420]);
        %Make figure to show pre-post removal of averages
        h=figure('visible','off');
        subplot(1,2,1);
        plot(input.ppm,real(input.specs(:,:)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Before','FontSize',12);
        subplot(1,2,2);
        plot(output.ppm,real(output.specs(:,:)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('After','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        % saveas(h,['plots' filesep ident '_rmBadAvg_prePostFig'],'jpg');
        % saveas(h,['plots' filesep ident '_rmBadAvg_prePostFig'],'fig');
        close(h);
        
        %figure('position',[0 550 560 420]);
        h=figure('visible','off');
        plot([1:length(metric{1})],metric{1},'.r',[1:length(metric{iter-1})],metric{iter-1},'x','MarkerSize',16);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Deviation Metric','FontSize',10);
        legend('Before rmBadAv','After rmBadAv');
        legend boxoff;
        title('Deviation Metric','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 10]);
        % saveas(h,['plots' filesep ident '_rmBadAvg_scatterFig'],'png');
        % saveas(h,['plots' filesep ident '_rmBadAvg_scatterFig'],'fig');
        close(h);
        
        %sat1=input('are you satisfied with the removal of bad averages? ','s');
        sat='y';
    end

end


function output = subDriftCorrection(input, ident, opt);
    sat='n';
    out_rm2=input;
    while sat=='n' || sat=='N'
        fsPoly=100;
        phsPoly=1000;
        fscum=zeros(out_rm2.sz(out_rm2.dims.averages),1);
        phscum=zeros(out_rm2.sz(out_rm2.dims.averages),1);
        iter=1;
        while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<opt.iterin
            iter=iter+1;
            close all
            tmax=0.05+0.01*randn(1);
            ppmmin=1.6+0.1*randn(1);
            ppmmaxarray=[3.5+0.1*randn(1,2),4+0.1*randn(1,3),5.5+0.1*randn(1,1)];
            ppmmax=ppmmaxarray(randi(6,1));
            switch opt.aaDomain
                case 't'
                    [out_aa,fs,phs]=op_alignAverages(out_rm2,opt.tmaxin);
                case 'f'
                    [out_aa,fs,phs]=op_alignAverages_fd(out_rm2,ppmmin,ppmmax,tmax,'y');
                otherwise
                    error('ERROR: avgAlignDomain not recognized!');
            end
            
            fsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',fs,1);
            phsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',phs,1);
            %iter
            
            fscum=fscum+fs;
            phscum=phscum+phs;
            
            output=out_aa;
        end
        h=figure('visible','off');
        subplot(1,2,1);
        plot(input.ppm,real(input.specs(:,:)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Before','FontSize',12);
        subplot(1,2,2);
        plot(out_aa.ppm,real(out_aa.specs(:,:)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('After','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        % saveas(h,['plots' filesep ident 'alignAvgs_prePostFig'],'jpg');
        % saveas(h,['plots' filesep ident 'alignAvgs_prePostFig'],'fig');
        close(h);
        
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],fscum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Frequency Drift [Hz]','FontSize',10);
        legend('Frequency Drift','Location','SouthEast');
        legend boxoff;
        title('Estimated Freqeuncy Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        % saveas(h,['plots' filesep ident  '_freqDriftFig'],'jpg');
        % saveas(h,['plots' filesep ident  '_freqDriftFig'],'fig');
        close(h);
        
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],phscum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Phase Drift [Deg.]','FontSize',10);
        legend('Phase Drift','Location','SouthEast');
        legend boxoff;
        title('Estimated Phase Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        % saveas(h,['plots' filesep ident '_phaseDriftFig'],'jpg');
        % saveas(h,['plots' filesep ident '_phaseDriftFig'],'fig');
        close(h);

        sat='y';
        % if sat=='n'
        %     iter=0;
        %     p1=100;
        %     fscum=zeros(out_rm.sz(2:end));
        %     phscum=zeros(out_rm.sz(2:end));
        %     fs2cum=zeros(out_cc.sz(2:end));
        %     phs2cum=zeros(out_cc.sz(2:end));
        %     out_rm2=out_rm;
        %     out_cc2=out_cc;
        % end
        totalFreqDrift=mean(max(fscum)-min(fscum));
        totalPhaseDrift=mean(max(phscum)-min(phscum));
    end
end

