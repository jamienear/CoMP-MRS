%compMRS_DPproc.m

%Jamie Near, Sunnybrook Research Institute, 2025
%Diana Rotaru, Columbia University, 2025
%Thanh Phong Lê, CIBM Center for Biomedical Imaging and École polytechnique fédérale de Lausanne, 2025
%


% Input: 
% Output:

% USAGE:

% Add the 'code' folder to path then run this code from the data folder.

% [out,outw]=compMRS_DPproc.m


function [out, outw] = compMRS_DPproc(DPid)


mkdir('plots');

% Options for the processing
opt.doECCbeforeAvg = 1; % This is how it usually done with Bruker data
opt.predefCoilAmpl = 1; % Use the predefined coil scaling coefficients, when available
opt.rmBadAvg = 0;

% try
    [check]=compMRS_DPcheck(DPid);
    if check.allSame
        [in, inw, inw_auto] = compMRS_DPload(DPid);

        %Loop through subjects and sessions.
        for m = 1:check.nSubj
            for n = 1:check.nSes(m)
    
                disp(['Processing ' DPid ' sub-' num2str(m) ' ses-' num2str(n)])
    
                % If separate water scan is available
                if ~isempty(inw) && ~isempty(inw{m, n})
                    ident = [DPid '_sub-' num2str(m) '_ses-' num2str(n) '_sepWater'];
                    [out{m, n}, outw{m, n}] = compMRS_DPproc_sub(in{m, n}, inw{m, n}, ident, opt);
                
                % If automatic water scan is available
                elseif ~isempty(inw_auto) && ~isempty(inw_auto{m, n})
                    ident = [DPid '_sub-' num2str(m) '_ses-' num2str(n) '_autoWater'];
                    [out{m, n}, outw{m, n}] = compMRS_DPproc_sub(in{m, n}, inw_auto{m, n}, ident, opt);
                end
            end
        end
    end
% catch
%     disp([DPid ' error'])
% end

end

function [out, outw] = compMRS_DPproc_sub(in_mn, inw_mn, ident, opt)
    % do ECC before averaging (yes/no)
    if opt.doECCbeforeAvg
        [out_mn, outw_mn]=op_eccKlose(in_mn,inw_mn);
    end

    % Do coil combination WITHOUT averaging (if applicable)
    % Here we use the water scan to compute the coefficients
    if ~out_mn.flags.addedrcvrs
        % The code below is mostly from op_combineRcvrs
        
        %first find the weights using the water unsuppressed data:
        coilcombos_to_apply=op_getcoilcombos(outw_mn,2,'h');

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
        
        out_mn=op_addrcvrs(out_mn,2,'h',coilcombos_to_apply);
        outw_mn=op_addrcvrs(outw_mn,2,'h',coilcombos_to_apply);
 
        % Perform some scaling for Bruker, if applicable
        % In PV 360 the channels are averaged, this should be done in case
        % we need to compare to already combined ref scans.
        if contains(in_mn.version, ["PV 360", "PV-360"])
            % divide by number of channels to achieve averaging
            out_mn = op_ampScale(out_mn, 1.0/length(coilcombos_to_apply.ph));
            outw_mn = op_ampScale(ref_proc, 1.0/length(coilcombos_to_apply.ph));
        end
    end

    % do bad averages removal (yes/no) - not implemented yet
    if opt.rmBadAvg
        % op_rmbadaverages
    end

    % Perform partial averaging
    
    % Find the possible blocks sizes for partial averaging
    av_block_sizes = divisors(out_mn.averages);
    
    % Some data are already partially averaged (Varian datasets)
    av_eff_block_sizes = av_block_sizes*out_mn.rawAverages/out_mn.averages;

    % Remove block sizes bigger than 32 (does not really make sense to go higher
    % than that)
    av_block_sizes    =av_block_sizes(av_eff_block_sizes<=32);
    av_eff_block_sizes=av_eff_block_sizes(av_eff_block_sizes<=32);
    
    % Iterate along the list of block sizes. The first element is without
    % block averaging
    for kk=1:length(av_block_sizes)
        out_part_avg = op_blockAvg(out_mn,av_block_sizes(kk));
        
        % do drift correction (if applicable)
        [out_part_avg, fs]=op_freqAlignAverages(out_part_avg, 1, 'n');
        
        % do averaging (if applicable)
        out_part_avg=op_averaging(out_part_avg);
        
        % do EDC if not already done
        if ~opt.doECCbeforeAvg
            [out_part_avg, ~]=op_eccKlose(out_part_avg, outw_mn);
        end
        
        % Compute the quality metrics
        % Get LW (of NAA) and SNR
        out_part_avg = op_autophase(out_part_avg, 1.7, 2.3);
        [FWHM_NAA] = op_getLW(out_part_avg, 1.8, 2.2, 8, 1);
        [SNR]=op_getSNR(out_part_avg,1.8,2.2,-2, 0, 1);
        
        out_part_avg.SNR = SNR;
        out_part_avg.LW = FWHM_NAA;
        out_part_avg.SNR_LW_ratio = SNR/FWHM_NAA;
        out_part_avg.block_size = av_eff_block_sizes(kk);
        out_all{kk}=out_part_avg;
    end

    % Search for the best output
    SNR_LW_ratios = zeros(length(out_all), 0);
    for kk=1:length(out_all)
        SNR_LW_ratios(kk) = out_all{kk}.SNR_LW_ratio;
    end
    [~, index] = max(SNR_LW_ratios);
    disp(['The best result is with block size ' num2str(av_eff_block_sizes(index)) '.'])

    % Output the output with best SNR/LW
    out = out_all{index};
    outw= outw_mn;

    % Plot and save the result to check
    plotlegend = {};
    for ii=1:length(out_all)
        plotlegend{ii} = [num2str(out_all{ii}.block_size) ' blocks'];
    end
    f=figure ('name', ident);
    hold on
    box on
    set(gca, 'XDir','reverse')
    xlabel('Chemical shift (ppm)')
    ylabel('Signal amplitude (a.u.)')

    xlim([0 4.5])

    for ii=1:length(out_all)
        plot(out_all{ii}.ppm, real(out_all{ii}.specs))
    end
    legend(plotlegend)
    saveas(f, ['plots/' ident], 'fig')
    print(f, ['plots/' ident], '-dpng', '-r600');
end