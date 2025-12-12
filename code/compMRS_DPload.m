%compMRS_DPload.m
%Jamie Near, Sunnybrook Research Institute, 2025
%Diana Rotaru, Columbia University, 2025
%20251121 updates made to accommodate Windows (i.e. filepaths) and to use only
%individually saved functions (i.e. compMRS_parseBrukerFormat)
% USAGE:
% [out,outw]=compMRS_DPload(DPid);
%
% DESCRIPTION:  
% Simple script to load an entire DP into MATLAB/FID-A.  The function
% returns a m x n cell array where m is the number of subjects in the DP,
% and n is the number of sessions in the DP.  Each element of the cell
% array is a FID-A data structure. 
%
% INPUTS:
% DPid:     Data Packet ID (i.e. DP01, DP02, etc.)
% 
% OUTPUTS:
% out:      m x n cell array where m is the number of subjects in the DP,
%           and n is the number of sessions in the DP.  Each element of the
%           cell array is a water suppressed FID-A data struct. 
% outw:     m x n cell array where m is the number of subjects in the DP,
%           and n is the number of sessions in the DP.  Each element of the
%           cell array is a water unsuppressed FID-A data struct.
% outw_auto:m x n cell array where m is the number of subjects in the DP,
%           and n is the number of sessions in the DP.  Each element of the
%           cell array is a water unsuppressed FID-A data struct acquired
%           during the adjustment.
% All of these include the following structure fields:
%           - flags
%           - fids, specs, sz, ppm, t, averages, subSpecs
%           - spectralwidth, dwelltime, txfrq, dims, BO, pointsToLeftshift
%           - seq, te, tr
%           - date, version, filepath
%

function [out, outw, outw_auto]=compMRS_DPload(DPid)

outw=[];
outw_auto = [];

%First check the vendor to be the same for all subj using DPcheck:
check = compMRS_DPcheck(DPid);

%Loop through subjs and sess.  If bruker, load using compMRS_loadspecBruker (initially, io_loadspec_bruk_new.m).  
% If varian, load using compMRS_loadspecVarian (initially,
% io_loadspec_varian.m).
if strcmp(check.vendor(1),'BRUKER')

    for m = 1:check.nSubj
        subjs=dir([DPid filesep 'sub*']);
        for n = 1:check.nSes(m)
            sess=dir([DPid filesep subjs(m).name filesep 'ses*']);
            
            %Find the MRS data path and the REF data path:
            svspath = dir([DPid filesep subjs(m).name filesep sess(n).name filesep 'mrs' filesep '*svs']);
            refpath = dir([DPid filesep subjs(m).name filesep sess(n).name filesep 'mrs' filesep '*mrsref']);
            
            % need to find out whether the ref data is acquired separately
            % or included with the metabolite data
            % For PV6 and higher where only one directory is expected to 
            % exist if water data was acquired automatically, or two 
            % directories if the water data was acquired separately
            
            % The syntax is different for ParaVision 7: it is usually
            % 'PV-7' while in version 6 it is PV 6 with a space...
            if contains(check.version{1}, ["PV 6", "PV 7", "PV-7", "PV-360", "PV 360"])

				% refscan may be acquired automatically
				[out{m,n}, outw_auto{m,n}]=compMRS_loadspecBruker([svspath(length(svspath)).folder filesep svspath(length(svspath)).name],'y');
				
                % If there was a refscan acquired separatly, load it as well
                if ~isempty(svspath) && ~isempty(refpath) % refscan acquired separately
                    [outw{m,n}]=compMRS_loadspecBruker([refpath(length(refpath)).folder filesep refpath(length(refpath)).name],'y');
                end

            %If PV5, then the water reference scans must be collected
            %separately:
            elseif str2num(extractBetween(check.version{m,n}, 'PV ', '.')) == 5
                if exist([svspath(length(svspath)).folder filesep svspath(length(svspath)).name]) && exist([refpath(length(refpath)).folder filesep refpath(length(refpath)).name]) % refscan acquired separately
                    [out{m,n}]=compMRS_loadspecBruker([svspath(length(svspath)).folder filesep svspath(length(svspath)).name],'y');
                    [outw{m,n}]=compMRS_loadspecBruker([refpath(length(refpath)).folder filesep refpath(length(refpath)).name],'y');
                else 
                    error(['ERROR:  For PV5, must have both svs and ref directories!! Aborting!! '])
                end
            end
        end
    end
elseif strcmp(check.vendor(1),'VARIAN')

    for m = 1:check.nSubj
        subjs=dir([DPid filesep 'sub*']);
        for n = 1:check.nSes(m)
            sess=dir([DPid filesep subjs(m).name filesep 'ses*']);
            %Find the MRS data path and the REF data path:
            svspath = dir([DPid filesep subjs(m).name filesep sess(n).name filesep 'mrs' filesep '*svs']);
            refpath = dir([DPid filesep subjs(m).name filesep sess(n).name filesep 'mrs' filesep '*mrsref']);

            [out{m,n}]=compMRS_loadspecVarian([svspath(length(svspath)).folder filesep svspath(length(svspath)).name]);
            [outw{m,n}]=compMRS_loadspecVarian([refpath(length(refpath)).folder filesep refpath(length(refpath)).name]);
        end
    end
end

   

