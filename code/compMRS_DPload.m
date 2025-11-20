%compMRS_DPload.m
%Jamie Near, Sunnybrook Research Institute, 2025
%Diana Rotaru, Columbia University, 2025
%
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


function [out, outw]=compMRS_DPload(DPid)

%First check the vendor using DPcheck:
check = compMRS_DPcheck(DPid);

%Loop through subjs and sess.  If bruker, load using io_loadspec_bruk_new.m.  
% If varian, load using io_loadspec_varian.m
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
            if str2num(extractBetween(check.version{m,n}, 'PV ', '.')) >= 6 

                if ~isempty(svspath) && ~isempty(refpath) % refscan acquired separately
                    [out{m,n}]=io_loadspec_bruk_new([svspath.folder filesep svspath(length(svspath)).name],'y');
                    [outw{m,n}]=io_loadspec_bruk_new([refpath.folder filesep refpath(length(refpath)).name],'y');
               
                elseif ~isempty(svspath) && isempty(refpath) % refscan acquired automatically
                    [out{m,n},outw{m,n}]=io_loadspec_bruk_new([svspath(length(svspath)).folder filesep svspath(length(svspath)).name],'y');
                
                end

            %If PV5, then the water reference scans must be collected
            %separately:
            elseif str2num(extractBetween(check.version{m,n}, 'PV ', '.')) == 5
                if exist([svspath.folder filesep svspath(length(svspath)).name]) && exist([refpath.folder filesep refpath(length(refpath)).name]) % refscan acquired separately
                    [out{m,n}]=io_loadspec_bruk_new([svspath(length(svspath)).folder filesep svspath(length(svspath)).name],'y');
                    [outw{m,n}]=io_loadspec_bruk_new([refpath(length(refpath)).folder filesep refpath(length(refpath)).name],'y');
               
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

            [out{m,n}]=io_loadspec_varian([svspath.folder filesep svspath.name]);
            [outw{m,n}]=io_loadspec_varian([refpath.folder filesep refpath.name]);
        end
    end
end

   

