% compMRS_DPoverview.m
% Diana Rotaru, Columbia University & Medical University of Vienna, 2025
% Emma Van Praagh, Columbia University & Medical University of Vienna, 2025
%
% DESCRIPTION:
% Gathers an inventory of data information. Outputs all folders saved for
% all DPs, as well as information on vendor, sequence, etc.
%
% USAGE:
% out = compMRS_DPoverview(allDPs_parentfolder)
%
% INPUTS:
% allDPs_parentfolder = directory containing all data packets
%
% OUTPUTS:
% out = excel file containing all data details
%

function out = compMRS_DPoverview(allDPs_parentfolder)

% initialize storage
DPoverview = {};

DPidPath = dir(fullfile(allDPs_parentfolder,'DP*'));
counter = 0;
% for every DP, get DP folder path
for k = 1:numel(DPidPath)
    DP = DPidPath(k).name;

    % for every subject (usually subj 1-8), get subj folder path
    SUBJid = dir(fullfile(allDPs_parentfolder,DP, 'sub*'));
    for kk = 1:numel(SUBJid)
        subj = SUBJid(kk).name;

        % for every session (usually 1 or 1-2), get session folder path
        SESid = dir(fullfile(allDPs_parentfolder,DP,subj,'ses-*'));
        for kkk = 1:numel(SESid)
            ses = SESid(kkk).name;

            % for every scan (usually anat, mrs, other), get scan folder path
            SCANid = dir(fullfile(allDPs_parentfolder,DP,subj,ses));
            SCANid=SCANid([SCANid.isdir]);
            SCANid=SCANid(~ismember({SCANid.name},{'.','..'}));

            for kkkk = 1:numel(SCANid)
                scan = SCANid(kkkk).name;

                % for every acquisition (usually mrsref and svs), get
                % folder path
                ACQid = dir(fullfile(allDPs_parentfolder,DP,subj,ses,scan));
                ACQid = ACQid([ACQid.isdir] & ...
                    ~ismember({ACQid.name}, {'.','..'}));

                % populate table with DP, SUBJ, SES, SCAN, ACQ info for all
                % files received
                for kkkkk = 1:numel(ACQid)
                    acq = ACQid(kkkkk).name;

                    % initialize default values
                    acq_vendor = '';
                    acq_B0 = '';
                    acq_PVversion = '';
                    acq_voxel = '';
                    acq_dataformat = '';
                    acq_seq = '';
                    acq_TE = '';
                    acq_TR = '';
                    acq_NAvg = '';
                    acq_voxel_size = '';
                    acq_voxel_volume = '';

                    % for every mrs scan, get acquisition information and
                    % metabolite and unsuppressed water folders

                    % derive parameters based on data format: Bruker (if 'method' file exists)
                    % or Varian (if 'procpar' file exists)

                    acqPath = fullfile(allDPs_parentfolder,DP,subj,ses, scan, acq);

                    methodPath = fullfile(allDPs_parentfolder,DP,subj,ses, scan, acq, 'method');
                    disp(methodPath)

                    if exist(methodPath, 'file')
                        % Bruker dataset
                        acq_vendor = 'Bruker';

                        % Get field strength
                        methodString = fileread(methodPath);
                        if contains(methodString,'##$PVM_FrqRef=')
                            startIdx = strfind(methodString, '##$PVM_FrqRef=');
                            endIdx = strfind(methodString, '##$PVM_FrqWorkOffset=');
                            ref_freq = methodString(startIdx:endIdx(1));
                            ref_freq = str2double(extractBefore(extractBetween(ref_freq, '##$PVM_FrqRef=( 8 )',' 0 0 0 0 0 0 0'),'.'));
                            acq_B0 = str2num(char(string(round(ref_freq / 42.577,1))));
                        else
                            % For DP1, DP3, DP36 which use PV5 that doesnt
                            % store the reference frequency in the method
                            % file
                            acqpString = fileread(fullfile(allDPs_parentfolder,DP,subj,ses, scan, acq,'acqp'));
                            if contains(acqpString,'##$SFO1=')
                                startIdx = strfind(acqpString, '##$SFO1=');
                                endIdx = strfind(acqpString, '##$O1');
                                ref_freq = acqpString(startIdx:endIdx(1));
                                ref_freq = str2double(extractBetween(ref_freq, '=','.'));
                                acq_B0 = str2num(char(string(round(ref_freq / 42.577,1))));
                            else
                                acq_B0 = 'reference frequency not found';
                            end
                        end

                        % Get PV version
                        methodString = fileread(methodPath);
                        if contains(methodString,'$$ /opt/PV')
                            startIdx = strfind(methodString, '$$ /opt/PV');
                            endIdx = strfind(methodString, '/data/');
                            acq_PVversion = methodString(startIdx:endIdx(1));
                            acq_PVversion = extractAfter(acq_PVversion, '/opt/');
                            acq_PVversion = acq_PVversion(1:end-1);
                        elseif contains(methodString,'$$ process /opt/')
                            startIdx = strfind(methodString, '$$ process /opt/');
                            endIdx = strfind(methodString, '##$Method=');
                            acq_PVversion = methodString(startIdx:endIdx);
                            acq_PVversion = char(extractBetween(acq_PVversion, '/opt/', '/prog/'));
                        else
                            acq_PVversion = 'PV version not found';
                        end

                        % Get pulse sequence name
                        methodString = fileread(methodPath);
                        if contains(methodString,'##$Method=')
                            startIdx = strfind(methodString, '##$Method=');
                            endIdx = strfind(methodString, '##$PVM');
                            acq_seq = methodString(startIdx:endIdx(1)-2);
                            acq_seq = extractAfter(acq_seq,'##$Method=');
                            if contains (acq_seq,'<')
                                acq_seq = extractBetween(acq_seq,'<Bruker:','>');
                            end
                        else
                            acq_seq = 'sequence name not found';
                        end

                        % Get voxel location (as per folder naming) for
                        % both svs and mrsref
                        if isequal(SCANid(kkkk).name, 'mrs')
                            if contains(lower(methodPath), 'hipp')
                                acq_voxel = 'hippocampus';
                            elseif contains(lower(methodPath), 'str')
                                acq_voxel = 'striatum';
                            end


                            % Get data format (rawdata.job0 or fid)
                            rawdataPath = fullfile(acqPath, 'rawdata.job0');
                            if exist(rawdataPath, 'file')
                                acq_dataformat = 'rawdata.job0';
                            elseif exist(fullfile(acqPath, 'fid'))
                                acq_dataformat = 'fid';
                            else
                                acq_dataformat = 'not found';
                            end

                            % Get TE
                            methodString = fileread(methodPath);
                            if contains(methodString,'##$PVM_EchoTime=')
                                startIdx = strfind(methodString, '##$PVM_EchoTime=');
                                endIdx = strfind(methodString, '##$PVM_RepetitionTime=');
                                acq_TE = methodString(startIdx:endIdx(1));
                                acq_TE = round(str2num(char(string(extractBetween(acq_TE,'##$PVM_EchoTime=', '#')))),2);
                            else
                                acq_TE = 'TE not found';
                            end

                            % Get TR
                            methodString = fileread(methodPath);
                            if contains(methodString,'##$PVM_RepetitionTime=')
                                startIdx = strfind(methodString, '##$PVM_RepetitionTime=');
                                endIdx = strfind(methodString, '##$PVM_NAverages=');
                                acq_TR = methodString(startIdx:endIdx(1));
                                acq_TR = char(string(extractBetween(acq_TR,'##$PVM_RepetitionTime=','#')));
                            else
                                acq_TR = 'TR not found';
                            end

                            % Get number of averages
                            methodString = fileread(methodPath);
                            if contains(methodString,'##$PVM_NAverages=')
                                startIdx = strfind(methodString, '##$PVM_NAverages=');
                                endIdx = strfind(methodString, '##$PVM_ScanTime'); % ##$PVM_ScanTime takes care of DP08>sub-01>ses-01>mrs>svs which has fieldmap info saved in the method file
                                if sum(endIdx) == 0
                                    endIdx = strfind(methodString, '##$AverageList=');
                                end
                                acq_NAvg = methodString(startIdx:endIdx(1));
                                acq_NAvg = char(string(extractBetween(acq_NAvg,'##$PVM_NAverages=','#')));
                                if contains(['DP14', 'DP20', 'DP24'],DP) % handle DP14, DP20, DP24 differently
                                    acq_NAvg = char(string(extractBefore(acq_NAvg,'$$ @vis=')));
                                end
                            else
                                acq_NAvg = 'averages number not found';
                            end

                            % Get voxel size and volume
                            methodString = fileread(methodPath);
                            if contains(methodString,'##$PVM_VoxArrSize=')
                                startIdx = strfind(methodString, '##$PVM_VoxArrSize=');
                                endIdx = strfind(methodString, '##$PVM_VoxArrPosition=');
                                acq_voxel_size_temp = methodString(startIdx:endIdx(1));
                                acq_voxel_size_temp = char(string(extractBetween(acq_voxel_size_temp,'##$PVM_VoxArrSize=( 1, 3 )','#')));
                                acq_voxel_size_temp = acq_voxel_size_temp(2:end-1);
                                if contains(['DP14', 'DP18', 'DP20', 'DP24'],DP) % handle DP14, DP18, DP20, DP24 differently
                                    acq_voxel_size_temp = char(string(extractBefore(acq_voxel_size_temp,'$$ @vis=')));
                                end
                                acq_voxel_size = [char(string(extractBefore(acq_voxel_size_temp,' '))) ' x ' char(string(extractBetween(acq_voxel_size_temp,' ',' '))) ' x ' char(string(extractAfter(extractAfter(acq_voxel_size_temp,' '),' ')))];
                            else
                                acq_voxel_size = 'voxel size not found';
                            end
                            acq_voxel_volume = round(str2num(extractBefore(acq_voxel_size_temp,' ')) * str2num(char(string(extractBetween(acq_voxel_size_temp,' ',' ')))) * str2num(extractAfter(extractAfter(acq_voxel_size_temp,' '),' ')),1);

                            % Update spreadsheet info
                            DPoverview = [DPoverview; {
                                DP,subj, ses, scan, acq,...
                                acq_vendor, acq_B0, acq_PVversion, acq_voxel, ...
                                acq_dataformat, acq_seq, acq_TE, acq_TR, acq_NAvg, acq_voxel_size, acq_voxel_volume}];

                        end
                    else
                        % Varian dataset
                        procparPath = fullfile(allDPs_parentfolder,DP,subj,ses, scan, acq, 'procpar');
                        % Varian
                        if exist(procparPath, 'file')
                            disp(procparPath)
                            acq_vendor = 'Varian';
                            acq_B0 = '';
                            acq_PVversion = '';
                            acq_voxel = '';
                            acq_dataformat = '';
                            acq_seq = '';
                            acq_TE = '';
                            acq_TR = '';
                            acq_NAvg = '';
                            acq_voxel_size = '';
                            acq_voxel_volume = '';

                            fidPath = fullfile(acqPath, 'fid');
                            if exist(fidPath, 'file')
                                acq_dataformat = 'fid';
                            else
                                acq_dataformat = 'not found';
                            end
                        end
                        
                        % Update spreadsheet info
                        DPoverview = [DPoverview; {
                            DP,subj, ses, scan, acq,...
                            acq_vendor, acq_B0, acq_PVversion, acq_voxel, ...
                            acq_dataformat, acq_seq, acq_TE, acq_TR, acq_NAvg, acq_voxel_size, acq_voxel_volume}];
                    end
     
                end
            end
        end
    end
end

% Create spreadsheet
excelFile = cell2table(DPoverview, 'VariableNames', { ...
    'DP','subj','ses', 'acq', 'scan', 'vendor', 'B0', 'PV version', 'region', 'data format', ...
    'sequence', 'TE', 'TR', 'nr. averages', 'voxel size', 'voxel volume'});
excelFileName = fullfile('A:\CoMP-MRS\sourcecode\code', 'DPoverview.xlsx');
writetable(excelFile, excelFileName, 'FileType', 'spreadsheet');

out = excelFile;
end