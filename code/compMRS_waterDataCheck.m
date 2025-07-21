% compMRS_waterDataCheck.m
% Emma Van Praagh, Columbia University & Medical University of Vienna, 2025
% Diana Rotaru, Columbia University & Medical University of Vienna, 2025
%
% DESCRIPTION: 
% Gathers an inventory of which water data formats are available in CoMP-MRS
% data packets.
% In the metabolite and unsuppressed water folders, respectively, finds the
    % pulse sequence
    % existence of rawdata.job0 or fid file
    % existence of fid.refscan
    % existence of fid.ref
    % existence of reference printed out in the method file (PVM_RefScan = (Receivers,NA))
    % existence of reference indicated in the method file (PVM_RefScanYN=Yes/No)
% Functionality is limited to Bruker data packets as of 20250701. Varian
% data packet handling is to be added.

% USAGE:
% out = compMRS_waterDataCheck(DPPath,writeToExcelPath)
% Example: compMRS_waterDataCheck('/Users/emmavanpraagh/Documents/CoMP-MRS/DP','/Users/emmavanpraagh/Documents/CoMP-MRS/writeToExcel')
% 
% INPUTS:
% DPPath = directory containing data packets
    % Example: DPPath = '/Users/emmavanpraagh/Documents/CoMP-MRS/DP';
% writeToExcelPath = directory for writing Excel File
    % Example: writeToExcelPath = '/Users/emmavanpraagh/Documents/CoMP-MRS/writeToExcel';
% 
% OUTPUTS:
% out = when succeeds, saves a summary spreadsheet to writeToExcelPath and
%       prints a summary of available and preferred water data files


function out = compMRS_waterDataCheck(DPPath,writeToExcelPath)

% initialize storage
excelStorage = {};

DPidPath = dir(fullfile(DPPath,'DP*'));

% for every DP
for k = 1:numel(DPidPath)
    DPid = DPidPath(k).name;
    
    % for every subject
    subjList = dir(fullfile(DPPath,DPid, 'sub*'));
    for kk = 1:numel(subjList)
        subj = subjList(kk).name;
        
        % for every session, get mrs folder path
        sesList = dir(fullfile(DPPath,DPid,subj,'ses-*'));
        for kkk = 1:numel(sesList)
            ses = sesList(kkk).name;
            mrsPath = fullfile(DPPath,DPid,subj,ses,'mrs');
            if ~isfolder(mrsPath)
                continue;
            end

            allDirs = dir(mrsPath);
            mrsSubExpPaths = allDirs( ...
                [allDirs.isdir] & ...
                ~ismember({allDirs.name}, {'.','..'}) ...
                );
            
            % initialize default values
            refScanYNMet = 'not found';
            refScanYNWat = 'not found';
            rawdataJob0Met = 'not found';
            rawdataJob0Wat = 'not found';
            refScanNAMet = 'not found';
            refScanNAWat = 'not found';
            vendor = 'not found';
            sequenceMet = 'not found';
            sequenceWat = 'not found';
            fidrefscanMet = 'not found';
            fidrefscanWat = 'not found';
            fidrefMet = 'not found';
            fidrefWat = 'not found';
            refPrintedMet = 'not found';
            refPrintedWat = 'not found';
            
            % for every mrs folder, get metabolite and unsuppressed water folders
            for kkkk = 1:numel(mrsSubExpPaths)
                mrsSubExp = mrsSubExpPaths(kkkk).name;
                methodPath = fullfile(mrsPath, mrsSubExp, 'method');
                if exist(methodPath, 'file')
                    % Bruker
                    vendor = 'Bruker';
                    
                    % Check if the folder is svs (metab) or mrsref (water)
                    if contains(lower(mrsSubExp), 'mrsref') || contains(lower(mrsSubExp), 'eddycurrent') % handles DP10
                        %disp(mrsSubExp)
                        % check if rawdata.job0 or fid exists
                        rawdataPath = fullfile(mrsPath, mrsSubExp, 'rawdata.job0');
                        if ~exist(rawdataPath, 'file')
                            rawdataPath = fullfile(mrsPath, mrsSubExp, 'fid');
                        end
                        rawdataExists = exist(rawdataPath, 'file');                        
                        
                        if rawdataExists
                            rawdataJob0Wat = 'Yes';
                        else
                            rawdataJob0Wat = 'No';
                        end
                        
                        % get pulse sequence for water   
                        methodString = fileread(methodPath);
                        if contains(methodString,'##$Method=')
                            startIdx = strfind(methodString, '##$Method=');
                            endIdx = strfind(methodString, '##$PVM');
                            sequenceWat = methodString(startIdx:endIdx);
                        else
                            sequenceWat = 'did not find one'; 
                        end    

                        % start search for fidrefscan in water folder
                        fidrefscanWat = fullfile(mrsPath,mrsSubExp,'fid.refscan');
                        if exist(fidrefscanWat, 'file')
                            fidrefscanWat = 'Yes';
                            %fidrefscanWat = sprintf('Yes: %s', extractAfter(fidrefscanWat,'DP'));
                        else
                           fidrefscanWat = 'No'; 
                        end

                        % start search for fidref in water folder
                        fidrefWat = fullfile(mrsPath,mrsSubExp,'fid.ref');
                        if exist(fidrefWat, 'file')
                           fidrefWat = 'Yes';
                           %fidrefWat = sprintf('Yes: %s', extractAfter(fidrefWat,'DP'));
                        else
                           fidrefWat = 'No'; 
                        end

                        % ref scan indicated yes/no
                        if contains(methodString,'##$PVM_RefScanYN=Yes')
                            refScanYNWat = 'Yes';
                        elseif contains(methodString,'##$PVM_RefScanYN=No')
                            refScanYNWat = 'No';
                        end
    
                        % start search for ref scan printed in water method
                        if contains(methodString, '##$PVM_RefScan=(')
                            refPrintedWat = 'Yes';
                            %refPrintedWat = sprintf('Yes: %s', extractAfter(methodPath,'DP'));
                        else
                            refPrintedWat = 'No';
                        end
                        
                        if contains(methodString,'##$PVM_RefScanNA=')
                            startSearch = strfind(methodString, '##$PVM_RefScanNA=') + length('##$PVM_RefScanNA=');
                            endSearch = strfind(methodString(startSearch:end), '#') + startSearch - 2;
                            if ~isempty(endSearch)
                                refScanNAWat = str2double(methodString(startSearch:endSearch-1));
                            end
                        end
                    else
                    %this might be too limiting based on naming conventions: elseif contains(lower(mrsSubExp), 'svs')
                        
                        % check if rawdata.job0 or fid exists
                        rawdataPath = fullfile(mrsPath, mrsSubExp, 'rawdata.job0');
                        if ~exist(rawdataPath, 'file')
                            rawdataPath = fullfile(mrsPath, mrsSubExp, 'fid');
                        end
                        rawdataExists = exist(rawdataPath, 'file');
                    
                        if rawdataExists
                            rawdataJob0Met = 'Yes';
                        else
                            rawdataJob0Met = 'No';
                        end

                        % get pulse sequence for metab   
                        methodString = fileread(methodPath);
                        if contains(methodString,'##$Method=')
                            startIdx = strfind(methodString, '##$Method=');
                            endIdx = strfind(methodString, '##$PVM');
                            sequenceMet = methodString(startIdx:endIdx);
                        end  

                        % start search for fidrefscan in metab folder
                        fidrefscanMet = fullfile(mrsPath,mrsSubExp,'fid.refscan');
                        if exist(fidrefscanMet, 'file')
                            fidrefscanMet = 'Yes';
                            %fidrefscanMet = sprintf('Yes: %s', extractAfter(fidrefscanMet,'DP'));
                        else
                           fidrefscanMet = 'No'; 
                        end
                        
                        % start search for fidref in metab folder
                        fidrefMet = fullfile(mrsPath,mrsSubExp,'fid.ref');
                        if exist(fidrefMet, 'file')
                            fidrefMet = 'Yes';
                            %fidrefMet = sprintf('Yes: %s', extractAfter(fidrefMet,'DP'));
                        else
                           fidrefMet = 'No'; 
                        end
                        
                        % start search for ref scan printed in metab method
                        if contains(methodString, '##$PVM_RefScan=(')
                            refPrintedMet = 'Yes';
                            %refPrintedMet = sprintf('Yes: %s', extractAfter(methodPath,'DP'));
                        else
                            refPrintedMet = 'No';
                        end
                        
                        % indicated yes/no in method
                        if contains(methodString,'##$PVM_RefScanYN=Yes')
                            refScanYNMet = 'Yes';
                        elseif contains(methodString,'##$PVM_RefScanYN=No')
                            refScanYNMet = 'No';
                        end
                        
                        % averages
                        if contains(methodString,'##$PVM_RefScanNA=')
                            startSearch = strfind(methodString, '##$PVM_RefScanNA=') + length('##$PVM_RefScanNA=');
                            endSearch = strfind(methodString(startSearch:end), '#') + startSearch - 2;
                            if ~isempty(endSearch)
                                refScanNAMet = str2double(methodString(startSearch:endSearch-1));
                            end
                        end  
                    end
                else
                    procParPath = fullfile(mrsPath, mrsSubExp, 'procpar');
                    % Varian
                    if exist(procParPath, 'file')
                        vendor = 'Varian';
                        refScanYNMet = 'Not searched, Varian';
                        refScanYNWat = 'Not searched, Varian';
                        rawdataJob0Met = 'Not searched, Varian';
                        rawdataJob0Wat = 'Not searched, Varian';
                        refScanNAMet = 'Not searched, Varian';
                        refScanNAWat = 'Not searched, Varian';
                        sequenceMet = 'Not searched, Varian';
                        sequenceWat = 'Not searched, Varian';
                        fidrefscanMet = 'Not searched, Varian';
                        fidrefscanWat = 'Not searched, Varian';
                        fidrefMet = 'Not searched, Varian';
                        fidrefWat = 'Not searched, Varian';
                    end
                end
                
            end
            
            % available water summary
            availableWater = {};
            preferredWater = {};

            if contains(vendor, 'Bruker')
                if contains(rawdataJob0Wat, 'Yes')
                    if isempty(preferredWater)
                        preferredWater{end+1} = 'Unsup.Water.Acq.YN';
                    end
                    availableWater{end+1} = 'Unsup.Water.Acq.YN';
                end
                if contains(fidrefscanMet, 'Yes')
                    % pick it up here
                    if isempty(preferredWater)
                        preferredWater{end+1} = 'Metab.FidRefScan';
                    end
                    availableWater{end+1} = 'Metab.FidRefScan';
                end
                if contains(fidrefMet, 'Yes')
                    if isempty(preferredWater)
                        preferredWater{end+1} = 'Metab.FidRef';
                    end
                    availableWater{end+1} = 'Metab.FidRef';
                end
                if contains(fidrefscanWat, 'Yes')
                    if isempty(preferredWater)
                        preferredWater{end+1} = 'Unsup.Water.FidRefScan';
                    end                   
                    availableWater{end+1} = 'Unsup.Water.FidRefScan';
                end
                if contains(fidrefWat, 'Yes')
                    if isempty(preferredWater)
                        preferredWater{end+1} = 'Unsup.Water.FidRef';
                    end                    
                    availableWater{end+1} = 'Unsup.Water.FidRef';
                end
                if contains(refPrintedMet, 'Yes')
                    if isempty(preferredWater)
                        preferredWater{end+1} = 'Metab.Method.Printed.RefScan';
                    end
                    availableWater{end+1} = 'Metab.Method.Printed.RefScan';
                end
                if contains(refPrintedWat, 'Yes')
                    if isempty(preferredWater)
                        preferredWater{end+1} = 'Unsup.Water.Printed.RefScan';
                    end                    
                    availableWater{end+1} = 'Unsup.Water.Printed.RefScan';
                end
            else
                availableWater{end+1} = 'Not searched, Varian';
                preferredWater{end+1} = 'Not searched, Varian';
            end

            if isempty(availableWater)
                availableWaterList = 'None';
                preferredWaterList = 'None';
            else
                availableWaterList = strjoin(availableWater, ', ');
                preferredWaterList = strjoin(preferredWater, ', ');
            end

            % store row in Excel table
            excelStorage = [excelStorage; {
                DPid, vendor, subj, ses, rawdataJob0Met, sequenceMet, ...
                fidrefscanMet, fidrefMet, refScanYNMet, refPrintedMet, ...
                refScanNAMet, rawdataJob0Wat, sequenceWat, ...
                fidrefscanWat, fidrefWat, refScanYNWat, ...
                refPrintedWat, refScanNAWat, availableWaterList, preferredWater ...
                }];
        end       
    end
end

% create table with new columns
excelTable = cell2table(excelStorage, 'VariableNames', { ...
    'DP ID', 'Vendor','Subj','Ses', 'Metab.Acq.YN', 'Metab.Acq.Seq', 'Metab.FidRefScan','Metab.FidRef','Metab.Method.RefScan.YN', ...
    'Metab.Method.Printed.RefScan','Metab.Method.NA', 'Unsup.Water.Acq.YN', 'Unsup.Water.Acq.Seq', ...
    'Unsup.Water.FidRefScan','Unsup.Water.FidRef','Unsup.Water.Method.RefScan.YN','Unsup.Water.Printed.RefScan','Unsup.Water.Method.NA', ...
    'AvailableWater', 'PreferredWater' ...
    });

% save to Excel
excelFileName = fullfile(writeToExcelPath, '46_compMRS_waterDataCheck.xlsx');
writetable(excelTable, excelFileName, 'FileType', 'spreadsheet');
out = fprintf('Spreadsheet saved to %s\n', excelFileName);

% Available Water Data
[availableWaterData, ~, availableWaterIdx] = unique(excelTable.AvailableWater);
availableWaterCounts = accumarray(availableWaterIdx, 1);

disp('Variations of Available Water and their frequencies (in number of sessions):');
for i = 1:numel(availableWaterData)
    fprintf('%s: %d\n', string(availableWaterData(i)), availableWaterCounts(i));
end
fprintf('Total number of variations: %d\n\n', numel(availableWaterData));

% Preferred Water Data
[preferredWaterData, ~, preferredWaterIdx] = unique(excelTable.PreferredWater);
preferredWaterCounts = accumarray(preferredWaterIdx, 1);

disp('Variations of Preferred Water and their frequencies (in number of sessions):');
for i = 1:numel(preferredWaterData)
    fprintf('%s: %d\n', string(preferredWaterData(i)), preferredWaterCounts(i));
end
fprintf('Total number of variations: %d\n', numel(preferredWaterData));
