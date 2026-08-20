% compMRS_DPprocStats.m
% Diana Rotaru, Columbia University & Medical University of Vienna, 2025
% Mark Mikkelsen, Weill Cornell Medicine, 2025

clear; clc;

cd('A:\CoMP-MRS\CoMP-MRS-DGR\stats\ProcSpectralQuality\data\input_table\v3_20260415')

%% File names
xlsxFile       = 'CoMP-MRS_participantSpreadsheet.xlsx';
matFile        = 'CoMP-MRS-final.mat';
siteFile       = 'CoMP-MRS-sites.xlsx';
outCsv         = 'CoMP_MRS_Rstats_input.csv';
outSummaryXlsx = 'CoMP_MRS_summary_table.xlsx';

%% ------------------------------------------------------------------------
% Read Excel file
%% ------------------------------------------------------------------------
opts = detectImportOptions(xlsxFile, 'NumHeaderLines', 9);
opts.VariableNamesRange = '10:10';
opts.DataRange = '11:10000';

fullTable = readtable(xlsxFile, opts);

%% ------------------------------------------------------------------------
% Extract selected columns
%% ------------------------------------------------------------------------
% in 'CoMP-MRS_participantSpreadsheet.xlsx'
% 1  = data.packet.ID
% 3  = animal.ID
% 6  = animal.species
% 7  = animal.strain
% 9  = animal.ses-01.age
% 11 = animal.sex
% 12 = animal.ses-01.weight
% 28 = MRI.vendor
% 29 = MRI.field.strength
% 32 = MRI.coil.setup
% 35 = MRS.pulse.sequence
% 36 = MRS.brain.region
% 37 = MRS.VoI.size
% 44 = MRS.n.ave
% 62 = MRS.shim.method
selectedCols = [1 3 6 7 9 11 12 28 29 32 35 36 37 44 62];
selectedData = fullTable(:, selectedCols);

%% ------------------------------------------------------------------------
% Fill DP labels
%% ------------------------------------------------------------------------
rawDP = fullTable{:,1};

if iscell(rawDP)
    dpLabels = string(rawDP);
elseif isstring(rawDP)
    dpLabels = rawDP;
else
    dpLabels = string(rawDP);
end

dpLabels(dpLabels == "") = missing;
dpLabels = fillmissing(dpLabels, 'previous');
dpLabels = regexprep(dpLabels, '\s*\(.*?\)', '');
dpLabels = strtrim(dpLabels);

%% ------------------------------------------------------------------------
% Keep only DP01-DP36
%% ------------------------------------------------------------------------
nRows = height(fullTable);
dpNumAll = nan(nRows,1);

for i = 1:nRows
    token = regexp(dpLabels(i), 'DP\s*0*(\d+)', 'tokens', 'once');
    if ~isempty(token)
        dpNumAll(i) = str2double(token{1});
    end
end

validIdx = ~isnan(dpNumAll) & dpNumAll <= 36;

fullTable    = fullTable(validIdx,:);
selectedData = selectedData(validIdx,:);
dpLabels     = dpLabels(validIdx);

nRows = height(fullTable);

%% ------------------------------------------------------------------------
% Subject numbering within DP (internal only)
%% ------------------------------------------------------------------------
subjectNumberWithinDP = zeros(nRows,1);

currentDP = "";
countInCurrentDP = 0;

for i = 1:nRows
    if dpLabels(i) ~= currentDP
        currentDP = dpLabels(i);
        countInCurrentDP = 1;
    else
        countInCurrentDP = countInCurrentDP + 1;
    end
    subjectNumberWithinDP(i) = countInCurrentDP;
end

%% ------------------------------------------------------------------------
% Build participant table
%% ------------------------------------------------------------------------
participantTable = table();

participantTable.PacketID      = local_get_selected_column(selectedData, selectedCols, 1);
participantTable.AnimalID      = string(local_get_selected_column(selectedData, selectedCols, 3));
participantTable.AnimalSpecies = local_get_selected_column(selectedData, selectedCols, 6);
participantTable.AnimalStrain  = local_get_selected_column(selectedData, selectedCols, 7);
participantTable.AnimalAge     = local_get_selected_column(selectedData, selectedCols, 9);
participantTable.AnimalSex     = local_get_selected_column(selectedData, selectedCols, 11);
participantTable.AnimalWeight  = local_get_selected_column(selectedData, selectedCols, 12);

participantTable.MRvendor      = local_get_selected_column(selectedData, selectedCols, 28);
participantTable.MRfield       = local_get_selected_column(selectedData, selectedCols, 29);
participantTable.MRcoil        = local_get_selected_column(selectedData, selectedCols, 32);
participantTable.MRsequence    = local_get_selected_column(selectedData, selectedCols, 35);
participantTable.MRbrainregion = local_get_selected_column(selectedData, selectedCols, 36);
participantTable.MRVoxSize       = local_get_selected_column(selectedData, selectedCols, 37);
participantTable.MRaverages    = local_get_selected_column(selectedData, selectedCols, 44);

%% ------------------------------------------------------------------------
% Add extra fields from fullTable by header name
%% ------------------------------------------------------------------------
participantTable.MRsoftwareversion = local_get_fulltable_column(fullTable, "MRI.software.version");
participantTable.MRSsw             = local_get_fulltable_column(fullTable, "MRS.sw");
participantTable.MRSnpts           = local_get_fulltable_column(fullTable, "MRS.n.pts");
participantTable.MRSTE             = local_get_fulltable_column(fullTable, "MRS.TE");
participantTable.MRSTR             = local_get_fulltable_column(fullTable, "MRS.TR");
participantTable.MRSshim           = local_get_fulltable_column(fullTable, "MRS.shim.method");

%% ------------------------------------------------------------------------
% Standardize / normalize field types and values
%% ------------------------------------------------------------------------
participantTable.AnimalSpecies     = string(participantTable.AnimalSpecies);
participantTable.AnimalStrain      = string(participantTable.AnimalStrain);
participantTable.AnimalSex         = string(participantTable.AnimalSex);
participantTable.MRvendor          = string(participantTable.MRvendor);
participantTable.MRcoil            = string(participantTable.MRcoil);
participantTable.MRsequence        = string(participantTable.MRsequence);
participantTable.MRbrainregion     = string(participantTable.MRbrainregion);
participantTable.MRVoxSize           = string(participantTable.MRVoxSize);
participantTable.MRsoftwareversion = string(participantTable.MRsoftwareversion);
participantTable.MRSshim           = string(participantTable.MRSshim);

%% ------------------------------------------------------------------------
% Add Cryoprobe column (yes/no based on MRcoil)
%% ------------------------------------------------------------------------
coilStr = lower(string(participantTable.MRcoil));

isCryo = contains(coilStr, "cryo") | contains(coilStr, "cryoprobe");

Cryoprobe = repmat("FALSE", height(participantTable), 1);
Cryoprobe(isCryo) = "TRUE";

participantTable.Cryoprobe = Cryoprobe;

participantTable = movevars(participantTable, 'Cryoprobe', 'After', 'MRcoil');

% MRbrainregion
participantTable.MRbrainregion = strrep(participantTable.MRbrainregion, "Right hippocampus", "Rhippocampus");
participantTable.MRbrainregion = strrep(participantTable.MRbrainregion, "Left hippocampus", "Lhippocampus");
participantTable.MRbrainregion = strrep(participantTable.MRbrainregion, "Right striatum", "Rstriatum");
participantTable.MRbrainregion = strrep(participantTable.MRbrainregion, "Left striatum", "Lstriatum");

% MRvendor
participantTable.MRvendor = strrep(participantTable.MRvendor, "Agilent/Varian", "Varian");
participantTable.MRvendor = strrep(participantTable.MRvendor, "Varian/Agilent", "Varian");

% AnimalStrain
participantTable.AnimalStrain = strrep(participantTable.AnimalStrain, "C57BL/6J", "C57BL-6J");
participantTable.AnimalStrain = strrep(participantTable.AnimalStrain, "C57BL/6", "C57BL-6");
participantTable.AnimalStrain = strrep(participantTable.AnimalStrain, "balb/c", "BALB-c");
participantTable.AnimalStrain = strrep(participantTable.AnimalStrain, "BALB/c", "BALB-c");
participantTable.AnimalStrain = strrep(participantTable.AnimalStrain, "B6129SF2/J", "B6129SF2-J");
participantTable.AnimalStrain = strrep(participantTable.AnimalStrain, "FVB/N", "FVB-N");

% MRSshim
participantTable.MRSshim = strrep(participantTable.MRSshim, "FASTMAP/FASTESTMAP", "FASTMAP-FASTESTMAP");

%% ------------------------------------------------------------------------
% Standardize MRfield
%% ------------------------------------------------------------------------
participantTable.MRfield = str2double(string(participantTable.MRfield));
idxField = participantTable.MRfield >= 9 & participantTable.MRfield < 9.4;
participantTable.MRfield(idxField) = 9.4;

%% ------------------------------------------------------------------------
% Compute MRvoxelvolume
%% ------------------------------------------------------------------------
VoISizeProduct = nan(nRows,1);

for i = 1:nRows
    val = participantTable.MRVoxSize(i);

    if iscell(val)
        val = string(val{1});
    else
        val = string(val);
    end

    if strlength(val) == 0 || ismissing(val)
        continue
    end

    parts = regexp(val, '[xX\*\s]+', 'split');
    nums = str2double(parts);
    nums = nums(~isnan(nums));

    if ~isempty(nums)
        VoISizeProduct(i) = prod(nums);
    end
end

participantTable = addvars(participantTable, VoISizeProduct, ...
    'After', 'MRVoxSize', ...
    'NewVariableNames', 'MRvoxelvolume');


%% ------------------------------------------------------------------------
% Load MAT file
%% ------------------------------------------------------------------------
matData = load(matFile);

if ~isfield(matData, 'out')
    error('MAT file does not contain variable "out"');
end

out = matData.out;

if isfield(matData, 'out_auto')
    out_auto = matData.out_auto;
else
    out_auto = [];
end

nDPinMat = size(out,2);

%% ------------------------------------------------------------------------
% Extract LW / SNR / SNR_LW_ratio
%% ------------------------------------------------------------------------
LW = nan(nRows,1);
LW_hz = nan(nRows,1);
SNR = nan(nRows,1);
SNR_LW_ratio = nan(nRows,1);

CompCheck = repmat("pass", nRows, 1);

for i = 1:nRows
    token = regexp(dpLabels(i), 'DP\s*0*(\d+)', 'tokens', 'once');
    if isempty(token)
        continue
    end

    dpNumber = str2double(token{1});
    subjNum  = subjectNumberWithinDP(i);

    if isnan(dpNumber) || dpNumber < 1 || dpNumber > nDPinMat
        continue
    end

    % Exceptions
    if (dpNumber == 16 && subjNum == 2) || (dpNumber == 17 && subjNum == 1)
        CompCheck(i) = "fail";
        continue
    end

    matSubjNum = subjNum;
    if (dpNumber == 16) && subjNum > 2
        matSubjNum = subjNum - 1;
    end

    metricStruct = [];

    try
        dpCell = out{1, dpNumber};

        if matSubjNum <= size(dpCell,1)
            subjEntry = dpCell{matSubjNum,1};
            candidateStruct = local_extract_struct(subjEntry);

            if isstruct(candidateStruct)
                metricStruct = candidateStruct;
            end
        end
    catch
    end

    lw_val = NaN;
    lwhz_val = NaN;
    snr_val = NaN;
    ratio_val = NaN;

    if ~isempty(metricStruct)
        if isfield(metricStruct,'LW') && ~isempty(metricStruct.LW), lw_val = metricStruct.LW; end
        if isfield(metricStruct,'LW_hz') && ~isempty(metricStruct.LW_hz), lwhz_val = metricStruct.LW_hz; end
        if isfield(metricStruct,'SNR') && ~isempty(metricStruct.SNR), snr_val = metricStruct.SNR; end
        if isfield(metricStruct,'SNR_LW_ratio') && ~isempty(metricStruct.SNR_LW_ratio), ratio_val = metricStruct.SNR_LW_ratio; end
    end

    needsFallback = isnan(lw_val) || isnan(lwhz_val) || isnan(snr_val) || isnan(ratio_val);

    if needsFallback && ~isempty(out_auto) && dpNumber <= size(out_auto,2)
        try
            dpCell_auto = out_auto{1, dpNumber};

            if matSubjNum <= size(dpCell_auto,1)
                subjEntry_auto = dpCell_auto{matSubjNum,1};
                candidateStruct_auto = local_extract_struct(subjEntry_auto);

                if isstruct(candidateStruct_auto)
                    if isnan(lw_val) && isfield(candidateStruct_auto,'LW') && ~isempty(candidateStruct_auto.LW), lw_val = candidateStruct_auto.LW; end
                    if isnan(lwhz_val) && isfield(candidateStruct_auto,'LW_hz') && ~isempty(candidateStruct_auto.LW_hz), lwhz_val = candidateStruct_auto.LW_hz; end
                    if isnan(snr_val) && isfield(candidateStruct_auto,'SNR') && ~isempty(candidateStruct_auto.SNR), snr_val = candidateStruct_auto.SNR; end
                    if isnan(ratio_val) && isfield(candidateStruct_auto,'SNR_LW_ratio') && ~isempty(candidateStruct_auto.SNR_LW_ratio), ratio_val = candidateStruct_auto.SNR_LW_ratio; end
                end
            end
        catch
        end
    end

    LW(i) = lw_val;
    LW_hz(i) = lwhz_val;
    SNR(i) = snr_val;
    SNR_LW_ratio(i) = ratio_val;
end

%% ------------------------------------------------------------------------
% Add DP + metrics
%% ------------------------------------------------------------------------
participantTable.DP = dpLabels;
participantTable.SubjectInDP = subjectNumberWithinDP; % internal only
participantTable.LW = LW;
participantTable.LW_hz = LW_hz;
participantTable.SNR = SNR;
participantTable.SNR_LW_Ratio = SNR_LW_ratio;

participantTable = movevars(participantTable, {'DP','SubjectInDP'}, 'Before', 1);

%% ------------------------------------------------------------------------
% Add CompID
%% ------------------------------------------------------------------------
participantTable.CompID = strcat("compMRS", compose("%03d", (1:nRows)'));

%% ------------------------------------------------------------------------
% Add SiteID (robust DP ↔ site mapping)
%% ------------------------------------------------------------------------
siteTable = readtable(siteFile);

% Extract DP numbers from BOTH tables (robust to formatting)
siteDP_raw = string(siteTable{:,1});
participantDP_raw = string(participantTable.DP);

siteDP_num = nan(height(siteTable),1);
for i = 1:height(siteTable)
    token = regexp(siteDP_raw(i), 'DP\s*0*(\d+)', 'tokens', 'once');
    if ~isempty(token)
        siteDP_num(i) = str2double(token{1});
    end
end

participantDP_num = nan(height(participantTable),1);
for i = 1:height(participantTable)
    token = regexp(participantDP_raw(i), 'DP\s*0*(\d+)', 'tokens', 'once');
    if ~isempty(token)
        participantDP_num(i) = str2double(token{1});
    end
end

% Get site names
siteNames = string(siteTable{:,2});

% Convert site names → standardized SiteID (S01, S02, …)
[uniqueSites, ~, siteIdx] = unique(siteNames, 'stable');
standardSiteIDs = "S" + compose("%02d", (1:numel(uniqueSites))');

siteLookup = standardSiteIDs(siteIdx);

% Assign SiteID by DP number (SAFE MATCHING)
SiteID = strings(height(participantTable),1);

for i = 1:height(participantTable)
    dpNum = participantDP_num(i);
    matchIdx = find(siteDP_num == dpNum, 1);

    if ~isempty(matchIdx)
        SiteID(i) = siteLookup(matchIdx);
    else
        SiteID(i) = ""; % leave empty if no match
    end
end

participantTable.SiteID = SiteID;

%% ------------------------------------------------------------------------
% Sort rows
%% ------------------------------------------------------------------------
participantTable = sortrows(participantTable, {'DP','SubjectInDP'});

%% ------------------------------------------------------------------------
% Safe cleaning
%% ------------------------------------------------------------------------
for k = 1:width(participantTable)
    if iscell(participantTable{:,k})
        try
            participantTable.(k) = string(participantTable.(k));
        catch
        end
    end
end

%% ------------------------------------------------------------------------
% Add CompCheck
%% ------------------------------------------------------------------------
participantTable = addvars(participantTable, CompCheck, ...
    'NewVariableNames', 'CompCheck');

%% ------------------------------------------------------------------------
% Create summary table
% PI, location, n_fail, comments removed
%% ------------------------------------------------------------------------
uniqueDPs = unique(string(participantTable.DP), 'stable');
nDP = numel(uniqueDPs);

summaryTable = table();
summaryTable.dataset = uniqueDPs;
summaryTable.site = strings(nDP,1);
summaryTable.species = strings(nDP,1);
summaryTable.animal_strain = strings(nDP,1);
summaryTable.age_mean = nan(nDP,1);
summaryTable.age_SD = nan(nDP,1);
summaryTable.sex = strings(nDP,1);
summaryTable.vendor = strings(nDP,1);
summaryTable.field_strength = nan(nDP,1);
summaryTable.software_version = strings(nDP,1);
summaryTable.coil_setup = strings(nDP,1);
summaryTable.sequence = strings(nDP,1);
summaryTable.sw = nan(nDP,1);
summaryTable.n_pts = nan(nDP,1);
summaryTable.TE = nan(nDP,1);
summaryTable.TR = nan(nDP,1);
summaryTable.shim_method = strings(nDP,1);
summaryTable.number_of_averages = nan(nDP,1);
summaryTable.region = strings(nDP,1);
summaryTable.voxel_size = strings(nDP,1);
summaryTable.voxel_volume = nan(nDP,1);
summaryTable.LW_mean = nan(nDP,1);
summaryTable.LW_hz_mean = nan(nDP,1);
summaryTable.SNR_mean = nan(nDP,1);
summaryTable.SNR_LW_ratio_mean = nan(nDP,1);

for d = 1:nDP
    dp = uniqueDPs(d);
    subT = participantTable(string(participantTable.DP) == dp,:);

    summaryTable.site(d) = local_first(subT.SiteID);
    summaryTable.species(d) = local_first(subT.AnimalSpecies);
    summaryTable.animal_strain(d) = local_first(subT.AnimalStrain);
    summaryTable.vendor(d) = local_first(subT.MRvendor);
    summaryTable.software_version(d) = local_first(subT.MRsoftwareversion);
    summaryTable.coil_setup(d) = local_first(subT.MRcoil);
    summaryTable.sequence(d) = local_first(subT.MRsequence);
    summaryTable.shim_method(d) = local_first(subT.MRSshim);
    summaryTable.region(d) = local_first(subT.MRbrainregion);
    summaryTable.voxel_size(d) = local_first(subT.MRVoxSize);

    ageVals = str2double(string(subT.AnimalAge));
    ageVals = ageVals(~isnan(ageVals));
    if ~isempty(ageVals)
        summaryTable.age_mean(d) = mean(ageVals);
        if numel(ageVals) > 1
            summaryTable.age_SD(d) = std(ageVals,0);
        else
            summaryTable.age_SD(d) = 0;
        end
    end

    sexVals = upper(strtrim(string(subT.AnimalSex)));
    nF = sum(sexVals == "F");
    nM = sum(sexVals == "M");
    if nF > 0 && nM > 0
        summaryTable.sex(d) = sprintf('%dF, %dM', nF, nM);
    elseif nF > 0
        summaryTable.sex(d) = sprintf('%dF', nF);
    elseif nM > 0
        summaryTable.sex(d) = sprintf('%dM', nM);
    else
        summaryTable.sex(d) = "";
    end

    fieldVals = subT.MRfield;
    fieldVals = fieldVals(~isnan(fieldVals));
    if ~isempty(fieldVals), summaryTable.field_strength(d) = fieldVals(1); end

    vvVals = subT.MRvoxelvolume;
    vvVals = vvVals(~isnan(vvVals));
    if ~isempty(vvVals), summaryTable.voxel_volume(d) = vvVals(1); end

    lwVals = subT.LW;
    lwVals = lwVals(~isnan(lwVals));
    if ~isempty(lwVals), summaryTable.LW_mean(d) = mean(lwVals); end

    lwhzVals = subT.LW_hz;
    lwhzVals = lwhzVals(~isnan(lwhzVals));
    if ~isempty(lwhzVals), summaryTable.LW_hz_mean(d) = mean(lwhzVals); end

    snrVals = subT.SNR;
    snrVals = snrVals(~isnan(snrVals));
    if ~isempty(snrVals), summaryTable.SNR_mean(d) = mean(snrVals); end

    ratioVals = subT.SNR_LW_Ratio;
    ratioVals = ratioVals(~isnan(ratioVals));
    if ~isempty(ratioVals), summaryTable.SNR_LW_ratio_mean(d) = mean(ratioVals); end

    swVals = str2double(string(subT.MRSsw));
    swVals = swVals(~isnan(swVals));
    if ~isempty(swVals), summaryTable.sw(d) = swVals(1); end

    nptsVals = str2double(string(subT.MRSnpts));
    nptsVals = nptsVals(~isnan(nptsVals));
    if ~isempty(nptsVals), summaryTable.n_pts(d) = nptsVals(1); end

    teVals = str2double(string(subT.MRSTE));
    teVals = teVals(~isnan(teVals));
    if ~isempty(teVals), summaryTable.TE(d) = teVals(1); end

    trVals = str2double(string(subT.MRSTR));
    trVals = trVals(~isnan(trVals));
    if ~isempty(trVals), summaryTable.TR(d) = trVals(1); end

    navVals = str2double(string(subT.MRaverages));
    navVals = navVals(~isnan(navVals));
    if ~isempty(navVals), summaryTable.number_of_averages(d) = navVals(1); end
end

%% ------------------------------------------------------------------------
% Remove columns from final subject-level output
%% ------------------------------------------------------------------------
participantTable.PacketID = [];
participantTable.SubjectInDP = [];

%% ------------------------------------------------------------------------
% Put CompID and SiteID first
%% ------------------------------------------------------------------------
participantTable = movevars(participantTable, {'CompID','SiteID'}, 'Before', 1);

%% ------------------------------------------------------------------------
% Format summary table numeric columns to 2 decimals
%% ------------------------------------------------------------------------
for k = 1:width(summaryTable)
    if isnumeric(summaryTable{:,k})
        summaryTable{:,k} = round(summaryTable{:,k}, 2);
    end
end

%% ------------------------------------------------------------------------
% Reorder FULL summary table columns (core first, rest after)
%% ------------------------------------------------------------------------

coreCols = { ...
    'dataset','site','species','strain','sex','vendor','B0',...
    'software_version','sequence','coil','shim','averages',...
    'region','voxel_size','voxel_volume','sw','TE','TR','n_points' ...
};

% Keep only those that actually exist (safety)
coreCols = coreCols(ismember(coreCols, summaryTable.Properties.VariableNames));

% Get remaining columns automatically
otherCols = setdiff(summaryTable.Properties.VariableNames, coreCols, 'stable');

% Final order: core first, then everything else
summaryTable = summaryTable(:, [coreCols, otherCols]);

%% ------------------------------------------------------------------------
% Create SHORT summary table (selected columns only)
%% ------------------------------------------------------------------------
%% This section of the code needs updates to match the formats and entries#
%% in the table in figure 4 (boxplots)
summaryShort = table();

summaryShort.dataset           = summaryTable.dataset;
summaryShort.site              = summaryTable.site;
summaryShort.species           = summaryTable.species;
summaryShort.animal_strain     = summaryTable.animal_strain;
summaryShort.sex               = summaryTable.sex;
summaryShort.vendor            = summaryTable.vendor;

% Rename field_strength → B0
summaryShort.B0                = summaryTable.field_strength;

summaryShort.software_version  = summaryTable.software_version;
summaryShort.coil_setup        = summaryTable.coil_setup;
summaryShort.sequence          = summaryTable.sequence;
summaryShort.sw                = summaryTable.sw;
summaryShort.n_pts             = summaryTable.n_pts;
summaryShort.TE                = summaryTable.TE;
summaryShort.TR                = summaryTable.TR;
summaryShort.shim_method       = summaryTable.shim_method;

% Rename number_of_averages → averages
summaryShort.averages          = summaryTable.number_of_averages;

summaryShort.region            = summaryTable.region;
summaryShort.voxel_size        = summaryTable.voxel_size;

%% ------------------------------------------------------------------------
% Export
%% ------------------------------------------------------------------------
writetable(participantTable, outCsv);
writetable(summaryTable, outSummaryXlsx, 'Sheet','Summary');

writetable(summaryShort, outSummaryXlsx, 'Sheet','Summary_short');

disp('Done.')

%% ------------------------------------------------------------------------
% Helpers
%% ------------------------------------------------------------------------
function col = local_get_selected_column(selectedData, selectedCols, originalColNumber)
idx = find(selectedCols == originalColNumber, 1);
if isempty(idx)
    error('Requested original column %d is not present in selectedCols.', originalColNumber);
end
col = selectedData{:, idx};
end

function col = local_get_fulltable_column(tbl, candidateName)
varNames = string(tbl.Properties.VariableNames);
normVars = local_normalize_names(varNames);
normTarget = local_normalize_names(string(candidateName));

idx = find(normVars == normTarget, 1);
if isempty(idx)
    error('Could not find column "%s" in fullTable.', string(candidateName));
end

col = tbl{:,idx};
end

function out = local_normalize_names(x)
out = lower(string(x));
out = regexprep(out, '[^a-z0-9]', '');
end

function s = local_extract_struct(entry)
s = [];

if isstruct(entry)
    s = entry;
    return
end

if iscell(entry)
    if isempty(entry)
        return
    end

    firstVal = entry{1};

    if isstruct(firstVal)
        s = firstVal;
        return
    end

    if iscell(firstVal) && ~isempty(firstVal)
        innerVal = firstVal{1};
        if isstruct(innerVal)
            s = innerVal;
            return
        end
    end
end
end

function val = local_first(x)
x = string(x);
x = x(x ~= "" & ~ismissing(x));
if isempty(x)
    val = "";
else
    val = x(1);
end
end