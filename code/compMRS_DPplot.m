% compMRS_DPplot.m
% Jamie Near, Sunnybrook Research Institute, 2026
% Diana Rotaru, Medical University of Vienna, 2026
%
% USAGE:
% [out_as, mu, sigma]=compMRS_DPplot(inDatastruct,DPnum);
%
% DESCRIPTION:
% This function accepts a cell array of processed data from the CoMP-MRS
% project and returns a plot showing the mean spectrum with shaded region
% corresponding to +/- SD, across all subjects within the DP. Before
% computing the average and SD, the spectra are scaled and aligned.
%
% INPUTS:
% inDatastruct: path to .mat file
% DPnum: integer ID of the DP to be plotted
%
% OUTPUTS:
% out_as: processed, aligned spectra
% mu: mean spectrum
% sigma: standard deviation spectrum

function [out_as, mu, sigma] = compMRS_DPplot(inDatastruct, DPnum)

out = cell(1,8);
out_sc = cell(1,8);

S = load(inDatastruct);

hasOut = isfield(S,'out') && ~isempty(S.out);
hasOutAuto = isfield(S,'out_auto') && ~isempty(S.out_auto);

if ~hasOut && ~hasOutAuto
    error('Neither out nor out_auto exists in the input file.');
end

% First extract the desired DP scans from the input structure
for n = 1:8
    tmp = [];

    % First try standard layout: out{1,DPnum}{n,1}{1,1}
    if hasOut && DPnum <= size(S.out,2) && ~isempty(S.out{1,DPnum})
        if n <= size(S.out{1,DPnum},1) && ~isempty(S.out{1,DPnum}{n,1})
            candidate = S.out{1,DPnum}{n,1};

            if iscell(candidate)
                if ~isempty(candidate{1,1})
                    tmp = candidate{1,1};
                end
            elseif isstruct(candidate)
                tmp = candidate;
            end
        end
    end

    % Fallback for DP 4 / 9 or any empty standard entry:
    % out_auto{1,DPnum}{n,1}
    if isempty(tmp) && hasOutAuto && DPnum <= size(S.out_auto,2) && ~isempty(S.out_auto{1,DPnum})
        if n <= size(S.out_auto{1,DPnum},1) && ~isempty(S.out_auto{1,DPnum}{n,1})
            candidate = S.out_auto{1,DPnum}{n,1};

            if iscell(candidate)
                if ~isempty(candidate{1,1})
                    tmp = candidate{1,1};
                end
            elseif isstruct(candidate)
                tmp = candidate;
            end
        end
    end

    if isempty(tmp)
        out{n} = [];
        continue;
    end

    if ~isfield(tmp,'flags') || ~isstruct(tmp.flags)
        tmp.flags = struct();
    end
    tmp.flags.averaged = 1;

    out{n} = tmp;
end

% Scale each spectrum so that the Creatine peak has amplitude 1.0
for n = 1:8
    if isempty(out{n})
        out_sc{n} = [];
        continue;
    end

    out_sc{n} = op_ppmref(out{n},1.8,2.2,2.01);
    pk = op_getPeakHeight(out_sc{n},2.9,3.1);

    if isempty(pk) || ~isfinite(pk) || pk == 0
        out_sc{n} = [];
        continue;
    end

    out_sc{n} = op_ampScale(out_sc{n},1/pk);
end

out_sc = out_sc(~cellfun(@isempty,out_sc));

if isempty(out_sc)
    out_as = {};
    mu = [];
    sigma = [];
    return;
end

% Align all scans in the DP
out_as = op_alignAllScans_fd(out_sc,0.2,4.2,0.1);

if isempty(out_as)
    mu = [];
    sigma = [];
    return;
end

% Mean +/- SD plot values
[~, mu, sigma] = op_plotMeanSD(out_as,0.2,5.5,'','','',1.96);

end