%compMRS_checkBlockSize.m

% Thanh Phong Lê, CIBM Center for Biomedical Imaging and École polytechnique fédérale de Lausanne, 2026
%
% USAGE:compMRS_checkBlockSize()
%
% DESCRIPTION: Simple code to make an excel file with the calculated SNR/LW
% for data processed with separate and auto water scans, as well as the
% block size (to determine whether the output was better with or without
% frequency drift correction.
% 
% Input: none
% Output: none
 
function compMRS_checkBlockSize()
    pat = '(DP\d+).*?(sub-\d+).*?(ses-\d+)';
    filename = [datestr(datetime("now"),'yyyymmdd_HHMMss') 'LW_and_BlockSize.xlsx'];
    res = dir(['plots' filesep 'DP*.mat']);

    for ii=1:length(res)
        thisDP = erase(res(ii).name, '.mat');
        load(fullfile('plots', res(ii).name));

        subjs = size(out, 1);
        
        thisSNR_LW_ratio = zeros(subjs,1);
        thisBlock_size = zeros(subjs,1);
        thisSNR_LW_ratio_auto = zeros(subjs,1);
        thisBlock_size_auto = zeros(subjs,1);
        thisSubjSes = {};

        for jj=1:subjs
            if ~isempty(out{jj})
                thisSNR_LW_ratio(jj) = out{jj}{1}.SNR_LW_ratio;
                thisBlock_size(jj) = out{jj}{1}.block_size;
                tk = regexp(out{jj}{1}.filepath, pat, 'tokens', 'once');
                thisSubjSes{jj} = [tk{1} ' ' tk{2} ' ' tk{3}];
            
            end
            if ~isempty(out_auto{jj})
                thisSNR_LW_ratio_auto(jj) = out_auto{jj}{1}.SNR_LW_ratio;
                thisBlock_size_auto(jj) = out_auto{jj}{1}.block_size;
                tk = regexp(out_auto{jj}{1}.filepath, pat, 'tokens', 'once');
                thisSubjSes{jj} = [tk{1} ' ' tk{2} ' ' tk{3}];
            end


        end
        T = table(thisSubjSes', thisSNR_LW_ratio,thisBlock_size,thisSNR_LW_ratio_auto,thisBlock_size_auto);
        writetable(T,filename,'sheet', thisDP)
    end
end    