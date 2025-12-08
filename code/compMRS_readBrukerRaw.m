%compMRS_readBrukerRaw.m
%Georg Oeltzschner, Johns Hopkins University 2025
%Thanh Phong Le, EPFL, 2025
%Diana Rotaru, Medical University of Vienna, 2025
%
% USAGE:
% [fids_raw] = compMRS_readBrukerRaw(fileRaw, formatRaw)
%
% DESCRIPTION:
% This subroutine reads the complex data from 'fileRaw' in the format 
% 'formatRaw' and returns a complex time-domain FID array
%
% INPUTS:
% fileRaw     = String variable specifying the path to the scan number
%             directory.  
% formatRaw   = int (usually PV5), int32 (usually PV6 and later versions) 
%               or float64 (usually PV360 and later versions)
%
% OUTPUTS:
% fids_raw   = Complex metabolite scan data 

function fids_raw = compMRS_readBrukerRaw(fileRaw, formatRaw)

% Open and read
data     = fopen(fileRaw);
fid_data = fread(data, formatRaw);

% Make complex
real_fid = fid_data(1:2:length(fid_data));
imag_fid = fid_data(2:2:length(fid_data));
fids_raw = (real_fid+1i*imag_fid);

% Close
fclose(data);

end