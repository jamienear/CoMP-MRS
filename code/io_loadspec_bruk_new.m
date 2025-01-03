%io_loadspec_bruk.m
%Jamie Near, Sunnybrook Research Institute, 2023
%Wendy Oakden, Sunnybrook Research Institute, 2020
%Chathura Kumaragamage, McGill University 2016.
%Georg Oeltzschner, Johns Hopkins University 2025
%
% USAGE:
% [out,ref]=io_loadspec_bruk_new(inDir, rawData);
%
% DESCRIPTION:
% Reads in Bruker MRS data (fid.raw, fid.ref).
%
% io_loadspec_bruk_new reads a Bruker data file and outputs a data structure in
% FID-A format (with fields corresponding to time scale, fids, frequency
% scale, spectra, and header fields containing information about the
% acquisition.  The resulting matlab structure can be operated on by the
% other FID-A functions.
%
% This function gives you the option to read in either the raw data (with
% individual averages stored separately), or the combined data.
%
% For Bruker versions PV5 and earlier, the water suppressed and water
% unsuppressed data were acquired separately.  In such cases,
% io_loadspec_bruk would need to be run separately on the water suppressed
% and water unsuppressed scan folders, and the 2nd output 'ref' refers to
% any separately acquired navigator echoes.  For Bruker versions PV6 and later,
% it became common for the water suppressed data to be acquired
% concurrently with the water unsuppressed data (within a single scan ID).
% For such datasets, this function reads both the water suppressed data
% (1st output, "out"), and the water unsuppressed data (2nd output
% "ref").
%
% This function automatically reads the numbers of points before the echo
% from the $GRPDLY header field.
%
% INPUTS:
% inDir     = String variable specifying the path to the scan number
%             directory.  Alternatively, inDir can be an integer specifying
%             number of the scan directory to analyze (assuming it is in the
%             current directory).
% rawData   = 'y' or 'n' (default = 'y')
%             'y' - load the individually stored fids
%             'n' - load the Bruker combined data
%
% OUTPUTS:
% out = Input dataset in FID-A structure format.
% ref = The Reference scan data (Navigator echoes in PV5, water reference
%       data in later versions) in FID-A structure format, if applicable.

function [out,ref]=io_loadspec_bruk_new(inDir, rawData)

%Allow the user to pass the input directory as either a string or an
%integer.
if isnumeric(inDir)
    if ~rem(inDir,1)==0
        error('ERROR: only integer values of ''inDir'' are allowed.  ABORTING');
    else
        inDir=num2str(inDir);
    end
end

% Populate the header information from the ACQP file
acqpFile        = fullfile(inDir, 'acqp');
headerACQP      = parseBrukerFormat(acqpFile);

% Populate the header information from the Method file
methodFile      = fullfile(inDir, 'method');
headerMethod    = parseBrukerFormat(methodFile);

% Populate the header information from the ACQUS file
acqusFile      = fullfile(inDir, 'acqus');
if isfile(acqusFile)
    headerACQUS    = parseBrukerFormat(acqusFile);
end

% Get a few important bits
version         = headerACQP.ACQ_sw_version;
if contains(version, ["PV 5", "PV 6", "PV 7", "PV-7"])
    rawDataPoints = headerMethod.PVM_DigNp;
elseif contains(version,'PV-360')
    rawDataPoints = headerMethod.PVM_SpecMatrix;
end
rawAverages     = headerMethod.PVM_NAverages;
Nrcvrs          = headerMethod.PVM_EncNReceivers;
if Nrcvrs>1
    multiRcvrs=true;
else
    multiRcvrs=false;
end
if contains(version, ["PV 5", "PV 6", "PV 7", "PV-7"])
    leftshift       = round(headerACQP.GRPDLY);
elseif contains(version,'PV-360')
    leftshift       = round(headerACQUS.GRPDLY);
end

if contains(version, ["PV 5", "PV 6", "PV 7", "PV-7"])
    spectralwidth = headerMethod.PVM_DigSw;
elseif contains(version,'PV-360')
    spectralwidth = headerMethod.PVM_SpecSWH;
end
txfrq = headerACQP.BF1*1e6;
Bo=txfrq/42577000;
spectralwidthppm=spectralwidth/(txfrq/1e6);
te = headerMethod.PVM_EchoTime;
tr = headerMethod.PVM_RepetitionTime;
sequence = headerMethod.Method;


%Now load in the data.  Either raw or averaged data, depending on what was
%requested via the 'rawData' parameter. 
if strcmp(rawData,'y') || strcmp(rawData,'Y')
    if contains(version, ["PV 6", "PV 7", "PV-7", "PV-360"])
        data = fopen(fullfile(inDir, 'rawdata.job0')); % WO - changed for PV6.0
        fid_data = fread(data,'int32');
    elseif contains(version,'PV 5')
        data = fopen(fullfile(inDir, 'fid.raw'));
        fid_data = fread(data,'int');
    end
    real_fid = fid_data(2:2:length(fid_data));
    imag_fid = fid_data(1:2:length(fid_data));
    fids_raw=real_fid-1i*imag_fid;
    fclose(data); %WO - added fclose

    averages=rawAverages;  %since these data are uncombined;
    out.flags.averaged=0; %make the flags structure

    %Reshape to put the averages along a 2nd dimension
    fids_raw=reshape(fids_raw,[],rawAverages);

    %If there are multiple receivers *I think* that these always get stored
    %separately by default in the fid.raw file.  Therefore, at this stage, 
    %if this is a PV360 dataset with multiple receivers, we need to reshape 
    %the dataset again:
    if ~contains(version,'PV 5') && multiRcvrs
        fids_raw=reshape(fids_raw,rawDataPoints,Nrcvrs,rawAverages);
        %Permute so that time is along 1st dimension, averages is along 2nd 
        %dimension, and coils is along 3rd dimension:
        fids_raw=permute(fids_raw,[1,3,2]); 
    end

elseif strcmp(rawData,'n') || strcmp(rawData,'N')
    %REQUEST PROCESSED DATA ONLY:  Use the FID file instead of fid.raw.
    if contains(version, ["PV-360"])
        data = fopen(fullfile(inDir, 'pdata', '1', 'fid_proc.64'));
        fid_data=fread(data,'float64');
    else
        data = fopen(fullfile(inDir, 'fid'));
        fid_data=fread(data,'int');
    end
    real_fid = fid_data(2:2:length(fid_data));
    imag_fid = fid_data(1:2:length(fid_data));
    fids_raw=real_fid-1i*imag_fid;
    fclose(data); %WO

    averages=1; %since we requested the combined data.
    out.flags.averaged=1; %make the flags structure


    %***JN***
    % REMOVED THIS SECTION.  DON'T THINK IT IS NECESSARY.
    %     if mod(size(real_fid,1),rawDataPoints) ~=0
    %         display 'number of repetitions cannot be accurately found';
    %     end
    %     fids_raw=reshape(fids_raw,rawDataPoints,[]);
    %***JN***
else
    error('ERROR:  rawData variable not recognized.  Options are ''y'' or ''n''.');
end

%Perform left-shifting to remove points before the echo
fids_trunc=fids_raw(leftshift+1:end,:,:);

%replace the left-shifted points with zeros at the end
fids=padarray(fids_trunc, [leftshift,0],'post');

%Do the fourier transform
specs=fftshift(ifft(fids,[],1),1);
sz=size(specs); %size of the array

%specify the dims:
%Time dimension:
dims.t=1;%the time dimension is always the 1st dimension

%Coils dimension:
%As far as I am aware, the RF coils are only stored separately in PV360.
if ~contains(version,'PV 5') && multiRcvrs
    %Coils dimension should normally be after the averages dimension,
    %unless there are no averages, in which case the coild dimension will
    %be after the time dimension.  
    if rawAverages==1
        dims.coils=2;
    elseif rawAverages>1
        dims.coils=3;
    end
else
    dims.coils=0;
end

if strcmp(rawData,'y') || strcmp(rawData,'Y')
    dims.averages=2;
elseif strcmp(rawData,'n') || strcmp(rawData,'N')
    dims.averages=0;
end
%I have not encountered any Bruker datasets so far where there are
%subSpectra or extras dimensions.  
dims.subSpecs=0;
dims.extras=0;


%NOW TRY LOADING IN THE REFERENCE SCAN DATA (IF IT EXISTS)
isRef=false;
if contains(version,'PV 5')
    isRef=exist(fullfile(inDir, 'fid.ref'));
elseif contains(version,'PV 6') || contains(version,'PV 7') || contains(version,'PV-360.2')
    isRef=exist(fullfile(inDir, 'fid.refscan'));
elseif contains(version,'PV-360.3')
    isRef=exist(fullfile(inDir, 'pdata', '1', 'fid_refscan.64'));
end

if isRef
    if contains(version,'PV 5')
        data = fopen(fullfile(inDir, 'fid.ref'));
        ref_data=fread(data,'int16');
    elseif contains(version,'PV 6') || contains(version,'PV 7') || contains(version,'PV-360.2')
        data = fopen(fullfile(inDir, 'fid.refscan'));
        ref_data=fread(data,'int32');
    elseif contains(version,'PV-360.3')
        data = fopen(fullfile(inDir, 'pdata', '1', 'fid_refscan.64'));
        ref_data=fread(data,'float64');
    end
    real_ref = ref_data(2:2:length(ref_data));
    imag_ref = ref_data(1:2:length(ref_data));
    fids_ref=real_ref-1i*imag_ref;
    fclose(data); %WO - added fclose

    %Find the number of averages in the ref dataset:
    if contains(version,'PV 5')
        rawAverages_ref=rawAverages;
        averages_ref=rawAverages_ref;
    else
        rawAverages_ref = headerMethod.PVM_RefScanNA;
        averages_ref=rawAverages_ref;
    end

    fids_ref=reshape(fids_ref,[],rawAverages_ref);
    fids_ref_trunc=fids_ref(leftshift+1:end,:);
    reffids=padarray(fids_ref_trunc, [leftshift,0],'post');
    refspecs=fftshift(ifft(reffids,[],1),1);
    sz_ref=size(refspecs);
    ref.flags.averaged=0;
    %specify the dims
    refdims.t=1;
    refdims.coils=0;
    if rawAverages_ref>1
        refdims.averages=2;
    else
        refdims.averages=0;
    end    
    refdims.subSpecs=0;
    refdims.extras=0;
else
    %Ref scans not found.  Print warning:
    disp('WARNING REFERENCE SCANS NOT FOUND.  RETURNING EMPTY REF STRUCTURE.');
end

%Specify the number of subspecs.  For now, this will always be one.
subspecs=1;
rawSubspecs=1;

%calculate the ppm scale
ppm=[4.65+(spectralwidthppm/2):-spectralwidthppm/(length(specs)-1):4.65-(spectralwidthppm/2)];

%calculate the dwelltime:
dwelltime=1/spectralwidth;

%calculate the time scale
t=[0:dwelltime:(sz(1)-1)*dwelltime];


%FILLING IN DATA STRUCTURE FOR THE FID.RAW DATA
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;  
out.t=t;    
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
out.txfrq=txfrq;
out.date=date;
out.dims=dims;
out.Bo=Bo;
out.averages=averages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=sequence;
out.te=te;
out.tr=tr;
out.pointsToLeftshift=0;
out.version=version;


%FILLING IN THE FLAGS FOR THE FID.RAW DATA
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=1;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;

if multiRcvrs && dims.coils
    out.flags.addedrcvrs=0;
else
    out.flags.addedrcvrs=1;
end
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.avgNormalized=0;
if out.dims.subSpecs==0
    out.flags.isFourSteps=0;
else
    out.flags.isFourSteps=(out.sz(out.dims.subSpecs)==4);
end


if isRef
    %FILLING IN DATA STRUCTURE FOR THE FID.REF DATA
    ref.fids=reffids;
    ref.specs=refspecs;
    ref.sz=sz_ref;
    ref.ppm=ppm;
    ref.t=t;
    ref.spectralwidth=spectralwidth;
    ref.dwelltime=dwelltime;
    ref.txfrq=txfrq;
    ref.date=date;
    ref.dims=refdims;
    ref.Bo=Bo;
    ref.averages=averages_ref;
    ref.rawAverages=rawAverages_ref;
    ref.subspecs=subspecs;
    ref.rawSubspecs=rawSubspecs;
    ref.seq=sequence;
    ref.te=te;
    ref.tr=tr;
    ref.pointsToLeftshift=0;
    ref.version=version;
    
    
    %FILLING IN THE FLAGS FOR THE FID.REF DATA
    ref.flags.writtentostruct=1;
    ref.flags.gotparams=1;
    ref.flags.leftshifted=1;
    ref.flags.filtered=0;
    ref.flags.zeropadded=0;
    ref.flags.freqcorrected=0;
    ref.flags.phasecorrected=0;
    
    ref.flags.addedrcvrs=1;
    ref.flags.subtracted=0;
    ref.flags.writtentotext=0;
    ref.flags.downsampled=0;
    ref.flags.avgNormalized=0;
    if ref.dims.subSpecs==0
        ref.flags.isFourSteps=0;
    else
        ref.flags.isFourSteps=(ref.sz(ref.dims.subSpecs)==4);
    end
else
    %REF NOT FOUND.  RETURN EMPTY STRUCTURE.
    ref = [];
end


end




function header = parseBrukerFormat(inputFile)
% This subroutine uses regular expressions and case differentiations to
% extract all relevant information from a Bruker-formatted header file
% (acqp, method, etc.)

% Open file
fid = fopen(inputFile);

% Get first line
tline = fgets(fid);

% Loop over subsequent lines
while ~feof(fid)

    % First, get the parameters without a $
    [tokens, matches] = regexp(tline,'##([\w\[\].]*)\s*=\s*([-\(\w\s.\"\\:\.,\)]*)','tokens','match');

    % When a matching string is found, parse the results into a struct
    if length(tokens) == 1

        fieldname = regexprep(tokens{1}{1}, '\[|\]',''); % delete invalid characters

        % Convert numbers to doubles, leave strings & empty lines alone
        if ~isnan(str2double(tokens{1}{2}))
            value = str2double(tokens{1}{2});
        else
            value = strtrim(tokens{1}{2});
        end

        % Convert char to string
        if ischar(value)
            value = string(value);
        end

        % Store
        header.(fieldname) = value;

        % Get next line
        tline = fgets(fid);
        continue

    else

        % If not a match, get the parameters with a $
        [tokens, ~] = regexp(tline,'##\$([\w\[\].]*)\s*=\s*([-\(\w\s.\"\\:\.,\)]*)','tokens','match');


        % When a matching string is found, parse the results into a struct
        if length(tokens) == 1

            fieldname = regexprep(tokens{1}{1}, '\[|\]',''); % delete invalid characters

            % Determine if the value indexes an array (signaled by a number
            % inside a double bracket, e.g. ##$PULPROG=( 32 )), or a single
            % value (signaled by just a string, e.g. ##$ACQ_user_filter_mode=Special)
            [tokensValue, ~] = regexp(tokens{1}{2},'\( (.*) \)','tokens','match');

            % If there's a match, we need to parse the subsequent lines
            % which contain the array
            if length(tokensValue) == 1

                % Arrays can span more than one line, unfortunately, so we
                % need to do some clever pattern matching - basically, we
                % want to extract lines until they lead with ## or $$
                % again:
                endOfBlock = 0;
                multiLine  = '';
                while endOfBlock ~=1

                    % Get next line
                    tline = fgets(fid);

                    if contains(string(tline), ["$$", "##"])
                        endOfBlock = 1;
                    else
                        multiLine = [multiLine, tline];
                    end

                end

                % If the line is bracketed by <>, store that contents as
                % one
                contents = {};
                [tokensBrackets, ~] = regexp(multiLine,'<(.*)>\n','tokens','match');
                if length(tokensBrackets) == 1
                    contents{1} = tokensBrackets{1}{1};

                    % Convert numbers to doubles, leave strings & empty lines alone
                    if ~isnan(str2double(contents{1}))
                        value = str2double(contents{1});
                    else
                        value = strtrim(contents{1});
                    end

                    % Convert char to string
                    if ischar(value)
                        value = string(value);
                    end

                else
                    % If not, it's an array.
                    % Sometimes this array can even contain vectors, for example TPQQ
                    % In this case, let's look for recurring parentheses
                    % again:
                    multiLine = erase(multiLine, newline); % remove new line characters
                    [tokensParentheses, ~] = regexp(multiLine,'\(([^\)]+)\)','tokens','match');

                    if ~isempty(tokensParentheses)
                        for rr = 1:length(tokensParentheses)
                            test = textscan(tokensParentheses{1,rr}{1}, '%s', 'Delimiter', ',');
                            contents{rr} = test{1};
                        end
                    else
                        % use textscan to convert space-delimited vectors to cell array
                        test = textscan(multiLine, '%s');
                        contents = test{1};
                    end

                    % Convert numbers to doubles, leave strings & empty lines alone
                    if ~isnan(str2double(contents))
                        value = str2double(contents);
                    else
                        value = strtrim(contents);
                    end
                    
                    % Convert char to string
                    if ischar(value)
                        value = string(value);
                    end

                    if iscell(value) && length(value) == 1
                        value = value{1};
                    end

                end

                % Store
                header.(fieldname) = value;

                continue

            else

                % Convert numbers to doubles, leave strings & empty lines alone
                if ~isnan(str2double(tokens{1}{2}))
                    value = str2double(tokens{1}{2});
                else
                    value = strtrim(tokens{1}{2});
                end

                % Convert char to string
                if ischar(value)
                    value = string(value);
                end

                if iscell(value) && length(value) == 1
                    value = value{1};
                end

                % Store
                header.(fieldname) = value;

            end

            % Get next line
            tline = fgets(fid);
            continue

        else

            % Get next line
            tline = fgets(fid);
            continue

        end

    end

end

fclose(fid);

end
