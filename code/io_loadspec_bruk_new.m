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
% and water unsuppressed scan folders. When run on the water-suppressed 
% folders, the 2nd output 'ref' will be empty, and the 3rd output 'nav'
% will contain the acquired navigator echoes.  
% 
% For Bruker versions PV6 and later,
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
% out        = Metabolite scan data in FID-A structure format.
% ref        = Reference scan data in FID-A structure format, if applicable.
% nav        = Navigator echo data in FID-A structure format, if applicable.
% coilcombos = Structure containing two fields:
%              ph:  Vector of coil phases (in [degrees]) used for alignment.
%              sig: Vector of coil weights.

function [out,ref,nav,coilcombos]=io_loadspec_bruk_new(inDir, rawData)

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
% Software version
version         = headerACQP.ACQ_sw_version;

% Number of spectral points
if contains(version, ["PV 5", "PV 6", "PV 7", "PV-7"])
    rawDataPoints = headerMethod.PVM_DigNp;
elseif contains(version,'PV-360')
    rawDataPoints = headerMethod.PVM_SpecMatrix;
end

% Number of transients
rawAverages     = headerMethod.PVM_NAverages;

% Number of repetitions
if isfield(headerMethod, 'PVM_NRepetitions')
    rawRepetitions  = headerMethod.PVM_NRepetitions;
else
    rawRepetitions = 1;
end

% Number of receivers used
Nrcvrs          = headerMethod.PVM_EncNReceivers;
if Nrcvrs>1
    multiRcvrs = true;
end

% A few other fields seem to be relevant here
if isfield(headerMethod, 'PVM_ArrayPhase')
    arrayPhase = headerMethod.PVM_ArrayPhase;
else
    arrayPhase = zeros(Nrcvrs,1);
end
% At least one dataset has only one channel but stores the correct
% phase in headerMethod.PVM_EncChanPhase
if isfield(headerMethod, 'PVM_EncChanPhase')
    encChanPhase = headerMethod.PVM_EncChanPhase;
else
    encChanPhase = zeros(Nrcvrs,1);
end
% At least one dataset (PV5) has what appears to be channel scaling coefficients
% stored, even if PV5 does not store separate channel data.
if isfield(headerMethod, 'PVM_EncChanScaling')
    encChanScaling = headerMethod.PVM_EncChanScaling;
else
    encChanScaling = ones(size(Nrcvrs));
end

if Nrcvrs>1
    % If multi-phase, get coil phases
    multiRcvrs=true;
    if isfield(headerMethod, 'PVM_ArrayPhase')
        rcvrPhases = headerMethod.PVM_ArrayPhase;
    else
        rcvrPhases = zeros(Nrcvrs,1);
    end
else
    multiRcvrs=false;
    
    if isfield(headerMethod, 'PVM_EncChanPhase')
        rcvrPhases = headerMethod.PVM_EncChanPhase;
    else
        rcvrPhases = zeros(Nrcvrs,1);
    end
end

% Group delay (i.e. left shift)
if contains(version, ["PV 5", "PV 6", "PV 7", "PV-7"])
    leftshift       = round(headerACQP.GRPDLY);
elseif contains(version,'PV-360')
    % I'm not sure this differentiation is necessary but Jessie is getting
    % GRPDELY from the ACQUS file - might remove if not needed
    if exist('headerACQUS', 'var')
        leftshift       = round(headerACQUS.GRPDLY);
    else
        leftshift       = round(headerACQP.GRPDLY);
    end
    
    % Found a couple instances of GRPDLY being -1, in that case, set it to
    % the PV360 default which is 77 (from Brayan Alves' script)
    if leftshift == -1
        leftshift = 77;
    end
end

% Spectral width
if contains(version, ["PV 5", "PV 6", "PV 7", "PV-7"])
    spectralwidth = headerMethod.PVM_DigSw;
elseif contains(version,'PV-360')
    spectralwidth = headerMethod.PVM_SpecSWH;
end

% TX freqs for PV6 onwards are explicitly stored separately for metabolite
% and ref scans
if contains(version, "PV 5")
    txfrq       = headerACQP.BF1*1e6;
    txfrq_ref   = txfrq;
elseif contains(version, ["PV 6", "PV 7", "PV-7", "PV-360"])
    txfrq       = headerMethod.PVM_FrqWork(1)*1e6;
    txfrq_ref   = headerMethod.PVM_FrqRef(1)*1e6;
end
Bo=txfrq/42577000;

% spectralwidthppm is stored in the header BUT it is slightly different 
% than if we calculate it from txfrq and spectralwidth - this suggests we 
% may have to correct txfrq itself with one of the freq offsets?
% Note: We are now calculating the ppm axis via a frequency axis (as in the
% other FID-A functions), so this variable is probably not even needed
% anymore - keeping it for now for sanity checks)
spectralwidthppm=spectralwidth/(txfrq_ref/1e6); % old method
spectralwidthppm=headerMethod.PVM_SpecSW; % better method?

% Center frequency is *explicitly* stored in the headers from PV6 onwards -
% mostly it's been 4.7 so we'll assume that for PV5 too. 
% For compatibility with the other FID-A functions, we should probably
% store centerfreq in the FID-A header?
if contains(version, ["PV 5", "PV 6", "PV 7", "PV-7"])
    centerfreq      = 4.7;
    centerfreq_ref  = centerfreq;
elseif contains(version,'PV-360')
    centerfreq      = headerMethod.PVM_FrqRefPpm(1);
end

% Receiver gain which we seem to need when reconciling raw and processed data
% According to the manual, this appears to be stored in ACQ_jobs (in fifth
% position), *NOT* the regular ACQP RG parameter.
if isfield(headerACQP, 'ACQ_jobs')
    if iscell(headerACQP.ACQ_jobs{1})
        receiverGain = str2double(headerACQP.ACQ_jobs{1}{5});
    else
        receiverGain = str2double(headerACQP.ACQ_jobs{5});
    end
else
    receiverGain = 1;
end


% TE/TR
te = headerMethod.PVM_EchoTime;
tr = headerMethod.PVM_RepetitionTime;

% Sequence string
pulprog = headerACQP.PULPROG;
% Parse it
if contains(pulprog, 'press', 'IgnoreCase', true)
    sequence = 'PRESS';
elseif contains(pulprog, ["slaser", "semilaser"], 'IgnoreCase', true)
    sequence = 'sLASER';
elseif contains(pulprog, 'steam', 'IgnoreCase', true)
    sequence = 'STEAM';
elseif contains(pulprog, 'laser', 'IgnoreCase', true)
    sequence = 'LASER';
elseif contains(pulprog, 'special', 'IgnoreCase', true)
    sequence = 'SPECIAL';
else
    sequence = 'Unknown';
end

%Now load in the data.  Either raw or averaged data, depending on what was
%requested via the 'rawData' parameter. 
if strcmpi(rawData,'y')
    if contains(version, ["PV 6", "PV 7", "PV-7", "PV-360"])
        fileRaw     = fullfile(inDir, 'rawdata.job0');
        fids_raw    = readBrukerRaw(fileRaw, 'int32');
    elseif contains(version,'PV 5')
        fileRaw     = fullfile(inDir, 'fid.raw');
        fids_raw    = readBrukerRaw(fileRaw, 'int');
    end

    averages=rawAverages;  %since these data are uncombined;
    out.flags.averaged=0; %make the flags structure

    %If there are multiple receivers *I think* that these always get stored
    %separately by default in the fid.raw file.  Therefore, at this stage, 
    %if this is a PV360 dataset with multiple receivers, we need to reshape 
    %the dataset again:
    if ~contains(version,'PV 5') && multiRcvrs
        % Reshape into a Npts x Ncoils x Naverages x Nrepetitions array
        fids_raw=reshape(fids_raw,rawDataPoints,Nrcvrs,rawAverages,rawRepetitions);
        %Permute so that time is along 1st dimension, averages is along 2nd 
        %dimension, and coils is along 3rd dimension:
        fids_raw=permute(fids_raw,[1,3,2,4]); 

    elseif contains(version,'PV 5') && multiRcvrs
        % PV 5 does not appear to store individual coil channels in the
        % job0 file, but still preserves individual transients

        % Reshape into a Npts x Naverages x Nrepetitions array
        fids_raw=reshape(fids_raw,rawDataPoints,rawAverages,rawRepetitions);

    elseif ~contains(version,'PV 5') && ~multiRcvrs
        % Found a dataset with only one receiver and one repetition
        fids_raw=reshape(fids_raw,rawDataPoints,rawAverages,rawRepetitions);

    end

    % Normalize by the receiver gain for multi-channel data
    % This appears to be necessary but I can't quite figure out the
    % condition for that
    % Scale by receiver gain (if that has been previously extracted)
    if receiverGain
        fids_raw = fids_raw ./ receiverGain;
    end


elseif strcmpi(rawData,'n')
    %REQUEST PROCESSED DATA ONLY:  Use the FID file instead of fid.raw.
    if contains(version, ["PV-360"])
        fileRaw     = fullfile(inDir, 'pdata', '1', 'fid_proc.64');
        fids_raw    = readBrukerRaw(fileRaw, 'float64');
    else
        fileRaw     = fullfile(inDir, 'fid');
        fids_raw    = readBrukerRaw(fileRaw, 'int');
    end

    % Combined data have the averages combined, but the repetitions separate
    fids_raw = reshape(fids_raw, rawDataPoints, 1, rawRepetitions);

    % Remove singletons
    fids_raw = squeeze(fids_raw);

    averages=1; %since we requested the combined data.
    out.flags.averaged=1; %make the flags structure

else
    error('ERROR:  rawData variable not recognized.  Options are ''y'' or ''n''.');
end


%Perform left-shifting to remove points before the echo
fids_trunc=fids_raw(leftshift+1:end,:,:,:);

%replace the left-shifted points with zeros at the end
fids=padarray(fids_trunc, [leftshift,0],'post');

sz=size(fids); %size of the array

%calculate the dwelltime:
dwelltime=1/spectralwidth;

%calculate the time scale
t=[0:dwelltime:(rawDataPoints-1)*dwelltime];

%calculate the ppm scale
fmax=spectralwidth/2;
f=[fmax:-2*fmax/(rawDataPoints-1):-fmax];
ppm=f/(txfrq/1e6)+centerfreq;

% Apply ref freq shift (the difference between txfrq and txfrq_ref)
if contains(version, ["PV-360"])
    tmat=repmat(t',[1 sz(2:end)]);
    fids=fids.*exp(-1i*tmat*(txfrq_ref-txfrq)*2*pi);
end

%Do the fourier transform
specs=fftshift(ifft(fids,[],1),1);

%specify the dims:
%Time dimension:
dims.t=1;%the time dimension is always the 1st dimension
dims.extras=0; %so far, no data files with extra dimensions

%Coils dimension:
%As far as I am aware, the RF coils are only stored separately in PV360.
if ~contains(version,'PV 5') && multiRcvrs && strcmpi(rawData,'y')
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

if strcmpi(rawData,'y')
    dims.averages=2;
    if rawRepetitions>1
        dims.subSpecs=4;
    else
        dims.subSpecs=0;
    end
elseif strcmpi(rawData,'n')
    dims.averages=0;
    if rawRepetitions>1
        dims.subSpecs=2;
    else
        dims.subSpecs=0;
    end
end

%Specify the number of subspecs. If NRepetitions > 1, subSpecs will hold 
%the Repetitions dimensions (since it is the outer loop per Bruker manual)
subspecs=rawRepetitions;
rawSubspecs=rawRepetitions;


%% NOW TRY LOADING IN THE REFERENCE SCAN DATA (IF IT EXISTS)
% Reference data only exist for PV6 upwards
isRef=false;
if contains(version, ["PV 6", "PV 7", "PV-7", "PV-360.2"])
    isRef=exist(fullfile(inDir, 'fid.refscan'));
elseif contains(version,'PV-360.3')
    isRef=exist(fullfile(inDir, 'pdata', '1', 'fid_refscan.64'));
end

if isRef
    if contains(version, ["PV 6", "PV 7", "PV-7", "PV-360.2"])
        fileRef     = fullfile(inDir, 'fid.refscan');
        fids_ref    = readBrukerRaw(fileRef, 'int32');
    elseif contains(version,'PV-360.3')
        fileRef     = fullfile(inDir, 'pdata', '1', 'fid_refscan.64');
        fids_ref    = readBrukerRaw(fileRef, 'float64');
    end

    %Find the number of averages in the ref dataset:
    if contains(version,'PV 5')
        % These are navigator echoes acquired every TR
        rawAverages_ref=rawAverages;
        averages_ref=rawAverages_ref;
        ref.flags.averaged=0;
    else
        % These are water reference scans that are averaged
        rawAverages_ref = headerMethod.PVM_RefScanNA;
        averages_ref=1;
        ref.flags.averaged=1;
    end

    fids_ref=reshape(fids_ref,[],averages_ref);
    fids_ref_trunc=fids_ref(leftshift+1:end,:);
    reffids=padarray(fids_ref_trunc, [leftshift,0],'post');

    % Apply ref freq shift (the difference between txfrq and txfrq_ref)
    % (but not if it's PV360)
    sz_ref=size(fids_ref); %size of the array
    if ~contains(version, ["PV-360"])
        tmat=repmat(t',[1 sz_ref(2:end)]);
        reffids=reffids.*exp(-1i*tmat*(txfrq_ref-txfrq)*2*pi);
    end

    refspecs=fftshift(ifft(reffids,[],1),1);

    %specify the dims
    refdims.t=1;
    refdims.coils=0;
    if averages_ref>1
        refdims.averages=2;
    else
        refdims.averages=0;
    end    
    refdims.subSpecs=0;
    refdims.extras=0;

    %Specify the number of subspecs. For ref data, I think these are always
    %1
    subspecs_ref=1;
    rawSubspecs_ref=1;

else
    %Ref scans not found.  Print warning:
    disp('WARNING REFERENCE SCANS NOT FOUND.  RETURNING EMPTY REF STRUCTURE.');
end

%% Now load navigator scans
% For now, I know they're exported in PV 5 in fid.ref, and also in
% rawdata.job1 in the others, but I haven't figured out that format yet, so
% I'll skip PV6+ for now
isNav = false;
if contains(version,'PV 5')
    isNav=exist(fullfile(inDir, 'fid.ref'));
% elseif contains(version, ["PV 6", "PV 7", "PV-7", "PV-360"])
%     isNav=exist(fullfile(inDir, 'rawdata.job1'));
end

if isNav
    if contains(version,'PV 5')
        fileNav     = fullfile(inDir, 'fid.ref');
        fids_nav    = readBrukerRaw(fileNav, 'int');

%     elseif contains(version,'PV 6') || contains(version,'PV 7') || contains(version,'PV-360.2')
%         % Get TX freq of the reference, otherwise the reference ppm axis
%         % will not be correct
%         data = fopen(fullfile(inDir, 'fid.refscan'));
%         ref_data=fread(data,'int32');
%     elseif contains(version,'PV-360.3')
%         data = fopen(fullfile(inDir, 'pdata', '1', 'fid_refscan.64'));
%         ref_data=fread(data,'float64');

    end

    %Find the number of averages in the ref dataset:
    if contains(version,'PV 5')
        rawAverages_nav=rawAverages;
        averages_nav=rawAverages_nav;
%     else
%         rawAverages_ref = 1;
%         averages_nav=rawAverages_nav;
    end

    fids_nav=reshape(fids_nav,[],rawAverages_nav);
    fids_nav_trunc=fids_nav(leftshift+1:end,:);
    navfids=padarray(fids_nav_trunc, [leftshift,0],'post');

%     % Apply ref freq shift (the difference between txfrq and txfrq_ref)
%     % (but not if it's PV360)
%     if ~contains(version, ["PV-360"])
%         tmat=repmat(t',[1 sz(2:end)]);
%         navfids=navfids.*exp(-1i*tmat*(txfrq_ref-txfrq)*2*pi);
%     end

    navspecs=fftshift(ifft(navfids,[],1),1);
    sz_nav=size(navspecs);
    nav.flags.averaged=0;
    %specify the dims
    navdims.t=1;
    navdims.coils=0;
    if rawAverages_nav>1
        navdims.averages=2;
    else
        navdims.averages=0;
    end    
    navdims.subSpecs=0;
    navdims.extras=0;
else
    %Ref scans not found.  Print warning:
    disp('WARNING NAVIGATOR SCANS NOT FOUND.  RETURNING EMPTY NAV STRUCTURE.');
end

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
    ref.subspecs=subspecs_ref;
    ref.rawSubspecs=rawSubspecs_ref;
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


if isNav
    %FILLING IN DATA STRUCTURE FOR THE FID.REF DATA
    nav.fids=navfids;
    nav.specs=navspecs;
    nav.sz=sz_nav;
    nav.ppm=ppm;
    nav.t=t;
    nav.spectralwidth=spectralwidth;
    nav.dwelltime=dwelltime;
    nav.txfrq=txfrq;
    nav.date=date;
    nav.dims=navdims;
    nav.Bo=Bo;
    nav.averages=averages_nav;
    nav.rawAverages=rawAverages_nav;
    nav.subspecs=subspecs;
    nav.rawSubspecs=rawSubspecs;
    nav.seq=sequence;
    nav.te=te;
    nav.tr=tr;
    nav.pointsToLeftshift=0;
    nav.version=version;
    
    
    %FILLING IN THE FLAGS FOR THE NAV DATA
    nav.flags.writtentostruct=1;
    nav.flags.gotparams=1;
    nav.flags.leftshifted=1;
    nav.flags.filtered=0;
    nav.flags.zeropadded=0;
    nav.flags.freqcorrected=0;
    nav.flags.phasecorrected=0;
    
    nav.flags.addedrcvrs=1;
    nav.flags.subtracted=0;
    nav.flags.writtentotext=0;
    nav.flags.downsampled=0;
    nav.flags.avgNormalized=0;
    if nav.dims.subSpecs==0
        nav.flags.isFourSteps=0;
    else
        nav.flags.isFourSteps=(nav.sz(nav.dims.subSpecs)==4);
    end
else
    %NAV NOT FOUND.  RETURN EMPTY STRUCTURE.
    nav = [];
end

%% Store coil combination coefficients

% Only return coil combination coefficients if raw data are requested.
if multiRcvrs && strcmpi(rawData,'y')

    % All phases in PVM_ArrayPhase are calculated relative to the first 
    % channel, so the phases we need to feed into the coil combination are
    % the negative of that
    if length(rcvrPhases) > 1
        rcvrPhases = rcvrPhases(1) - rcvrPhases;
    else
        % I've seen single-coil data with a non-zero value here that
        % perfectly corresponds to the phase of that single channel. I'll
        % save that, but won't apply it here - can do that in the next step
        % outside of this function
    end
    coilcombos.ph   = rcvrPhases;

    % I *believe* that Bruker simply adds up the channel signals without
    % any weighting. This also means that the normalization of the various
    % FID-A functions to combine the signal will afterwards need to be
    % undone (by scaling by the number of channels)
    coilcombos.sig  = ones(size(rcvrPhases));

else
    % If only one coil, store trivial output
    coilcombos.ph   = 0;
    coilcombos.sig  = 1;
end

end



function fids_raw = readBrukerRaw(fileRaw, formatRaw)
% This subroutine reads the complex data from 'fileRaw' in the format 
% 'formatRaw' and returns a complex time-domain FID array

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

function coilcombos_emp = getCoilCombos(in, mode)
% This subroutine empirically calculates the coil combination coefficients
% to compare them against the ones that may be stored in the headers.
% Only here for validation/testing purposes.

% Get first point
point = 1;

coilcombos_emp.ph=zeros(in.sz(in.dims.coils),1);
coilcombos_emp.sig=zeros(in.sz(in.dims.coils),1);

for n=1:in.sz(in.dims.coils)
    coilcombos_emp.ph(n)=phase(in.fids(point,n,1,1))*180/pi; %in [degrees]
    switch mode
        case 'w'
            coilcombos_emp.sig(n)=abs(in.fids(point,n,1,1));
        case 'h'
            S=abs(in.fids(point,n,1,1));
            N=std(in.fids(end-100:end,n,1,1));
            coilcombos_emp.sig(n)=(S/(N.^2));
    end
end

%Now normalize the coilcombos.sig so that the max amplitude is 1;
coilcombos_emp.sig=coilcombos_emp.sig/max(coilcombos_emp.sig);

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
