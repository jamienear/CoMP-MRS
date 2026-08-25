% compMRS_getSimPars.m
% Jamie Near and Diana Rotaru, 2026
%
% USAGE:
% simPars = compMRS_getSimPars(DPid);
%
% DESCRIPTION:
% This function retrieves the set of parameters needed for basis set
% simulation for the CoMP-MRS dataset.  The parameters are found from the 
% header information of the first subject in the DP.  Both Bruker and
% Varian/Agilent DPs are supported. 
%
% INPUTS:
% DPid          = the Data Packet ID number (e.g. 'DP01', 'DP02',etc.)
%
% OUTPUTS:
% simPars       = a structure variable containing the relevant parameters
%                 for basis set simulation. 
%

function simPars = compMRS_getSimPars(DPid);

% First check the vendor using DPcheck:
check = compMRS_DPcheck(DPid);

% Make a list of all subjects and sessions, and the svs data path.  Only the 
% first subject/session will be used (all subjects in a DP must be collected
% using the same acquisition).  
subjs=dir([DPid filesep 'sub*']);
sess=dir([DPid filesep subjs(1).name filesep 'ses*']);
svspath = dir([DPid filesep subjs(1).name filesep sess(1).name filesep 'mrs' filesep '*svs']);


% Process for getting the header info depends on vendor:
if strcmp(check.vendor(1),'BRUKER')
    %read the method file:
    method = compMRS_loadMethod([svspath.folder filesep svspath.name filesep 'method']);
    acqp = compMRS_parseBrukerFormat([svspath.folder filesep svspath.name filesep 'acqp']);

    %Now extract the relevant params from the method file:
    sequence=method.Method;
    
    %Field strength:
    simPars.Bfield = acqp.SFO1/42.577; %[in Tesla]

    %initialize the "isCalcWaveform" boolean variable to false:
    isCalcWaveform=false;

    %Initialize all pulse sequence booleans to false:
    isPRESS = false;
    isSTEAM = false;
    isSPECIAL = false;
    isLASER = false;

    %sequence timings and RF pulse waveforms:
    if contains(sequence,'press','IgnoreCase',true)
        isPRESS = true;
        simPars.seq = 'PRESS';
        simPars.tau1 = method.TE1;
        simPars.tau2 = method.TE2;

        if isfield(method,'VoxPul2Enum')
            simPars.rfName = [char(method.VoxPul2Enum) '.rfc'];
            VP2=char(method.VoxPul2);
        elseif isfield(method,'VoxPulse2Enum')
            simPars.rfName = [char(method.VoxPulse2Enum) '.rfc'];
            VP2=char(method.VoxPulse2);
        end

        VP2=VP2(2:end);
        VP2_values = split(VP2,',');
        simPars.refTp = str2double(VP2_values{1});
        simPars.flipAngle=str2double(VP2_values{3});

        if contains(simPars.rfName,'calculated','IgnoreCase',true)
            isCalcWaveform = true;
            sharpness = str2double(VP2_values{5});
            bwfac = str2double(VP2_values{6});
        end

    elseif contains(sequence,'steam','IgnoreCase',true)
        isSTEAM = true;
        simPars.seq = 'STEAM';
        simPars.tau1 = method.PVM_EchoTime;
        simPars.tau2 = method.StTM;

        if isfield(method,'VoxPul1Enum')
            simPars.rfName = [char(method.VoxPul1Enum) '.exc'];
            VP1=char(method.VoxPul1);
        elseif isfield(method,'VoxPulse1Enum')
            simPars.rfName = [char(method.VoxPulse1Enum) '.exc'];
            VP1=char(method.VoxPulse1);
        end

        VP1=VP1(2:end);
        VP1_values = split(VP1,',');
        simPars.Tp = str2double(VP1_values{1});
        simPars.flipAngle=str2double(VP1_values{3});

        if contains(simPars.rfName,'calculated','IgnoreCase',true)
            isCalcWaveform = true;
            sharpness = str2double(VP1_values{5});
            bwfac = str2double(VP1_values{6});
        end

    elseif contains(sequence,'laser','IgnoreCase',true)
        isLASER = true;
        simPars.seq = 'LASER';
        simPars.te = method.PVM_EchoTime;

        if isfield(method,'VoxPul2Enum')
            simPars.rfName = [char(method.VoxPul2Enum) '.rfc'];
            VP2=char(method.VoxPul2);
        elseif isfield(method,'VoxPulse2Enum')
            simPars.rfName = [char(method.VoxPulse2Enum) '.rfc'];
            VP2=char(method.VoxPulse2);
        end

        VP2=VP2(2:end);
        VP2_values = split(VP2,',');
        simPars.refTp = str2double(VP2_values{1});
        simPars.flipAngle=str2double(VP2_values{3});

        if contains(simPars.rfName,'calculated','IgnoreCase',true)
            isCalcWaveform = true;
            sharpness = str2double(VP2_values{5});
            bwfac = str2double(VP2_values{6});
        end

    elseif contains(sequence,'special','IgnoreCase',true)
        isSPECIAL = true;
        simPars.seq = 'SPECIAL';
        simPars.tau1 = method.PVM_EchoTime;

        if isfield(method,'VoxPul3Enum')
            simPars.rfName = [char(method.VoxPul3Enum) '.rfc'];
            VP3=char(method.VoxPul3);
        elseif isfield(method,'VoxPulse3Enum')
            simPars.rfName = [char(method.VoxPulse3Enum) '.rfc'];
            VP3=char(method.VoxPulse3);
        end

        VP3=VP3(2:end);
        VP3_values = split(VP3,',');
        simPars.refTp = str2double(VP3_values{1});
        simPars.flipAngle=str2double(VP3_values{3});

        if contains(simPars.rfName,'calculated','IgnoreCase',true)
            isCalcWaveform = true;
            sharpness = str2double(VP3_values{5});
            bwfac = str2double(VP3_values{6});
        end

    end

    %Check if the waveform is calculated.  If so, re-name the waveform
    %accordingly:
    if isCalcWaveform
        if isSTEAM
            if bwfac == 4200 && sharpness == 3
                simPars.rfName = 'brukerCalc_exc_sh3_bwf4200.txt';
            else
                error('ERROR: STEAM - No matching BWFAC and Sharpness values found.  ABORTING!');
            end
        elseif isLASER
            if floor(bwfac) == 27431 && sharpness == 3
                simPars.rfName = 'brukerCalc_HSinv_sh3_bwf27431.txt';
            else
                error('ERROR: sLASER - No matching BWFAC and Sharpness values found.  ABORTING!');
            end
        elseif isPRESS || isSPECIAL
            if bwfac == 3400 && sharpness ==3
                simPars.rfName = 'brukerCalc_ref_sh3_bw3400.txt';
            else
                error('ERROR: PRESS/SPECIAL - No matching BWFAC and Sharpness values found.  ABORTING!');
            end
        end
    end

elseif strcmp(check.vendor(1),'VARIAN')
    %read the procpar file:
    par = readprocpar2([svspath.folder filesep svspath.name filesep 'procpar']);

    %Now extract the relevant params from the method file:
    sequence=par.seqfil.value;
    
    %Field strength:
    simPars.Bfield = par.B0.value / 10000; %Converted from [G] to [Tesla]

    %Initialize all pulse sequence booleans to false:
    isPRESS = false;
    isSTEAM = false;
    isSPECIAL = false;
    isLASER = false;
    
    %sequence timings and RF pulse waveforms:
    if contains(sequence,'press','IgnoreCase',true)
        isPRESS = true;
        simPars.seq = 'PRESS';
        error('ERROR:  No Varian/Agilent DPs currently using PRESS sequence')

    elseif contains(sequence,'steam','IgnoreCase',true)
        isSTEAM = true;
        simPars.seq = 'STEAM';
        simPars.tau1 = par.te.value * 1000; %convert from [s] to [ms]
        simPars.tau2 = par.tm.value * 1000; %convert from [s] to [ms]        
        simPars.rfName = [par.p1pat.value{1} '.RF'];
        simPars.Tp = par.pw90.value / 1000; %convert from [us] to [ms]
        simPars.flipAngle= 90; %[degrees] (Hard coding for now until I can find flip angle in procpar).

    elseif contains(sequence,'laser','IgnoreCase',true)
        isLASER = true;
        simPars.seq = 'STEAM';
        simPars.seq = 'STEAM';
        simPars.te = par.te.value * 1000; %convert from [s] to [ms]
        simPars.rfName = [par.pat180Y.value{1} '.RF'];
        simPars.refTp = par.pw180.value / 1000; %convert from [us] to [ms]
        simPars.flipAngle=180; %[degrees] (Hard coding for now until I can find flip angle in procpar).

    elseif contains(sequence,'special','IgnoreCase',true) || contains(sequence,'isise','IgnoreCase',true)
        isSPECIAL = true;
        simPars.tau1 = par.te.value * 1000;
        simPars.rfName = [par.p2pat.value{1} '.RF'];
        simPars.refTp = par.pw180.value / 1000; % converting from [s] to [ms]
        simPars.flipAngle=180; %[degrees].  (Hard coding for now until I can find the flip angle in procpar).
    
    end

end

% Now, load the RF pulse waveform and replace the simPars RF waveform with 
% the resulting FID-A structure:
if isPRESS || isSPECIAL || isLASER
    RF=io_loadRFwaveform(simPars.rfName,'ref');
    simPars.refocWaveform = RF;
elseif isSTEAM
    RF=io_loadRFwaveform(simPars.rfName,'exc');
    simPars.excWaveform = RF;
end

%Fill out the remaining fields in the simPars structure variable.  
%Some fields can be hard-coded, some depend on field strength:
if simPars.Bfield<8
    simPars.lw = 1;         %Simulation linewidth in [Hz]
    simPars.Npts = 8192;    %Number of spectral points
    simPars.sw = 4000;      %Spectral width in [Hz]
elseif simPars.Bfield>=8 && simPars.Bfield<12
    simPars.lw = 1.5;       %Simulation linewidth in [Hz]
    simPars.Npts = 12288;   %Number of spectral points
    simPars.sw = 6000;      %Spectral width in [Hz]
elseif simPars.Bfield>=12
    simPars.lw = 2;         %Simulation linewidth in [Hz]
    simPars.Npts = 16384;   %Number of spectral points
    simPars.sw = 8000;      %Spectral width in [Hz]
end

simPars.thkX = 2;           %Voxel size in x-direction in [arb units]
simPars.thkY = 2;           %Voxel size in y-direction in [arb units]
simPars.fovX = 4;           %Simulation FOV in x-direction in [arb units]
simPars.fovY = 4;           %Simulation FOV in y-direction in [arb units]
simPars.nX = 32;            %Number of spatial points to simulate in x-direction
simPars.nY = 32;            %Number of spatial points to simulate in y-direction
simPars.centreFreq = 2.7;   %Centre frequency in [ppm]
simPars.lineshape = 'L';    %Lorentzian lineshape for metabolite simulations

%**************Done generating the simPars structure

