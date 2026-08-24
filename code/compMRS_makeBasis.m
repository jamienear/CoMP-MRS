% compMRS_makeBasis.m
% Jamie Near and Diana Rotaru, 2025/2026
%
% USAGE:
% basis = compMRS_makeBasis(DPid)
%
% DESCRIPTION:
% This function is used to automatically generate an LCModel basis set for 
% any data packet (DP) from the CoMP-MRS dataset.  The function will first 
% read the header information from the desired DP to record relevant params
% from the DP's MRS acquisition (i.e. field strength, pulse sequence used, 
% pulse sequence timings (i.e. TE, TM, etc.), RF pulse waveforms, RF pulse
% durations, etc.  Then it will call one of the built in simulation
% functions (compMRS_simPRESS, compMRS_simSTEAM, compMRS_simsLASER, or
% compMRS_simSPECIAL, as appropriate), using the input parameters
% determined from the data header.  It will generate basis functions for 
% each metabolite of interest, in LCModel .RAW format.  The function then 
% makes a 'makebasis.in' file to wrap the .RAW files, and then makes a
% system call to LCModel's 'makebasis' method to generate the basis set. 
%
% INPUTS:
% DPid       = the Data Packet ID number (e.g. 'DP01', 'DP02', etc.);
%
% OUTPITS:
% basis      = a struct variable whose fields are FID-A data struct
%             variables specifying the basis spectra for each individual
%             metabolite.

function basis = compMRS_makeBasis(DPid)

% First check the vendor using DPcheck:
check = compMRS_DPcheck(DPid);

% Make a list of all subjects and sessions, and the svs data path.  Only the 
% first subject/session will be used (all subjects in a DP must be collected
% using the same acquisition).  
subjs=dir([DPid filesep 'sub*']);
sess=dir([DPid filesep subjs(1).name filesep 'ses*']);
svspath = dir([DPid filesep subjs(1).name filesep sess(1).name filesep 'mrs' filesep '*svs']);

%Make a list of metabolites to include in basis set:
%metabs = {'Ala','Asp','Cr','GABA','Glc','Gln','Glu','GPC','GSH','Ins','Lac','NAA','NAAG','PCh','PCr','PE','Ser','Tau','Ref0ppm'};
metabs = {'Lac','Ref0ppm'};

%Make an inline function definition to make the MM spin systems:
makeMM = @(n,shift,scale) struct('J',0,'shifts',shift,'name',['MM' num2str(n)],'scaleFactor',scale);
MMshifts = [0.89, 1.20, 1.39, 1.66, 2.02, 2.26, 2.97, 3.18, 3.84];  %From Fowler et al. 2021
MMscales = [   3,    2,    2,    2,    2,    2,    2,    2,    2];  %
MMlws =    [  34,   27,   31,   61,   82,   23,   27,   29,   88];  %From Fowler et al. 2021

%Now make the macromolecule basis functions:
for n=1:9
    sysMM{n}.sys=makeMM(n,MMshifts(n),MMscales(n));
    sysMM{n}.lw = MMlws(n);
    sysMM{n}.name = ['MM' num2str(n)];
end

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
        simPars.tau1 = method.TE1;
        simPars.tau2 = method.TE2;

        if isfield(method,'VoxPul2Enum')
            simPars.refocWaveform = [char(method.VoxPul2Enum) '.rfc'];
            VP2=char(method.VoxPul2);
        elseif isfield(method,'VoxPulse2Enum')
            simPars.refocWaveform = [char(method.VoxPulse2Enum) '.rfc'];
            VP2=char(method.VoxPulse2);
        end

        VP2=VP2(2:end);
        VP2_values = split(VP2,',');
        simPars.refTp = str2double(VP2_values{1});
        simPars.flipAngle=str2double(VP2_values{3});

        if contains(simPars.refocWaveform,'calculated','IgnoreCase',true)
            isCalcWaveform = true;
            sharpness = str2double(VP2_values{5});
            bwfac = str2double(VP2_values{6});
        end

    elseif contains(sequence,'steam','IgnoreCase',true)
        isSTEAM = true;
        simPars.tau1 = method.PVM_EchoTime;
        simPars.tau2 = method.StTM;

        if isfield(method,'VoxPul1Enum')
            simPars.excWaveform = [char(method.VoxPul1Enum) '.exc'];
            VP1=char(method.VoxPul1);
        elseif isfield(method,'VoxPulse1Enum')
            simPars.excWaveform = [char(method.VoxPulse1Enum) '.exc'];
            VP1=char(method.VoxPulse1);
        end

        VP1=VP1(2:end);
        VP1_values = split(VP1,',');
        simPars.Tp = str2double(VP1_values{1});
        simPars.flipAngle=str2double(VP1_values{3});

        if contains(simPars.excWaveform,'calculated','IgnoreCase',true)
            isCalcWaveform = true;
            sharpness = str2double(VP1_values{5});
            bwfac = str2double(VP1_values{6});
        end

    elseif contains(sequence,'laser','IgnoreCase',true)
        isLASER = true;
        simPars.te = method.PVM_EchoTime;

        if isfield(method,'VoxPul2Enum')
            simPars.refocWaveform = [char(method.VoxPul2Enum) '.rfc'];
            VP2=char(method.VoxPul2);
        elseif isfield(method,'VoxPulse2Enum')
            simPars.refocWaveform = [char(method.VoxPulse2Enum) '.rfc'];
            VP2=char(method.VoxPulse2);
        end

        VP2=VP2(2:end);
        VP2_values = split(VP2,',');
        simPars.refTp = str2double(VP2_values{1});
        simPars.flipAngle=str2double(VP2_values{3});

        if contains(simPars.refocWaveform,'calculated','IgnoreCase',true)
            isCalcWaveform = true;
            sharpness = str2double(VP2_values{5});
            bwfac = str2double(VP2_values{6});
        end

    elseif contains(sequence,'special','IgnoreCase',true)
        isSPECIAL = true;
        simPars.tau1 = method.PVM_EchoTime;

        if isfield(method,'VoxPul3Enum')
            simPars.refocWaveform = [char(method.VoxPul3Enum) '.rfc'];
            VP3=char(method.VoxPul3);
        elseif isfield(method,'VoxPulse3Enum')
            simPars.refocWaveform = [char(method.VoxPulse3Enum) '.rfc'];
            VP3=char(method.VoxPulse3);
        end

        VP3=VP3(2:end);
        VP3_values = split(VP3,',');
        simPars.refTp = str2double(VP3_values{1});
        simPars.flipAngle=str2double(VP3_values{3});

        if contains(simPars.excWaveform,'calculated','IgnoreCase',true)
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
                simPars.excWaveform = 'brukerCalc_exc_sh3_bwf4200.txt';
            else
                error('ERROR: STEAM - No matching BWFAC and Sharpness values found.  ABORTING!');
            end
        elseif isLASER
            if floor(bwfac) == 27431 && sharpness == 3
                simPars.refocWaveform = 'brukerCalc_HSinv_sh3_bwf27431.txt';
            else
                error('ERROR: sLASER - No matching BWFAC and Sharpness values found.  ABORTING!');
            end
        elseif isPRESS || isSPECIAL
            if bwfac == 3400 && sharpness ==3
                simPars.refocWaveform = 'brukerCalc_ref_sh3_bw3400.txt'
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
        error('ERROR:  No Varian/Agilent DPs currently using PRESS sequence')

    elseif contains(sequence,'steam','IgnoreCase',true)
        isSTEAM = true;
        simPars.tau1 = par.te.value * 1000; %convert from [s] to [ms]
        simPars.tau2 = par.tm.value * 1000; %convert from [s] to [ms]        
        simPars.excWaveform = [par.p1pat.value{1} '.RF']
        simPars.Tp = par.pw90.value / 1000; %convert from [us] to [ms]
        simPars.flipAngle= 90; %[degrees] (Hard coding for now until I can find flip angle in procpar).

    elseif contains(sequence,'laser','IgnoreCase',true)
        isLASER = true;
        simPars.te = par.te.value * 1000; %convert from [s] to [ms]
        simPars.refocWaveform = [par.pat180Y.value{1} '.RF'];
        simPars.refTp = par.pw180.value / 1000; %convert from [us] to [ms]
        simPars.flipAngle=180; %[degrees] (Hard coding for now until I can find flip angle in procpar).

    elseif contains(sequence,'special','IgnoreCase',true) || contains(sequence,'isise','IgnoreCase',true)
        isSPECIAL = true;
        simPars.tau1 = par.te.value * 1000;
        simPars.refocWaveform = [par.p2pat.value{1} '.RF'];
        simPars.refTp = par.pw180.value / 1000; % converting from [s] to [ms]
        simPars.flipAngle=180; %[degrees].  (Hard coding for now until I can find the flip angle in procpar).
    
    end

end


% Now, load the RF pulse waveform and replace the simPars RF waveform with 
% the resulting FID-A structure:
if isPRESS || isSPECIAL || isLASER
    RF=io_loadRFwaveform(simPars.refocWaveform,'ref');
    simPars.refocWaveform = RF;
elseif isSTEAM
    RF=io_loadRFwaveform(simPars.excWaveform,'exc');
    simPars.excWaveform = RF;
end

%Fill out the remaining fields in the simPars structure variable.  
%Some fields can be hard-coded, some depend on field strength:
if simPars.Bfield<8
    simPars.lw = 1;         %Simulation linewidth in [Hz]
    %simPars.Npts = 8192;    %Number of spectral points
    simPars.Npts = 4096;
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


%Now run the simulations:
disp(['Simulating ' check.vendor{1} ' ' num2str(simPars.Bfield) ' Tesla ' char(sequence) ' Basis set for ' char(DPid) ' with following params:'])
simPars

%For the metabolites:
for n=1:length(metabs)
    disp(['Simulating ' metabs{n} '!']);
    eval(['load ' metabs{n}]);
    if isPRESS
        eval([metabs{n} '= compMRS_simPRESS(sys' metabs{n} ',simPars);']);
    elseif isSTEAM
        eval([metabs{n} '= compMRS_simSTEAM(sys' metabs{n} ',simPars);']);
    elseif isLASER
        eval([metabs{n} '= compMRS_simsLASER(sys' metabs{n} ',simPars);']);
    elseif isSPECIAL
        eval([metabs{n} '= compMRS_simSPECIAL(sys' metabs{n} ',simPars);']);
    end
end

%For the macromolecules:
simPars.lineshape = 'G'; %Gaussian lineshape for the MM simulations;
for n=1:length(sysMM)
    disp(['Simulating ' sysMM{n}.name '!']);
    simPars.lw=sysMM{n}.lw * simPars.Bfield / 7; %MM linewidths will scale with field strength
    if isPRESS
        eval([sysMM{n}.name '= compMRS_simPRESS(sysMM{n}.sys,simPars);']);
    elseif isSTEAM
        eval([sysMM{n}.name '= compMRS_simSTEAM(sysMM{n}.sys,simPars);']);
    elseif isLASER
        eval([sysMM{n}.name '= compMRS_simsLASER(sysMM{n}.sys,simPars);']);
    elseif isSPECIAL
        eval([sysMM{n}.name '= compMRS_simSPECIAL(sysMM{n}.sys,simPars);']);
    end
end
    
%Make an output directory in the DP folder
mkdir([DPid '/basis-set']);

%Now add the reference peak to all of the other metabolite and MM basis
%functions.  Then, shift them to be centered at 4.65 ppm.  Then write to
%.RAW format.  Also make an output cell array with all of the simulated
%spectra in it:
for n=1:length(metabs)-1
    eval([metabs{n} '=op_addScans(' metabs{n} ',' metabs{end} ');']);
    eval([metabs{n} '=op_movef0(' metabs{n} ',(4.65 -' num2str(simPars.centreFreq) ')*42.577*' num2str(simPars.Bfield) ');']);
    eval(['[~]=io_writelcmraw(' metabs{n} ',''' DPid '/basis-set/' metabs{n} '.RAW'',''' metabs{n} ''');']);
    eval(['basis{n}=' metabs{n} ';']);
end
for n=1:length(sysMM)
    eval([sysMM{n}.name '=op_addScans(' sysMM{n}.name ',' metabs{end} ');']);
    eval([sysMM{n}.name '=op_movef0(' sysMM{n}.name ',(4.65 -' num2str(simPars.centreFreq) ')*42.577*' num2str(simPars.Bfield) ');']);
    eval(['[~]=io_writelcmraw(' sysMM{n}.name ',''' DPid '/basis-set/' sysMM{n}.name '.RAW'',''' sysMM{n}.name ''');']);
    eval(['basis{length(metabs)-1 + n}=' sysMM{n}.name ';']);
end



    
