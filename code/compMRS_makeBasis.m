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
% basis      = a cell array whose fields are FID-A data struct
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
metabs = {'NAA','Lac','Ref0ppm'};  %Shorter list of metabolites for testing.  Uncomment line above for full list.

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


%Get the simpars structure:
simPars=compMRS_getSimPars(DPid);


%Now run the simulations:
disp(['Simulating ' check.vendor{1} ' ' num2str(simPars.Bfield) ' Tesla ' simPars.seq ' Basis set for ' char(DPid) ' with following params:'])
simPars

%For the metabolites:
for n=1:length(metabs)
    disp(['Simulating ' metabs{n} '!']);
    eval(['load ' metabs{n}]);
    if strcmp(simPars.seq,'PRESS')
        eval([metabs{n} '= compMRS_simPRESS(sys' metabs{n} ',simPars);']);
    elseif strcmp(simPars.seq,'STEAM')
        eval([metabs{n} '= compMRS_simSTEAM(sys' metabs{n} ',simPars);']);
    elseif strcmp(simPars.seq,'LASER')
        eval([metabs{n} '= compMRS_simsLASER(sys' metabs{n} ',simPars);']);
    elseif strcmp(simPars.seq,'SPECIAL')
        eval([metabs{n} '= compMRS_simSPECIAL(sys' metabs{n} ',simPars);']);
    end
end

%For the macromolecules:
simPars.lineshape = 'G'; %Gaussian lineshape for the MM simulations;
for n=1:length(sysMM)
    disp(['Simulating ' sysMM{n}.name '!']);
    simPars.lw=sysMM{n}.lw * simPars.Bfield / 7; %MM linewidths will scale with field strength
    if strcmp(simPars.seq,'PRESS')
        eval([sysMM{n}.name '= compMRS_simPRESS(sysMM{n}.sys,simPars);']);
    elseif strcmp(simPars.seq,'STEAM')
        eval([sysMM{n}.name '= compMRS_simSTEAM(sysMM{n}.sys,simPars);']);
    elseif strcmp(simPars.seq,'LASER')
        eval([sysMM{n}.name '= compMRS_simsLASER(sysMM{n}.sys,simPars);']);
    elseif strcmp(simPars.seq,'SPECIAL')
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



    
