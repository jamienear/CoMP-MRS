%compMRS_rfLoadFromMethod.m
%Jamie Near, Sunnybrook Research Institute, 2025
%
% USAGE:
% [out]=compMRS_rfLoadFromMethod(inputFile);
%
% DESCRIPTION:
% This script loads the RF waveform information from a Bruker Method file,
% including both the waveform (if saved) and the other relevant RF
% parameters.  The output of this function is a single structure with three
% fields:  1.  The RF pulse waveform, stored in FID-A rf pulse structure
% format;  2. The RF parameter info, stored as cell array; and 3. The RF
% pulse type (e.g. 'hermite', 'sinc', 'calcuated', etc.)
%
% INPUTS:
% inputFile:    The filename of the method file containing the RF pulse
%               info.
% 
% OUTPUTS:
% out:      N-element structure array (where N is the number of pulses, usually 3) 
%           with the following fields:
%               -RFobj:     RF pulse waveform structure for nth waveform, in FID-A structure format.
%               -RFparams: RF parameter info, stored as a cell array.
%               -PVver: Contains the paravision version information:
%               -RFtype: Type of RF pulse (e.g. 'hermite','sinc','calculated',etc.)


function [out] = compMRS_rfLoadFromMethod(inputFile)

%First load the method file using compMRS_loadMethod:
method = compMRS_loadMethod (inputFile);

%Now extract the RF pulse waveform and parameter info:
%Parameters first;
if isfield(method,'VoxPul1');
    params1=method.VoxPul1;
end
if isfield(method,'VoxPul2');
    params2=method.VoxPul2;
end
if isfield(method,'VoxPul3');
    params3=method.VoxPul3;
end

%Now waveforms
if isfield(method,'VoxPul1Shape')
    waveform1 = method.VoxPul1Shape;
end
if isfield(method,'VoxPul2Shape')
    waveform2 = method.VoxPul2Shape;
end
if isfield(method,'VoxPul3Shape')
    waveform3 = method.VoxPul3Shape;
end

%Waveforms are a vector with amplitude and phase interleaved.  Separate
%amplitude and phase as follows:
amplitude1=waveform1(1:2:end);
phase1=waveform1(2:2:end);

amplitude2=waveform2(1:2:end);
phase2=waveform2(2:2:end);

amplitude3=waveform3(1:2:end);
phase3=waveform3(2:2:end);

%Now redefine the RF waveforms as Nx3 arrays:
waveform1=zeros(length(amplitude1),3);
waveform1(:,1)=phase1;
waveform1(:,2)=amplitude1;
waveform1(:,3)=ones(length(amplitude1),1);

waveform2=zeros(length(amplitude2),3);
waveform2(:,1)=phase2;
waveform2(:,2)=amplitude2;
waveform2(:,3)=ones(length(amplitude2),1);

waveform3=zeros(length(amplitude3),3);
waveform3(:,1)=phase3;
waveform3(:,2)=amplitude3;
waveform3(:,3)=ones(length(amplitude3),1);


%Now convert all of those waveforms into FID-A objects using
%io_loadRFwaveform.m:

%DISPLAY THE PULSE SEQUENCE NAME:
disp(['The pulse sequence is' method.Method '.']);

%The first pulse will always be an excitation pulse:
disp('LOADING THE FIRST RF PULSE:')
RFobj1=io_loadRFwaveform(waveform1,'exc');

%The next two pulses will be excitation pulses if the sequences is STEAM,
%but will be refocusing pulses otherwise:

if contains(method.Method,'steam','IgnoreCase',true);
    disp('LOADING THE 2nd RF PULSE:')
    RFobj2=io_loadRFwaveform(waveform2,'exc');
    disp('LOADING THE 3rd RF PULSE:')
    RFobj3=io_loadRFwaveform(waveform3,'exc');
else
    disp('LOADING THE 2nd RF PULSE:')
    RFobj2=io_loadRFwaveform(waveform2,'ref');
    disp('LOADING THE 3rd RF PULSE:')
    RFobj3=io_loadRFwaveform(waveform3,'ref');
end

%Find the RF pulse types:
if isfield(method,'VoxPul1Enum');
    type1=method.VoxPul1Enum;
end
if isfield(method,'VoxPul2Enum');
    type2=method.VoxPul2Enum;
end
if isfield(method,'VoxPul3Enum');
    type3=method.VoxPul3Enum;
end

%Find the PV version:
if isfield(method,'TITLE');
    title=method.TITLE;
end


%Finally, save the output:
out(1).RFobj=RFobj1;
out(1).RFparams=params1;
out(1).RFtype=type1;
out(1).PVver=title;

out(2).RFobj=RFobj2;
out(2).RFparams=params2;
out(2).RFtype=type2;
out(2).PVver=title;

out(3).RFobj=RFobj3;
out(3).RFparams=params3;
out(3).RFtype=type3;
out(3).PVver=title;



