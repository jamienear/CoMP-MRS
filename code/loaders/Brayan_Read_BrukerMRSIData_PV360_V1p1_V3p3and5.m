%% Open and display Bruker CSI 
%JM - 09/02/2021 

clear;close all;clc;


%% open FID
expnb=18;% Bruker data are stored in folders with numbers based on order they were acquired. for teh data I sent you I modified the names of these folders 
expwater=19;
reconb = 2;

% put the name of the folder where the Bruker data is stored
namepath=" "; %Put name in between quotation marks

fullpath=namepath+num2str(expnb)+"\fid"; %Open files from PV360 V.1.1
path_water=namepath+num2str(expwater)+"\fid";      

fileid=fopen(fullpath,'r','ieee-le'); %read binary format
waterid=fopen(path_water,'r','ieee-le');
if fileid == -1
    fullpath=namepath+num2str(expnb)+"\pdata\"+reconb+"\fid_proc.64"; %Open files from PV360 V.3.3
    path_water=namepath+num2str(expwater)+"\pdata\"+reconb+"\fid_proc.64";
    fileid=fopen(fullpath,'r','ieee-le'); %read binary format
    waterid=fopen(path_water,'r','ieee-le');
    if fileid == -1
        disp('Cannot open file');
        return
    end
end

%% infos 
Nav=1;
MatSize= ReadMatrixSizeBruker(namepath+num2str(expnb));
if(length(MatSize)==3)
    Nslices=MatSize(3);
else
    Nslices=1;
end
FidPoints= ReadSpectPointsBruker(namepath+num2str(expnb));%1024;% %complex points
FidPoints_water= ReadSpectPointsBruker(namepath+num2str(expwater));%1024;%
NMRFreq =ReadNMRFreqBruker(namepath+num2str(expnb));
minFreq =-0.5*NMRFreq;
maxFreq =0.5*NMRFreq;
NbPtForWaterPhAmp = 5;
samplerate =ReadBandwidthBruker(namepath+num2str(expnb)); %Hz
AcqFreqShift = ReadAcqFreqShiftBruker(namepath+num2str(expnb));
receiveroffset_ppm=4.7;
grpdly=77;  % these first points of FID are not ok -need to be taken out

%ppm scale
t=[0:1/samplerate:(FidPoints-1)/samplerate]; 
fmax=samplerate/2;
f=[fmax:-2*fmax/(FidPoints-1):-fmax];
scale_ppm=f/NMRFreq+receiveroffset_ppm;

AcqDelay = ReadTEBruker(namepath+num2str(expnb)); % in ms
voxel_size = [ReadVoxelSizeBruker(namepath+num2str(expnb)),ReadSliceThicknessBruker(namepath+num2str(expnb))]; % in mm
voxel_volume = voxel_size(1)*voxel_size(2)*voxel_size(3);
AcqTime = ReadAcquisitionTimeBruker(namepath+num2str(expnb))*1E-3; %in seconds
RepTime = ReadRepetitionTimeBruker(namepath+num2str(expnb))*1E-3; %in seconds
flipangle = ReadFlipAngleBruker(namepath+num2str(expnb)); %in degrees
Fs = samplerate;
Time_t=[0 :(FidPoints-1)]'/Fs;
Freqshift_1t=exp(2*pi*1i*Time_t*AcqFreqShift); %-2pi to go from k -> r

%Receiver Gains values
RG_met = ReadRGValueBruker(namepath+num2str(expnb));
RG_wat = ReadRGValueBruker(namepath+num2str(expwater));

%% Read data		
buffer=fread(fileid,'double'); %note: CSI Bruker format is double 
buffer_c=buffer(1:2:end)+sqrt(-1)*buffer(2:2:end);

wasser=fread(waterid,'double'); %note: CSI Bruker format is double 
wasser_c=wasser(1:2:end)+sqrt(-1)*wasser(2:2:end);

fid_mat_c=reshape(buffer_c, [FidPoints,MatSize(1)*MatSize(2),Nslices]);
was_mat_c=reshape(wasser_c, [FidPoints_water,MatSize(1)*MatSize(2),Nslices]);


fid_mat_c_shift = [fid_mat_c(grpdly:end,:,:);0*fid_mat_c(1:grpdly-1,:,:)];
was_mat_c_shift = [was_mat_c(grpdly:end,:,:);0*was_mat_c(1:grpdly-1,:,:)];


fid_mat_tkkn = conj(fft(fft(reshape(fid_mat_c_shift,[FidPoints,MatSize(1),MatSize(2),Nslices]),[],2),[],3)); % Restructuring the data + going to the k-space
was_mat_tkkn = conj(fft(fft(reshape(was_mat_c_shift,[FidPoints,MatSize(1),MatSize(2),Nslices]),[],2),[],3));

fid_mat_tkkn = fid_mat_tkkn.*Freqshift_1t; % We apply the frequency shift to center the water peak on the receiveroffset (ppm scale)
was_mat_tkkn = was_mat_tkkn.*Freqshift_1t;

%Normalization to have the RG between water and metabolites

if(RG_met~=RG_wat)
    RG_corr = (RG_met/RG_wat);
    was_mat_tkkn = was_mat_tkkn*(RG_met/RG_wat);
end