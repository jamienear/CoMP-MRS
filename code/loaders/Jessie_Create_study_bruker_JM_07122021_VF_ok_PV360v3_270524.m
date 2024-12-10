
%% From a Bruker exp. folder, generates MATLAB structures with necessary 
%params + data (fid or ser file) + fid.refscan (fidrefscan) if it exists. 
%This matrix can then be read by Nico's GUI.
% JM - 17/02/2021
%includes: - normalization by RG/nrep, reading the files not opened in
%topspin, grpdly, offset
%JM - 27/5/2024
%adapted to read 360v3 MRS: adapted the location of ser file + format int32-->double 
% included inside the reading of rawdatajob0 - from
% Read_bruker_rawdatajob0.... 
% updated the display to show exactly what is saved - better to make sure
% the normalisation is correct
clear; close all;


%% Enter exp nb. + directory and run the program 

folderexp="G:\47-Preclinical_initiative_CC_JN\SHAM759\";
explisttot=[10,15]
% expnb=6;
timeacq='';%'_BDL727';%'' or am or pm
dateexp='2024'; 
for kkkk=1:length(explisttot)
    expnb=explisttot(kkkk);
    clear study;
%% Parameters - method file
study.path=char(folderexp+num2str(expnb));
methodfile = fileread(folderexp +num2str(expnb)+ "\method");

%date and time
startind_date=strfind(methodfile,['$$ ' dateexp]); %adapt the data here is needed
date=methodfile(startind_date+3:startind_date+21);
study.time=date;

%literal names = acquisition time (sec), spectral width (Hz), sequence, np
%(complex nb *2), nucleus, resonance frequency B0, ppmoffset, tr(sec), RG,
% voxel dimensions (mm), acq type (MRS), nrep, nav
BRUKERparamnames=["PVM_SpecAcquisitionTime=","##$PVM_SpecSWH=( 1 )", ...
    "##$Method=","##$PVM_SpecMatrix=( 1 )","##$PVM_Nucleus1Enum=", ...
    ["##$PVM_FrqRef=( 8 )" + sprintf('\n')], ...
    ["##$PVM_FrqWorkOffsetPpm=( 8 )" + sprintf('\n')], ...
    "##$PVM_RepetitionTime=","##$PVM_RgValue=", ...
    ["##$PVM_VoxArrSize=( 1, 3 )" + sprintf('\n')],"##$PVM_EncSpectroscopy=", ...
    "PVM_NRepetitions=","PVM_NAverages=","PVM_EchoTime=",["##$PVM_ArrayPhase=( 2 )" + sprintf('\n')]];


GUIparamnames=["acq_time","spectralwidth","sequence","np","nucleus",...
    "resfreq","ppm_ref","tr","gain","voxs","acqtype","nrep","nav","te","rxarrayphases"];
numerical=[true, true, false, true, false, true, true, true, true,false, ...
    false, true, true,true,true];
arraylistbol=[false, false, false, false, false, false, false, false, ...
    false, true, false, false, false,false,true];

for par=1:length(BRUKERparamnames)
    startind=strfind(methodfile,BRUKERparamnames(par));
    if isempty(startind)
    else
    startind=startind(1);
    param=methodfile(startind+strlength(BRUKERparamnames(par)));
    
    if arraylistbol(par)
        kparam=1;
        while methodfile(startind+strlength(BRUKERparamnames(par))+kparam+1)~='#'
            param=[param methodfile(startind+strlength(BRUKERparamnames(par))+kparam)];
            kparam=kparam+1;
        end
        param=strsplit(param);
    else
        kparam=1;
%         disp(methodfile(startind+strlength(BRUKERparamnames(par))+kparam))
        while methodfile(startind+strlength(BRUKERparamnames(par))+kparam)~=['#',' ','$',sprintf('\n')]
            param=[param methodfile(startind+strlength(BRUKERparamnames(par))+kparam)];
            kparam=kparam+1;
        end
    end
    if numerical(par)
        param=str2double(param);
    end
    varname=matlab.lang.makeValidName(GUIparamnames(par));
    study.(varname{1})=param; 
    end 
end 

%adjustements
study.acq_time=study.acq_time*10^-3; % in seconds
study.np=study.np*2; %real/imag
study.tr=study.tr/1000; %in seconds 
if study.acqtype=='Yes'
    study.acqtype='MRS';

else 
    disp("can't open this type of data")
end
study.format='Matlab';
study.params.sfrq=study.resfreq;
study.params.arraydim=1; %MRS
study.params.np=study.np;
study.params.sw=study.spectralwidth;
study.params.tr=study.tr;
study=rmfield(study,'tr');
study.params.te=study.te; 
stuy=rmfield(study,'te');
study.params.gain=study.gain;
study=rmfield(study,'gain');
study.params.vox1=str2double(study.voxs(1));
study.params.vox2=str2double(study.voxs(2));
study.params.vox3=str2double(study.voxs(3));
study=rmfield(study,'voxs');

%% Group Delay - acqus file
if isfile(folderexp +num2str(expnb)+ "\acqus")
    acqusfile = fileread(folderexp +num2str(expnb)+ "\acqus");
    startind_grpdly=strfind(acqusfile,"##$GRPDLY= ");
    startind_grpdly=startind_grpdly(1);
    grpdly=acqusfile(startind_grpdly+strlength("##$GRPDLY= "));
    kgrp=1;
    while acqusfile(startind_grpdly+strlength("##$GRPDLY= ")+kgrp+1)~='#'
        grpdly=[grpdly acqusfile(startind_grpdly+strlength("##$GRPDLY= ")+kgrp)];
        kgrp=kgrp+1;
    end
    grpdly=str2double(grpdly);

elseif isfile(folderexp +num2str(expnb)+ "\acqp")
    acqpfile = fileread(folderexp +num2str(expnb)+ "\acqp");
    startind_grpdly=strfind(acqpfile,"##$ACQ_RxFilterInfo=( 2 )" + sprintf('\n') + "(");
    startind_grpdly=startind_grpdly(1);
    grpdly=acqpfile(startind_grpdly+strlength("##$ACQ_RxFilterInfo=( 2 )" + sprintf('\n') + "("));
    kgrp=1;
    while acqpfile(startind_grpdly+strlength("##$ACQ_RxFilterInfo=( 2 )" + sprintf('\n') + "(")+kgrp+1)~='#'
        grpdly=[grpdly acqpfile(startind_grpdly+strlength("##$ACQ_RxFilterInfo=( 2 )" + sprintf('\n') + "(")+kgrp)];
        kgrp=kgrp+1;
    end
    grpdly=strsplit(grpdly);
    grpdly=str2double(grpdly{1});
end
study.params.grpdly=grpdly;

%% ADCoverflow? - acqp file

acqpfile = fileread(folderexp +num2str(expnb)+ "\acqp");
startind_overflow=strfind(acqpfile,"##$ACQ_adc_overflow=( 2 )" + sprintf('\n')); %2 for 2 receiver channels
if isempty(startind_overflow)==1
    study.params.adcoverflow='No';
else
startind_overflow=startind_overflow(1);
overflow=acqpfile(startind_overflow+strlength("##$ACQ_adc_overflow=( 2 )" + sprintf('\n')));
kov=1;
while acqpfile(startind_overflow+strlength("##$ACQ_adc_overflow=( 2 )" + sprintf('\n')) +kov+1)~='#'
    overflow=[overflow acqpfile(startind_overflow+strlength("##$ACQ_adc_overflow=( 2 )" + sprintf('\n'))+kov)];
    kov=kov+1;
end
overflow=strsplit(overflow);
if strcmp(overflow{1},'Yes') | strcmp(overflow{2},'Yes')
    study.params.adcoverflow='Yes';
    f=msgbox(['ADC overflow during acquisition E' num2str(expnb)]);
else
    study.params.adcoverflow='No';
end 
end 


%% 1) Read data ser/fid

disp(study.sequence)

%job0 as ser/fid
if isfile(folderexp +num2str(expnb)+ "\ser") % case PV360 v1 where there is a ser file in the main folder (2D acq opened in topspin)
    disp ("2D ser, " + "nav=" + num2str(study.nav) + " nrep=" + num2str(study.nrep) + " saved")
    fileid=fopen(folderexp +num2str(expnb)+ "\ser",'r','ieee-le'); %read binary format
    if fileid == -1
        disp('Cannot open file');
        return
    end
    buffer=fread(fileid,'int32'); %native format v1 for MRS is int32  
elseif isfile(folderexp +num2str(expnb)+ "\fid") % case PV360 v1 where there is a fid file in the main folder (2D acq NOT opened in topspin or 1D acq opened or not in topspin)
    if study.nrep==1
        dimtype="1D";
    else
        dimtype="2D";
    end 
    disp (dimtype + " fid, " + "nav=" + num2str(study.nav) + " nrep=" + num2str(study.nrep)+ " saved")
    fileid=fopen(folderexp +num2str(expnb)+ "\fid",'r','ieee-le'); %read binary format
    if fileid == -1
        disp('Cannot open file');
        return
    end
    buffer=fread(fileid,'int32'); %native format v1 for MRS is int32 
elseif isfile(folderexp +num2str(expnb)+ "\pv2tsdata\1\ser") % case PV360 v3 where there is a ser file in the subfolder pv2tsdata\1\ (2D acq opened in topspin)
    disp ("2D ser, " + "nav=" + num2str(study.nav) + " nrep=" + num2str(study.nrep) + " saved")
    fileid=fopen(folderexp +num2str(expnb)+ "\pv2tsdata\1\ser",'r','ieee-le'); %read binary format
    if fileid == -1
        disp('Cannot open file');
        return
    end
    buffer=fread(fileid,'double'); %native format v3 for MRS is double 
elseif  isfile(folderexp +num2str(expnb)+ "\pv2tsdata\1\fid") % case PV360 v3 where there is a fid file in the subfolder pv2tsdata\1\ (1D acq opened in topspin)
    disp ("1D fid, " + "nav=" + num2str(study.nav) + " nrep=" + num2str(study.nrep) + " saved")
    fileid=fopen(folderexp +num2str(expnb)+ "\pv2tsdata\1\fid",'r','ieee-le'); %read binary format
    if fileid == -1
        disp('Cannot open file');
        return
    end
    buffer=fread(fileid,'double'); %native format v3 for MRS is double  
    %% THE FOLLOWING CONDITION HAS NOT BEEN TESTED 27/5/24
elseif  isfile(folderexp +num2str(expnb)+ "\pdata\1\fid_proc.64") % case PV360 v3 where there is a fid file in the subfolder pdata\1\ (1D or 2D acq NOT opened in topspin)
    disp ("1D fid, " + "nav=" + num2str(study.nav) + " nrep=" + num2str(study.nrep) + " saved")
    fileid=fopen(folderexp +num2str(expnb)+ "\pdata\1\fid_proc.64",'r','ieee-le'); %read binary format
    if fileid == -1
        disp('Cannot open file');
        return
    end
    buffer=fread(fileid,'double'); %native format v3 for MRS is double  
end  

nbptsfid=length(buffer)/(2*study.nrep);
buffer_ser=zeros(study.nrep,nbptsfid*2);
for rep=1:study.nrep
    buffer_ser(rep,:)=buffer((rep-1)*(nbptsfid*2)+1:rep*(nbptsfid*2))';
end
ser_c=buffer_ser(:,1:2:end)+1i*buffer_ser(:,2:2:end);
fclose(fileid);

%offset 
if study.ppm_ref ~= 0
    offset_hz=study.ppm_ref*study.resfreq;
    dw=1/study.spectralwidth;
    t=[0:dw:(study.np/2-1)*dw];
    tmat=repmat(t,study.nrep,1);
    ser_c_shift=ser_c.*exp(1i.*(2*pi*offset_hz).*tmat);
    study.data.real(:,1,:)=real(ser_c_shift);
    study.data.imag(:,1,:)=-imag(ser_c_shift); %flips the spectrum
else
    study.data.real(:,1,:)=real(ser_c);
    study.data.imag(:,1,:)=-imag(ser_c); %flips the spectrum
end 

%filename and liststring
if study.nrep>1
    filename=['Bruker_' date(1:10)  '_'  num2str(expnb) timeacq '_ser.mat'];
    %nt
    study.params.nt=study.nrep*study.nav;

else
    filename=['Bruker_' date(1:10)  '_'  num2str(expnb) timeacq '_fid.mat'];
    %nt
    study.params.nt=study.nav;
end 

%normalize: xNA/RG/voxelsize
voxvol=study.params.vox1.*study.params.vox2.*study.params.vox3;
study.data.real=study.data.real.*study.nav./study.params.gain./voxvol;
study.data.imag=study.data.imag.*study.nav./study.params.gain./voxvol;

study.filename=filename;
study.liststring=char(folderexp + study.filename);

%multiplicity
study.multiplicity=study.nrep;

%process
study.process.lsfid=round(grpdly)+1;
study.process.apodizefct='exponential'; %default
study.process.apodparam1= zeros(1,study.nrep);%default
study.process.apodparam2=zeros(1,study.nrep);%default
study.process.transfsize=0; %default
study.process.appltoarray1=0;%default
study.process.appltoarray2=0;%default
study.process.phasecorr0=zeros(1,study.nrep);%default
study.process.phasecorr1=zeros(1,study.nrep);%default
study.process.DCoffset=0; %default


%% 2) read data rawdatajob0 - moved from the separate code Read_bruker_rawdatajob0.... to here - JM 27/5/24

nloops=2;
scalingloops=[100,100];
rxarrayphases=study.rxarrayphases;

pathrawdata=folderexp +num2str(expnb)+ "\rawdata.job0";

fileid=fopen(pathrawdata,'r','ieee-le'); %read binary format
if fileid == -1
    disp('Cannot open file');
    return
end

bufferrawdata=fread(fileid,'int32'); 
nbptsfid=length(bufferrawdata)/(2*study.nav*study.nrep*nloops); %2 is for real, imag
fclose(fileid);

fid_reorganized=zeros(nloops,study.nav*study.nrep,nbptsfid);

for avrep=1:study.nav*study.nrep
    for coil=1:nloops
        fid_reorganized(coil,avrep,:)=bufferrawdata((avrep-1)*nbptsfid*2*nloops+(coil-1)*nbptsfid*2+1:2:(avrep-1)*nbptsfid*2*nloops+(coil-1)*nbptsfid*2+nbptsfid*2)+1i*bufferrawdata((avrep-1)*nbptsfid*2*nloops+(coil-1)*nbptsfid*2+2:2:(avrep-1)*nbptsfid*2*nloops+(coil-1)*nbptsfid*2+nbptsfid*2);
    end 
end 

%
fid_comb_re=zeros(study.nav*study.nrep,nbptsfid);
fid_comb_im=zeros(study.nav*study.nrep,nbptsfid);

fid_reorganized_rephased_re=zeros(nloops,study.nav*study.nrep,nbptsfid);
fid_reorganized_rephased_im=zeros(nloops,study.nav*study.nrep,nbptsfid);

for coil=1:nloops
    fid_comb_re=fid_comb_re + ...
                    (squeeze(real(fid_reorganized(coil,:,:))).*cos(rxarrayphases(coil)/180*pi) ... 
                    -squeeze(imag(fid_reorganized(coil,:,:))).*sin(rxarrayphases(coil)/180*pi)).*scalingloops(coil);

    fid_reorganized_rephased_re(coil,:,:)=(squeeze(real(fid_reorganized(coil,:,:))).*cos(rxarrayphases(coil)/180*pi) ... 
                    -squeeze(imag(fid_reorganized(coil,:,:))).*sin(rxarrayphases(coil)/180*pi)).*scalingloops(coil);

    fid_comb_im= fid_comb_im+ ...
                    (squeeze(real(fid_reorganized(coil,:,:))).*sin(rxarrayphases(coil)/180*pi) ... 
                    +squeeze(imag(fid_reorganized(coil,:,:))).*cos(rxarrayphases(coil)/180*pi)).*scalingloops(coil);

    fid_reorganized_rephased_im(coil,:,:)=(squeeze(real(fid_reorganized(coil,:,:))).*sin(rxarrayphases(coil)/180*pi) ... 
                    +squeeze(imag(fid_reorganized(coil,:,:))).*cos(rxarrayphases(coil)/180*pi)).*scalingloops(coil);
end 


%fid before coil combination but with rephasing
fid_reorganized_rephased=fid_reorganized_rephased_re-1i*fid_reorganized_rephased_im;

%fid after coil combination
fid_comb=fid_comb_re-1i*fid_comb_im; 


study.datajob0.real=zeros(size(fid_comb,1),1,size(fid_comb,2));
study.datajob0.real(:,1,:)=real(fid_comb);

study.datajob0.imag=zeros(size(fid_comb,1),1,size(fid_comb,2));
study.datajob0.imag(:,1,:)=imag(fid_comb);

study.datajob0.real=study.datajob0.real./study.params.gain./voxvol.*study.nav/(nloops*scalingloops(1));
study.datajob0.imag=study.datajob0.imag./study.params.gain./voxvol.*study.nav/(nloops*scalingloops(1));

study=rmfield(study,'nav');
study=rmfield(study,'nrep');

save(folderexp + filename,'study')


%simple display of what is actually saved (study.data...)
figure; 
grpdly=round(study.params.grpdly)+1;
savedrawdatajob0=study.datajob0.real+1i.*study.datajob0.imag;
savedrawdatajob0=[savedrawdatajob0(:,grpdly:end),savedrawdatajob0(:,1:grpdly-1)];
ftsaved=fftshift(fft(savedrawdatajob0,[],2),2);
plot(real(sum(ftsaved).*exp(-i*0.5)),'DisplayName','rawdatajob0')
title(num2str(expnb)) 
legend()

figure;
savedjob0=study.data.real+1i.*study.data.imag;
savedjob0=[savedjob0(:,grpdly:end),savedjob0(:,1:grpdly-1)];
ftsavedjob0=fftshift(fft(savedjob0,[],2),2);
if size(study.data.real,1)>1
    plot(real(sum(ftsavedjob0)),'DisplayName','ser/fid')
else
    plot(real(ftsavedjob0),'DisplayName','ser/fid')
end
legend()
title(num2str(expnb)) 




%% Refscan

if isfile(folderexp + num2str(expnb)+ "\fid.refscan") %for 360v1

    %literal names = navrefscan,
    BRUKERparamnames_refscan=["PVM_RefScanNA=","##$PVM_RefScanRG="];
    localparamnames_refscan=["NArefscan","RGrefscan"];
    numerical=[true, true];

    for par=1:length(BRUKERparamnames_refscan)
        startind=strfind(methodfile,BRUKERparamnames_refscan(par));
        startind=startind(1);
        param=methodfile(startind+strlength(BRUKERparamnames_refscan(par)));
        kparam=1;
        while methodfile(startind+strlength(BRUKERparamnames_refscan(par))+kparam+1)~=['#',' ','$']
            param=[param methodfile(startind+strlength(BRUKERparamnames_refscan(par))+kparam)];
            kparam=kparam+1;
        end
        param=str2double(param);
        varname=matlab.lang.makeValidName(localparamnames_refscan(par));
        refscan.(varname{1})=param;    
    end 

    %store data
    disp('refscan saved')
    fileid=fopen(folderexp + num2str(expnb)+ "\fid.refscan",'r','ieee-le'); %read binary format
    if fileid == -1
        disp('Cannot open file');
        return
    end
    buffer=fread(fileid,'int32'); 
    nbptsfid=length(buffer)/2;
    buffer_c_ref=buffer(1:2:end)+1i*buffer(2:2:end); 
    fclose(fileid);

    study.data.real=zeros(1,1,study.np/2);
    study.data.imag=zeros(1,1,study.np/2);
    study.data.real(1,1,:)=real(buffer_c_ref);
    study.data.imag(1,1,:)=-imag(buffer_c_ref); %flips the spectrum

    %filename and liststring
    filename=['Bruker_' date(1:10)  '_'  num2str(expnb)  timeacq '_fidrefscan.mat'];
    study.filename=filename;
    study.liststring=char(folderexp + study.filename);


    %normalize: xNA/RG/voxelsize 
    study.data.real=study.data.real.*refscan.NArefscan./refscan.RGrefscan./voxvol;
    study.data.imag=study.data.imag.*refscan.NArefscan./refscan.RGrefscan./voxvol;

    %multiplicity
    study.multiplicity=1; %1D

    %process
    study.process.lsfid=round(grpdly)+1;
    study.process.apodizefct='exponential'; %default
    study.process.apodparam1= 0;%default
    study.process.apodparam2=0;%default
    study.process.transfsize=0; %default
    study.process.appltoarray1=0;%default
    study.process.appltoarray2=0;%default
    study.process.phasecorr0=0;%default
    study.process.phasecorr1=0;%default
    study.process.DCoffset=0; %default
    study.process.B0=zeros(1,study.np/2);

    %nt
    study.params.nt=refscan.NArefscan;

    save(folderexp + filename,'study');

    %simple display
    fidrefscan=squeeze(study.data.real)'+1i*squeeze(study.data.imag)'; 
    fidrefscan_sum=sum(fidrefscan,1);
    fidrefscan_sum=[fidrefscan_sum(round(grpdly)+1:end), zeros(1,round(grpdly))];
    ftrefscan=fftshift(fft(fidrefscan_sum./study.params.nt.*2.8e-4,[],2),2);
    figure; 
    plot(real(ftrefscan.*exp(i*3.14)),'DisplayName','refscan')
    legend()
    title(num2str(expnb))

end
   
x1=sum(real(fftshift(fft(squeeze(fid_reorganized_rephased(1,:,:)),[],2),2)));
x2=sum(real(fftshift(fft(squeeze(fid_reorganized_rephased(2,:,:)),[],2),2)));

x1s=[x1(77:4096),x1(1:76)];
x2s=[x2(77:4096),x2(1:76)];
end

