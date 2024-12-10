%io_loadspec_bruk_PV6.m
%Diana Rotaru, KCL, 2018, updated in 2019 for MEGAPRESS and HERMES data
%based on io_loadspec_bruk.m
%Chathura Kumaragamage, McGill University 2016.
%Jamie Near, McGill University 2016.
%
% THIS VERSION OF THE FUNCTION WORKS FOR BRUKER PV6!!!
% Version before PV6 output the ref scan named fid.ref whereas PV6 uses
% fid.refscan
%
% USAGE:
% [out,ref]=io_loadspec_bruk_PV6(filename);
%
% DESCRIPTION:
% Reads in Bruker MRS data (fid.raw, fid.refscan).
%
% op_loadspec_bruk_PV6 outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.
%
% INPUTS:
% inDir   = Path to the scan number directory that contains the 'pdata' folder.
%
% OUTPUTS:
% out = Input dataset in FID-A structure format.
% ref = The Reference scan data (navigator echoes) in FID-A structure
%       format, if applicable.



function [out,ref]=io_loadspec_bruk_PV6(inDir)

%%%%%%%%chathu mod starts

%get the original number of Repetitions (all repatitions from scanner is generated as
%a 1D vector. Need to split before further processing
averages=1; %Because Bruker does the averaging online.
ADC_OFFSET=70;      %(offset between ADC on and ADC acquire, empirically determined)
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$PVM_DigNp=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'$PVM_DigNp=');
end
equals_index=findstr(line,'=');
rawDataPoints=line(equals_index+1:end);
rawDataPoints=str2double(rawDataPoints);
fclose(method_fid);

%NOW LOAD IN THE RAW DATA.  FIRST TRY USING THE FID.RAW FILE.  IF THAT DOES
%NOT WORK, THEN USE THE REGULAR FID FID.
%%% Diana: get the method name and based on it distinguish non-editing vs.
%%% editing sequences

%Now get the sequence
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$Method=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'$Method=');
end
equals_index=findstr(line,'=');
sequence=line(equals_index+1:end);
sequence=strtrim(sequence);
fclose(method_fid);

if string(extractBetween(sequence,'<','>')) == 'User:mpress'
    sequence = 'MEGAPRESS/HERMES';
    subspecs=4;
    rawSubspecs=4;
else
    sequence = 'PRESS';
    subspecs=0;
    rawSubspecs=0;
end
%%% End Diana

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$PVM_NAverages=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'$PVM_NAverages=');
end
equals_index=findstr(line,'=');
rawAverages=line(equals_index+1:end);
rawAverages=str2double(rawAverages);
fclose(method_fid);

if strcmp(sequence,'PRESS')
    try
        fid_data=fread(fopen([inDir '/fid.raw']),'int');
        real_fid = fid_data(2:2:length(fid_data));
        imag_fid = fid_data(1:2:length(fid_data));
        fids_raw=real_fid-1i*imag_fid;
        
        fids_raw=reshape(fids_raw,[],rawAverages);
        fids_trunc=fids_raw(ADC_OFFSET:end,:);
        %         fids_trunc=fids_trunc';
        fids=padarray(fids_trunc, [ADC_OFFSET-1,0],'post');
        specs=fftshift(ifft(fids,[],1),1);
        sz=size(specs);
        out.flags.averaged=1;
        %specify the dims
        dims.t=1;
        dims.coils=0;
        dims.averages=0;
        dims.subSpecs=0;
        dims.extras=0;
    catch
        %NOW TRY USING FID FILE
        display('WARNING: /fid.raw not found. Using /fid ....')
        %specs=fread(fopen([inDir '/pdata/1/2dseq']),'int');
        if exist(([inDir '/fid']))
        fid_data=fread(fopen([inDir '/fid']),'int');
        real_fid = fid_data(2:2:length(fid_data));
        imag_fid = fid_data(1:2:length(fid_data));
        fids_raw=real_fid-1i*imag_fid;
        
        rawAverages=size(real_fid,1)./rawDataPoints;
        
        if mod(size(real_fid,1),rawDataPoints) ~=0
            display 'number of repetitions cannot be accurately found';
        end
        fids_raw=reshape(fids_raw,rawDataPoints,[]);
        
        fids_trunc=fids_raw(ADC_OFFSET-1:end,:);
        fids=padarray(fids_trunc, [ADC_OFFSET-2,0],'post');
        %convert back to time domain
        %if the length of Fids is odd, then you have to do a circshift of one to
        %make sure that you don't introduce a small frequency shift into the fids
        %vector.
        specs=fftshift(ifft(fids,[],1),1);
        
        %calculate the size;
        sz=size(specs);
        out.flags.averaged=1;
        %specify the dims
        dims.t=1;
        dims.coils=0;
        dims.averages=0;
        dims.subSpecs=0;
        dims.extras=0;
        else
        display('WARNING: /fid not found. Using /fid .... Aborting!')   
        end
        
    end
    %NOW TRY LOADING IN THE REFERENCE SCAN DATA (IF IT EXISTS)
    isRef=true;
    if exist([inDir '/fid.refscan'])
        isRef=true;
        rawAverages=1;
        ref_data=fread(fopen([inDir '/fid.refscan']),'int');
        real_ref = ref_data(2:2:length(ref_data));
        imag_ref = ref_data(1:2:length(ref_data));
        fids_ref=real_ref-1i*imag_ref;
        
        fids_ref=reshape(fids_ref,[],rawAverages);
        fids_ref_trunc=fids_ref(ADC_OFFSET:end,:);
        %         fids_trunc=fids_trunc';
        reffids=padarray(fids_ref_trunc, [ADC_OFFSET-1,0],'post');
        refspecs=fftshift(ifft(reffids,[],1),1);
        szreffids=size(refspecs);
        out.flags.averaged=0;
        
        
        %specify the dims
        refdims.t=1;
        refdims.coils=0;
        refdims.averages=0;
        refdims.subSpecs=0;
        refdims.extras=0;
        %%% Diana
        if size(sequence) == size('PRESS')
            addedrcvrs=1;
        else
            addedrcvrs=0;
        end
        %%% End Diana
    end
else
    
    % Get access to rawdata
    fid_data = fopen(['rawdata.job0'],'r');
    fids_raw = fread(fid_data,inf,'int32');
    fclose(fid_data);
    
    % Determine whether it is MEGAPRESS or HERMES sequence
    % Compare the editing frequencies and if any 2 are identical it is
    % MEGAPRESS, otherwise HERMES
    
    method_fid=fopen([inDir '/method']);
    line=fgets(method_fid);
    index=findstr(line,'$Ed1Freqppm=( 2 )');
    while isempty(index)
        line=fgets(method_fid);
        index=findstr(line,'$Ed1Freqppm=( 2 )');
    end
    Ed1Freqppm = fgets(method_fid,10);
    
    method_fid=fopen([inDir '/method']);
    line=fgets(method_fid);
    index=findstr(line,'$Ed2Freqppm=( 2 )');
    while isempty(index)
        line=fgets(method_fid);
        index=findstr(line,'$Ed2Freqppm=( 2 )');
    end
    Ed2Freqppm = fgets(method_fid,10);
    
    method_fid=fopen([inDir '/method']);
    line=fgets(method_fid);
    index=findstr(line,'$Ed3Freqppm=( 2 )');
    while isempty(index)
        line=fgets(method_fid);
        index=findstr(line,'$Ed3Freqppm=( 2 )');
    end
    Ed3Freqppm = fgets(method_fid,10);
    
    method_fid=fopen([inDir '/method']);
    line=fgets(method_fid);
    index=findstr(line,'$Ed4Freqppm=( 2 )');
    while isempty(index)
        line=fgets(method_fid);
        index=findstr(line,'$Ed4Freqppm=( 2 )');
    end
    Ed4Freqppm = fgets(method_fid,10);
    
    % Determine the number of data points
    method_fid=fopen([inDir '/method']);
    line=fgets(method_fid);
    index=findstr(line,'$PVM_SpecMatrix=( 1 )');
    while isempty(index)
        line=fgets(method_fid);
        index=findstr(line,'$PVM_SpecMatrix=( 1 )');
    end
    PVM_SpecMatrix = fgets(method_fid,4);
    PVM_SpecMatrix =str2num(PVM_SpecMatrix);
    
    % Determine the number of averages
    method_fid=fopen([inDir '/method']);
    line=fgets(method_fid);
    index=findstr(line,'$PVM_NAverages=');
    while isempty(index)
        line=fgets(method_fid);
        index=findstr(line,'$PVM_NAverages=');
    end
    equals_index=findstr(line,'=');
    PVM_NAverages =line(equals_index+1:end);
    PVM_NAverages =str2num(PVM_NAverages);
    
    if ~strcmp(Ed1Freqppm,Ed2Freqppm) && ~strcmp(Ed1Freqppm,Ed3Freqppm) && ~strcmp(Ed1Freqppm,Ed4Freqppm)
        seq_name = 'HERMES';
        
        % Reshape the 1D vector of data points into a 5D array
        fids_raw=reshape(fids_raw, [2 PVM_SpecMatrix 4 4 PVM_NAverages/4]);

        % Get the real an imaginary parts for each data point and reduce
        % the dimension by 1
        fids_raw_cmplx = squeeze(fids_raw(1,:,:,:,:)) + 1i*squeeze(fids_raw(2,:,:,:,:));
        
    else
        seq_name = 'MEGAPRESS';
        % Reshape the 1D vector of data points into a 5D array
        fids_raw=reshape(fids_raw, [2 PVM_SpecMatrix 4 2 PVM_NAverages/2]);
        % the reshape below doesn't produce better results and should be
        % accompanied by the lines after coilc comb in
        % megapressbrukerproc_auto_updated.m
        %fids_raw=reshape(fids_raw, [2 PVM_SpecMatrix 4 4 PVM_NAverages/4]);
        % Get the real an imaginary parts for each data point and reduce
        % the dimension by 1
        fids_raw_cmplx = squeeze(fids_raw(1,:,:,:,:)) + 1i*squeeze(fids_raw(2,:,:,:,:));
    end
    
    
    % Determine the number of points corresponding to a [group delay of
    % the filter], shift the data to left and null at same number of
    % points at the end of the data set
    method_fid=fopen([inDir '/method']);
    line=fgets(method_fid);
    index=findstr(line,'$PVM_DigShift=');
    while isempty(index)
        line=fgets(method_fid);
        index=findstr(line,'$PVM_DigShift=');
    end
    equals_index=findstr(line,'=');
    PVM_DigShift=line(equals_index+1:end);
    PVM_DigShift=str2num(PVM_DigShift);
    
    fids_raw_shifted=circshift(fids_raw_cmplx,-PVM_DigShift);
    fids_raw_shifted(end-PVM_DigShift+1:end,:)=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % FOR TESTING REASOSNS, the next lines show the accumulated averages
%     % per coil with the FIDs filtered using an exponential, the ppm
%     % scale calculated and specs plotted.
%     
%     % For security reasons, all averages will be accumulated and the
%     % signals per coil and per sub-acquisition will be shown
%     fids_raw_cmplx_mean = squeeze(mean(fids_raw_cmplx,4));
%     fids_raw_shifted_mean = squeeze(mean(fids_raw_shifted,4));
%     
%     dataset = extractBetween(inDir, 'Animal-data/','/');
%     dataset = strrep(dataset, '_','-');
%     
%     figure
%     subplot(221)
%     plot(real(fids_raw_cmplx_mean(:,1,2)));
%     subplot(222)
%     plot(real(fids_raw_cmplx_mean(:,2,2)));
%     subplot(223)
%     plot(real(fids_raw_cmplx_mean(:,3,2)));
%     subplot(224)
%     plot(real(fids_raw_cmplx_mean(:,4,2)));
%     %title([dataset, '(raw): accumulated averages - ch2 - acq A-D']);
%     
%     figure
%     subplot(221)
%     plot(real(fids_raw_shifted_mean(:,1,2)));
%     subplot(222)
%     plot(real(fids_raw_shifted_mean(:,2,2)));
%     subplot(223)
%     plot(real(fids_raw_shifted_mean(:,3,2)));
%     subplot(224)
%     plot(real(fids_raw_shifted_mean(:,4,2)));
%     %title([dataset, '(shifted): accumulated averages - ch2 - acq A-D']);
%     
%     % Create an exponential function with a fast decay and apply it
%     % to all FIDs
%     x=1:PVM_SpecMatrix;
%     y=exp(-x/100);
%     figure,
%     plot(y)
%     
%     if strcmp(seq_name,'HERMES')
%         
%         
%         for q=1:4
%             
%             fids_raw_shifted_mean(:,1,q)= y'.*fids_raw_shifted_mean(:,1,q);
%             fids_raw_shifted_mean(:,2,q)= y'.*fids_raw_shifted_mean(:,2,q);
%             fids_raw_shifted_mean(:,3,q)= y'.*fids_raw_shifted_mean(:,3,q);
%             fids_raw_shifted_mean(:,4,q)= y'.*fids_raw_shifted_mean(:,4,q);
%             for r=1:PVM_NAverages/4
%                 fids_raw_shifted_filt(:,1,q,r)= y'.*fids_raw_shifted(:,1,q,r);
%                 fids_raw_shifted_filt(:,2,q,r)= y'.*fids_raw_shifted(:,2,q,r);
%                 fids_raw_shifted_filt(:,3,q,r)= y'.*fids_raw_shifted(:,3,q,r);
%                 fids_raw_shifted_filt(:,4,q,r)= y'.*fids_raw_shifted(:,4,q,r);
%                 
%             end
%         end
%     else
%         
%         for q=1:2
%             fids_raw_shifted_mean(:,1,q)= y'.*fids_raw_shifted_mean(:,1,q);
%             fids_raw_shifted_mean(:,2,q)= y'.*fids_raw_shifted_mean(:,2,q);
%             fids_raw_shifted_mean(:,3,q)= y'.*fids_raw_shifted_mean(:,3,q);
%             fids_raw_shifted_mean(:,4,q)= y'.*fids_raw_shifted_mean(:,4,q);
%             
%             for r=1:PVM_NAverages/2
%                 fids_raw_shifted_filt(:,1,q,r)= y'.*fids_raw_shifted(:,1,q,r);
%                 fids_raw_shifted_filt(:,2,q,r)= y'.*fids_raw_shifted(:,2,q,r);
%                 fids_raw_shifted_filt(:,3,q,r)= y'.*fids_raw_shifted(:,3,q,r);
%                 fids_raw_shifted_filt(:,4,q,r)= y'.*fids_raw_shifted(:,4,q,r);
%                 
%             end
%         end
%     end
%     
%     figure
%     title([dataset, '(shifted,fitered-exp): accumulated averages - ch2 - acq A-D']);
%     subplot(221)
%     plot(real(fids_raw_shifted_mean(:,1,2)))
%     subplot(222)
%     plot(real(fids_raw_shifted_mean(:,2,2)))
%     subplot(223)
%     plot(real(fids_raw_shifted_mean(:,3,2)))
%     subplot(224)
%     plot(real(fids_raw_shifted_mean(:,4,2)))
%     
    % Determine the spectral width in Hz and the transmitter frequency and
    % calculate the spectaral width in ppm
    method_fid=fopen([inDir '/method']);
    line=fgets(method_fid);
    index=findstr(line,'$PVM_DigSw=');
    while isempty(index)
        line=fgets(method_fid);
        index=findstr(line,'$PVM_DigSw=');
    end
    equals_index=findstr(line,'=');
    spectralwidth=line(equals_index+1:end);
    spectralwidth=str2double(spectralwidth);
    fclose(method_fid);
    
    acqp_fid=fopen([inDir '/acqp']);
    line=fgets(acqp_fid);
    index=findstr(line,'$BF1=');
    while isempty(index)
        line=fgets(acqp_fid);
        index=findstr(line,'$BF1=');
    end
    equals_index=findstr(line,'=');
    txfrq=line(equals_index+1:end);
    txfrq=str2double(txfrq);
    txfrq=txfrq*1e6;
    fclose(acqp_fid);
    
    %B0
    Bo=txfrq/42577000;
    
    %Spectral width in PPM
    spectralwidthppm=spectralwidth/(txfrq/1e6);
    ppm=[4.65-(spectralwidthppm/2):spectralwidthppm/(PVM_SpecMatrix-1):4.65+(spectralwidthppm/2)];
    
%     % Plot the specs for each  sub-acqusition, ch2 only
%     figure
%     subplot(221)
%     plot(ppm,-real(fftshift(fft(fids_raw_shifted_mean(:,1,2)))));
%     set(gca,'XDir','reverse');
%     axis([0 4.5 -1000 4000]);
%     
%     subplot(222)
%     plot(ppm,-real(fftshift(fft(fids_raw_shifted_mean(:,2,2)))));
%     set(gca,'XDir','reverse');
%     axis([0 4.5 -1000 4000]);
%     subplot(223)
%     plot(ppm,-real(fftshift(fft(fids_raw_shifted_mean(:,3,2)))));
%     set(gca,'XDir','reverse');
%     axis([0 4.5 -1000 4000]);
%     subplot(224)
%     plot(ppm,-real(fftshift(fft(fids_raw_shifted_mean(:,4,2)))));
%     set(gca,'XDir','reverse');
%     axis([0 4.5 -1000 4000]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fids = fids_raw_shifted;
    % Update the structure fields the same way as in
    % io_loadspec_orchestra_hermes in order to use the same functions for coil
    % combination and averages accumulation
    fids = permute(fids,[1 2 4 3]);
    
    x=1:PVM_SpecMatrix;
    y=exp(-x/100);
    figure,plot(y)
    fids= y'.*fids(:,:,:,:);
    
    specs=fftshift(ifft(fids,[],1),1);
    
    sz=size(specs);
    out.flags.averaged=0;
    %specify the dims
    dims.t=1;
    dims.coils=2;
    dims.averages=3;
    dims.subSpecs=4;
    dims.extras=0;
    
    %NOW TRY LOADING IN THE REFERENCE SCAN DATA (IF IT EXISTS)
    % The unprocessed water data can be found in the method file, between
    % '##$PVM_RefScan=( 4, 4096 )'and '##$PVM_RefScanJob=0'
    % First, these two lines must be found and stored in different
    % variables, then they must be used to select the line between them and
    % save those as an array of numbers
    
    method_fid=fopen([inDir '/method']);
    line=fgets(method_fid);
    index=findstr(line,'##$PVM_RefScan=');
    while isempty(index)
        line=fgets(method_fid);
        index=findstr(line,'##$PVM_RefScan=');
    end
    
    % Save the start line for the unprocessed water data selection
    start_line = line(index:end-1);
    
    obraket_index=findstr(line,'( ');
    cbraket_index=findstr(line,' )');
    coma_index=findstr(line,',');
    refScanCoilsNum=line(obraket_index+1:coma_index-1);
    refScanCoilsNum=refScanCoilsNum(find(~isspace(refScanCoilsNum)));
    refScanCoilsNum=str2num(refScanCoilsNum);
    refScanPointsNum=line(coma_index+1:cbraket_index-1);
    refScanPointsNum=refScanPointsNum(find(~isspace(refScanPointsNum)));
    refScanPointsNum=str2num(refScanPointsNum);
    line=fgets(method_fid);
    
    method_fid=fopen([inDir '/method']);
    line=fgets(method_fid);
    index=findstr(line,'##$PVM_RefScanJob');
    while isempty(index)
        line=fgets(method_fid);
        index=findstr(line,'##$PVM_RefScanJob');
    end
    
    % Save the end line for the unprocessed water data selection
    end_line = line(index:end-1);
    
    % Extract data between start_line and end_line
    fid=fopen('method');  % open the file
    while ~feof(fid)  % loop to the end of file
        l=fgetl(fid);            % read a record of each line
        if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
            fids_ref=0;
            while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
                l=fgetl(fid);
                if fids_ref==0
                    fids_ref=l;
                else
                    fids_ref= horzcat(fids_ref,l);
                end
            end
            
        end
    end
    
    % Remove the last characters corresponding to the end_line
    fids_ref = (fids_ref(1:end-length(end_line)));
    % Convert variable to a number
    fids_ref = str2num(fids_ref);
    
    fclose(method_fid);
    
    % Reshape data as 4 (coils num.) x 2048 (data points num.)
    refscan = reshape(fids_ref, [], refScanCoilsNum);
    refscan_real = refscan(1:2:end,:);
    refscan_imag = refscan(2:2:end,:);
    
    if ~isequaln(size(refscan_real,1),size(refscan_imag,1))
        if size(refscan_real,1)>size(refscan_imag,1)
            refscan_real= refscan_real(size(refscan_imag,1),:);
        else 
            refscan_imag= refscan_imag(size(refscan_real,1),:);
        end
    end

    figure
    subplot(2,2,1)
    plot(refscan_real(:,1))
    hold on
    plot(refscan_imag(:,1))
    subplot(2,2,2)
    plot(refscan_real(:,2))
    hold on
    plot(refscan_imag(:,2))
    subplot(2,2,3)
    plot(refscan_real(:,3))
    hold on
    plot(refscan_imag(:,3))
    subplot(2,2,4)
    plot(refscan_real(:,4))
    hold on
    plot(refscan_imag(:,4))
    
    % Save complex data
    reffids = refscan_real(:,:)+1i*refscan_imag(:,:);
    
    % Remove digital filter delay points
    reffids=circshift(reffids,-PVM_DigShift);
    reffids(end-PVM_DigShift+1:end,:)=0;
    
    
    refspecs=fftshift(ifft(reffids,[],1),1);
    szreffids=size(refspecs);
    ref.flags.averaged=1;
    isRef=1;
    
    figure
    plot(ppm(1:size(reffids,1)),-real(fftshift(fft(reffids(:,1)))));
    set(gca,'XDir','reverse');
    
    %specify the dims
    refdims.t=1;
    refdims.coils=2;
    refdims.averages=0;
    refdims.subSpecs=0;
    refdims.extras=0;
    %%% Diana
    
    %%% End Diana
    end


%%%%%Chathu mod ends

%Now get some of the relevent spectral parameters

%First get the spectral width
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$PVM_DigSw=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'$PVM_DigSw=');
end
equals_index=findstr(line,'=');
spectralwidth=line(equals_index+1:end);
spectralwidth=str2double(spectralwidth);
fclose(method_fid);

%Now get the transmitter frequency
acqp_fid=fopen([inDir '/acqp']);
line=fgets(acqp_fid);
index=findstr(line,'$BF1=');
while isempty(index)
    line=fgets(acqp_fid);
    index=findstr(line,'$BF1=');
end
equals_index=findstr(line,'=');
txfrq=line(equals_index+1:end);
txfrq=str2double(txfrq);
txfrq=txfrq*1e6;
fclose(acqp_fid);

%B0
Bo=txfrq/42577000;

%Spectral width in PPM
spectralwidthppm=spectralwidth/(txfrq/1e6);

%Now get the TE
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$PVM_EchoTime=');
%%%%%CHATHU MOD for FIDS with Pulse acquire
loop_count=0;
TE_present=true;
while isempty(index) && TE_present
    line=fgets(method_fid);
    index=findstr(line,'$PVM_EchoTime=');
    loop_count=loop_count+1;
    if loop_count>10000
        TE_present=false;
    end
end

if TE_present
    equals_index=findstr(line,'=');
    te=line(equals_index+1:end);
    te=str2double(te);
else
    te=0;
end

%%%%%CHATHU MOD for FIDS with Pulse acquire END
fclose(method_fid);

%Now get the TR
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'$PVM_RepetitionTime=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'$PVM_RepetitionTime=');
end
equals_index=findstr(line,'=');
tr=line(equals_index+1:end);
tr=str2double(tr);
fclose(method_fid);

%calculate the ppm scale
ppm=[4.65+(spectralwidthppm/2):-spectralwidthppm/(length(specs)-1):4.65-(spectralwidthppm/2)];

%calculate the dwelltime:
dwelltime=1/spectralwidth;

%calculate the time scale
t=[0:dwelltime:(sz(1)-1)*dwelltime];

if isequaln(sequence,'PRESS')
    addedrcvrsrffids=1;
    addedrcvrsreffids=1;
else
    addedrcvrsrffids=0;
    addedrcvrsreffids=0;
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
out.averages=rawAverages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=sequence;
out.te=te;
out.tr=tr;
out.pointsToLeftshift=0;


%FILLING IN THE FLAGS FOR THE FID.RAW DATA
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;

out.flags.addedrcvrs=addedrcvrsrffids;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.avgNormalized=0;
out.flags.isISIS=0;


if isRef
    %FILLING IN DATA STRUCTURE FOR THE FID.REF DATA
    ref.fids=reffids;
    ref.specs=refspecs;
    ref.sz=szreffids;
    ref.ppm=ppm;
    ref.t=t;
    ref.spectralwidth=spectralwidth;
    ref.dwelltime=dwelltime;
    ref.txfrq=txfrq;
    ref.date=date;
    ref.dims=refdims;
    ref.Bo=Bo;
    ref.averages=rawAverages;
    ref.rawAverages=rawAverages;
    ref.subspecs=subspecs;
    ref.rawSubspecs=rawSubspecs;
    ref.seq=sequence;
    ref.te=te;
    ref.tr=tr;
    ref.pointsToLeftshift=0;
    
    
    %FILLING IN THE FLAGS FOR THE FID.REF DATA
    ref.flags.writtentostruct=1;
    ref.flags.gotparams=1;
    ref.flags.leftshifted=0;
    ref.flags.filtered=0;
    ref.flags.zeropadded=0;
    ref.flags.freqcorrected=0;
    ref.flags.phasecorrected=0;
    
    ref.flags.addedrcvrs=addedrcvrsreffids;
    ref.flags.subtracted=0;
    ref.flags.writtentotext=0;
    ref.flags.downsampled=0;
    ref.flags.avgNormalized=0;
    ref.flags.isISIS=0;
    %     if ref.dims.subSpecs==0
    %         ref.flags.isISIS=0;
    %     else
    %         ref.flags.isISIS=(ref.sz(ref.dims.subSpecs)==4);
    %     end
end
