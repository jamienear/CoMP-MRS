function [VoxPul1parameters, VoxPul2parameters, VoxPul3parameters, EdPulparameters, VoxPul1w1max, VoxPul2w1max,VoxPul3w1max, EdPulw1max] = RF_pulse_profile(inDir,homedir)
close all
% RUN IN HOMEDIR
addpath(genpath(homedir))
% Give dataset folder with specific exam folder to get RF pulses params

% Get file name
filename= extractAfter(inDir,'Data/');
filename= strrep(filename,'_','-');
 
mkdir(filename)
cd(filename)

%Now get sequence name
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



if sequence == 'MEGAPRESS/HERMES'
                if sum(strfind(Ed1Freqppm,'-1'))==0
                    sequence='HERMES';
                elseif sum(strfind(Ed2Freqppm,'-1'))==0
                    sequence='HERMES';
                elseif sum(strfind(Ed3Freqppm,'-1'))==0
                    sequence='HERMES';
                elseif sum(strfind(Ed4Freqppm,'-1'))==0
                    sequence='HERMES';
                else 
                    sequence='MEGAPRESS';
                end
end

if ~strcmp(Ed1Freqppm,Ed2Freqppm) && ~strcmp(Ed1Freqppm,Ed3Freqppm) && ~strcmp(Ed1Freqppm,Ed4Freqppm)
    sequence = 'HERMES';
else
    sequence = 'MEGAPRESS';
end

else
    sequence = 'PRESS';
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VoxPul1 - EXCITATION PULSE

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$VoxPul1Shape=( 600 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$VoxPul1Shape=( 600 )');
end

% Save the start line for the search to be started
start_line = line(index:end-1);

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$VoxPul2Shape=( 600 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$VoxPul2Shape=( 600 )');
end

% Save the end line for the search already done
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read the record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the data
        VoxPul1Shape=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);   % get line 
            if VoxPul1Shape==0  % if nothing added yet to the variable chosen, then initiate the assignment of the variable to the info detected on the first line
                VoxPul1Shape=l;
            else
                VoxPul1Shape= horzcat(VoxPul1Shape,l);  % else, if there is info stored by the chosen variable, then add the new info horizontally (keep the array format)
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
VoxPul1Shape = (VoxPul1Shape(1:end-length(end_line)));
% Convert variable to a number
VoxPul1Shape = str2num(VoxPul1Shape);

fclose(method_fid); % close the file

% Separate phase and magnitude 
VoxPul1=zeros(size(VoxPul1Shape,2)/2,2);
VoxPul1(:,2)= VoxPul1Shape(1:2:end);
VoxPul1(:,1)= VoxPul1Shape(2:2:end);

% Scale the magnitude component of the RF pulses by the maximum value
% detected
VoxPul1(:,2) = VoxPul1(:,2)/100;

figure,
% plot(VoxPul1(:,1),'Linewidth',2)
hold on
plot(VoxPul1(:,2),'Linewidth',2)
%legend('ph', 'amp')
title([filename '-' sequence '-VoxPul1Shape'])
saveas(gcf, ['VoxPul1-' sequence], 'jpg')
saveas(gcf, ['VoxPul1-' sequence], 'fig')
saveas(gcf, [homedir 'VoxPul1-' sequence], 'jpg')
%%% End VoxPul1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VoxPul2 - REFOCUSING PULSE
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$VoxPul2Shape=( 600 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$VoxPul2Shape=( 600 )');
end

% Save the start line 
start_line = line(index:end-1);

switch sequence
    case 'PRESS'
        method_fid=fopen([inDir '/method']);
        line=fgets(method_fid);
        index=findstr(line,'##$VoxPul3Shape=( 600 )');
        while isempty(index)
            line=fgets(method_fid);
            index=findstr(line,'##$VoxPul3Shape=( 600 )');
        end
    otherwise
        method_fid=fopen([inDir '/method']);
        line=fgets(method_fid);
        index=findstr(line,'$$ @vis= Ed4FreqHz MEGAWS MEGA EdEnum EdPul EdAmpl EdShape VoxPul1Shape');
        while isempty(index)
            line=fgets(method_fid);
            index=findstr(line,'$$ @vis= Ed4FreqHz MEGAWS MEGA EdEnum EdPul EdAmpl EdShape VoxPul1Shape');
        end
end
% Save the end line 
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read a record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
        VoxPul2Shape=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);
            if VoxPul2Shape==0
                VoxPul2Shape=l;
            else
                VoxPul2Shape= horzcat(VoxPul2Shape,l);
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
VoxPul2Shape = (VoxPul2Shape(1:end-length(end_line)));

% Convert variable to a number
VoxPul2Shape = str2num(VoxPul2Shape);

fclose(method_fid);

VoxPul2=zeros(size(VoxPul2Shape,2)/2,2);
VoxPul2(:,2)= VoxPul2Shape(1:2:end);
VoxPul2(:,1)= VoxPul2Shape(2:2:end);

% Scale the magnitude component of the RF pulses by the maximum value
% detected
VoxPul2(:,2) = VoxPul2(:,2)/100;

figure,
% plot(VoxPul2(:,1),'Linewidth',2)
hold on
plot(VoxPul2(:,2),'Linewidth',2)
%legend('ph', 'amp')
title([filename '-' sequence '-VoxPul2Shape'])
saveas(gcf, ['VoxPul2-' sequence], 'jpg')
saveas(gcf, ['VoxPul2-' sequence], 'fig')
saveas(gcf, [homedir 'VoxPul2-' sequence], 'jpg')

%%% End VoxPul2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VoxPul3 - REFOCUSING PULSE
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$VoxPul3Shape=( 600 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$VoxPul3Shape=( 600 )');
end

% Save the start line 
start_line = line(index:end-1);


method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$LockVoxPulAmpls=No');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$LockVoxPulAmpls=No');
end

% Save the end line 
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read a record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
        VoxPul3Shape=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);
            if VoxPul3Shape==0
                VoxPul3Shape=l;
            else
                VoxPul3Shape= horzcat(VoxPul3Shape,l);
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
VoxPul3Shape = (VoxPul3Shape(1:end-length(end_line)));
% Convert variable to a number
VoxPul3Shape = str2num(VoxPul3Shape);

fclose(method_fid);

VoxPul3=zeros(size(VoxPul3Shape,2)/2,2);
VoxPul3(:,2)= VoxPul3Shape(1:2:end);
VoxPul3(:,1)= VoxPul3Shape(2:2:end);

% Scale the magnitude component of the RF pulses by the maximum value
% detected
VoxPul3(:,2) = VoxPul3(:,2)/100;

figure,
% plot(VoxPul3(:,1),'Linewidth',2)
hold on
plot(VoxPul3(:,2),'Linewidth',2)
%legend('ph', 'amp')
title([filename '-' sequence '-VoxPul3Shape'])
saveas(gcf, ['VoxPul3-' sequence], 'jpg')
saveas(gcf, ['VoxPul3-' sequence], 'fig')
saveas(gcf, [homedir 'VoxPul3-' sequence], 'jpg')

%%% End VoxPul3

if strcmp(sequence, 'HERMES') || strcmp(sequence, 'MEGAPRESS')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EdShape - EDITING PULSE
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$EdShape=( 400 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$EdShape=( 400 )');
end

% Save the start line 
start_line = line(index:end-1);

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$VoxPul1Shape=( 600 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$VoxPul1Shape=( 600 )');
end

% Save the end line 
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read a record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
        EdShape=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);
            if EdShape==0
                EdShape=l;
            else
                EdShape= horzcat(EdShape,l);
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
EdShape = (EdShape(1:end-length(end_line)));

% Convert variable to a number
EdShape = str2num(EdShape);

fclose(method_fid);

Ed=zeros(size(EdShape,2)/2,2);
Ed(:,2)= EdShape(1:2:end);
Ed(:,1)= EdShape(2:2:end);

% Scale the magnitude component of the RF pulses by the maximum value
% detected
Ed(:,2) = Ed(:,2)/100;

figure,
% plot(Ed(:,1),'Linewidth',2)
hold on
plot(Ed(:,2),'Linewidth',2)
%legend('ph', 'amp')
title([filename '-' sequence '-EdShape'])
saveas(gcf, ['EdShape-' sequence], 'jpg')
saveas(gcf, ['EdShape-' sequence], 'fig')
saveas(gcf, [homedir 'EdShape-' sequence], 'jpg')

%%% End EdShape

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ed1Shape
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$Ed1Shape=( 400 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$Ed1Shape=( 400 )');
end

% Save the start line 
start_line = line(index:end-1);

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$Ed2Shape=( 400 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$Ed2Shape=( 400 )');
end

% Save the end line 
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read a record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
        Ed1Shape=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);
            if Ed1Shape==0
                Ed1Shape=l;
            else
                Ed1Shape= horzcat(Ed1Shape,l);
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
Ed1Shape = (Ed1Shape(1:end-length(end_line)));

% Convert variable to a number
Ed1Shape = str2num(Ed1Shape);

fclose(method_fid);

Ed1=zeros(size(Ed1Shape,2)/2,2);
Ed1(:,2)= Ed1Shape(1:2:end);
Ed1(:,1)= Ed1Shape(2:2:end);

% Scale the magnitude component of the RF pulses by the maximum value
% detected
Ed1(:,2) = Ed1(:,2)/100;

figure,
%plot(Ed1(:,1),'Linewidth',2)
hold on
plot(Ed1(:,2),'Linewidth',2)
%legend('ph', 'amp')
title([filename '-' sequence '-Ed1Shape'])
saveas(gcf, ['Ed1Shape' sequence], 'jpg')
saveas(gcf, ['Ed1Shape' sequence], 'fig')
saveas(gcf, [homedir 'Ed1Shape-' sequence], 'jpg')

%%% End Ed1Shape

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ed2Shape
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$Ed2Shape=( 400 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$Ed2Shape=( 400 )');
end

% Save the start line 
start_line = line(index:end-1);

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$Ed3Shape=( 400 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$Ed3Shape=( 400 )');
end

% Save the end line 
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read a record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
        Ed2Shape=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);
            if Ed2Shape==0
                Ed2Shape=l;
            else
                Ed2Shape= horzcat(Ed2Shape,l);
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
Ed2Shape = (Ed2Shape(1:end-length(end_line)));

% Convert variable to a number
Ed2Shape = str2num(Ed2Shape);

fclose(method_fid);

Ed2=zeros(size(Ed2Shape,2)/2,2);
Ed2(:,2)= Ed2Shape(1:2:end);
Ed2(:,1)= Ed2Shape(2:2:end);

% Scale the magnitude component of the RF pulses by the maximum value
% detected
Ed2(:,2) = Ed2(:,2)/100;

figure,
%plot(Ed2(:,1),'Linewidth',2)
hold on
plot(Ed2(:,2),'Linewidth',2)
%legend('ph', 'amp')
title([filename '-' sequence '-Ed2Shape'])
saveas(gcf, ['Ed2Shape' sequence], 'jpg')
saveas(gcf, ['Ed2Shape' sequence], 'fig')
saveas(gcf, [homedir 'Ed2Shape-' sequence], 'jpg')

%%% End Ed2Shape

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ed3Shape
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$Ed3Shape=( 400 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$Ed3Shape=( 400 )');
end

% Save the start line 
start_line = line(index:end-1);

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$Ed4Shape=( 400 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$Ed4Shape=( 400 )');
end

% Save the end line 
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read a record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
        Ed3Shape=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);
            if Ed3Shape==0
                Ed3Shape=l;
            else
                Ed3Shape= horzcat(Ed3Shape,l);
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
Ed3Shape = (Ed3Shape(1:end-length(end_line)));

% Convert variable to a number
Ed3Shape = str2num(Ed3Shape);

fclose(method_fid);

Ed3=zeros(size(Ed3Shape,2)/2,2);
Ed3(:,2)= Ed3Shape(1:2:end);
Ed3(:,1)= Ed3Shape(2:2:end);

% Scale the magnitude component of the RF pulses by the maximum value
% detected
Ed3(:,2) = Ed3(:,2)/100;

figure,
%plot(Ed3(:,1),'Linewidth',2)
hold on
plot(Ed3(:,2),'Linewidth',2)
%legend('ph', 'amp')
title([filename '-' sequence '-Ed3Shape'])
saveas(gcf, ['Ed3Shape' sequence], 'jpg')
saveas(gcf, ['Ed3Shape' sequence], 'fig')
saveas(gcf, [homedir 'Ed3Shape-' sequence], 'jpg')

% End Ed3Shape

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ed4Shape
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$Ed4Shape=( 400 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$Ed4Shape=( 400 )');
end

% Save the start line 
start_line = line(index:end-1);

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$Ed1Freqppm=( 2 )');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$Ed1Freqppm=( 2 )');
end

% Save the end line 
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read a record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
        Ed4Shape=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);
            if Ed4Shape==0
                Ed4Shape=l;
            else
                Ed4Shape= horzcat(Ed4Shape,l);
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
Ed4Shape = (Ed4Shape(1:end-length(end_line)));

% Convert variable to a number
Ed4Shape = str2num(Ed4Shape);

fclose(method_fid);

Ed4=zeros(size(Ed4Shape,2)/2,2);
Ed4(:,2)= Ed4Shape(1:2:end);
Ed4(:,1)= Ed4Shape(2:2:end);

% Scale the magnitude component of the RF pulses by the maximum value
% detected
Ed4(:,2) = Ed4(:,2)/100;

figure,
%plot(Ed4(:,1),'Linewidth',2)
hold on
plot(Ed4(:,2),'Linewidth',2)
%legend('ph', 'amp')
title([filename '-' sequence '-Ed4Shape'])
saveas(gcf, ['Ed4Shape' sequence], 'jpg')
saveas(gcf, ['Ed4Shape' sequence], 'fig')
saveas(gcf, [homedir 'Ed4Shape-' sequence], 'jpg')

%%% End Ed4Shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
end

% Determine localisation and editing pulses durations, bandwidths, FA,
% time-bandwidts products and powers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VoxPul1

% Determine the duration
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$VoxPul3Enum=<Calculated>');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$VoxPul3Enum=<Calculated>');
end

% Save the start line 
start_line = line(index:end-1);

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$VoxPul2=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$VoxPul2=');
end

% Save the end line 
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read a record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
        VoxPul1params=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);
            if VoxPul1params==0
                VoxPul1params=l;
            else
                VoxPul1params= horzcat(VoxPul1params,l);
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
VoxPul1params = (VoxPul1params(1:end-length(end_line)));
comasinthestring = strfind(VoxPul1params,',');
fclose(method_fid);
% Determine pulse duration (ms), bandwidth (Hz), FA, time-bandwidth product (Hz*ms), power  
VoxPul1parameters.Tp = str2num(VoxPul1params(strfind(VoxPul1params,'(')+1:comasinthestring(1)-1));
VoxPul1parameters.bw = str2num(VoxPul1params(comasinthestring(1)+2:comasinthestring(2)-1));
VoxPul1parameters.FA = str2num(VoxPul1params(comasinthestring(2)+2:comasinthestring(3)-1));
VoxPul1parameters.bwfac = str2num(VoxPul1params(comasinthestring(5)+2:comasinthestring(6)-1));
VoxPul1parameters.pow = str2num(VoxPul1params(comasinthestring(size(comasinthestring,2)-1)+2:comasinthestring(size(comasinthestring,2))-1));

%%% End VoxPul1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VoxPul2

% Determine the duration
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'<$VoxPul1Shape>)');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'<$VoxPul1Shape>)');
end

% Save the start line 
start_line = line(index:end-1);

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$VoxPul3=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$VoxPul3=');
end

% Save the end line 
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read a record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
        VoxPul2params=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);
            if VoxPul2params==0
                VoxPul2params=l;
            else
                VoxPul2params= horzcat(VoxPul2params,l);
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
VoxPul2params = (VoxPul2params(1:end-length(end_line)));
comasinthestring = strfind(VoxPul2params,',');
fclose(method_fid);
% Determine pulse duration (ms), bandwidth (Hz), FA, time-bandwidth product (Hz*ms), power  
VoxPul2parameters.Tp = str2num(VoxPul2params(strfind(VoxPul2params,'(')+1:comasinthestring(1)-1));
VoxPul2parameters.bw = str2num(VoxPul2params(comasinthestring(1)+2:comasinthestring(2)-1));
VoxPul2parameters.FA = str2num(VoxPul2params(comasinthestring(2)+2:comasinthestring(3)-1));
VoxPul2parameters.bwfac = str2num(VoxPul2params(comasinthestring(5)+2:comasinthestring(6)-1));
VoxPul2parameters.pow = str2num(VoxPul2params(comasinthestring(size(comasinthestring,2)-1)+2:comasinthestring(size(comasinthestring,2))-1));

%%% End VoxPul2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VoxPul2

% Determine the duration
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'<$VoxPul2Shape>)');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'<$VoxPul2Shape>)');
end

% Save the start line 
start_line = line(index:end-1);

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$VoxPul1Ampl');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$VoxPul1Ampl');
end

% Save the end line 
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read a record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
        VoxPul3params=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);
            if VoxPul3params==0
                VoxPul3params=l;
            else
                VoxPul3params= horzcat(VoxPul3params,l);
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
VoxPul3params = (VoxPul3params(1:end-length(end_line)));
comasinthestring = strfind(VoxPul3params,',');
fclose(method_fid);
% Determine pulse duration (ms), bandwidth (Hz), FA, time-bandwidth product (Hz*ms), power  
VoxPul3parameters.Tp = str2num(VoxPul3params(strfind(VoxPul3params,'(')+1:comasinthestring(1)-1));
VoxPul3parameters.bw = str2num(VoxPul3params(comasinthestring(1)+2:comasinthestring(2)-1));
VoxPul3parameters.FA = str2num(VoxPul3params(comasinthestring(2)+2:comasinthestring(3)-1));
VoxPul3parameters.bwfac = str2num(VoxPul3params(comasinthestring(5)+2:comasinthestring(6)-1));
VoxPul3parameters.pow = str2num(VoxPul3params(comasinthestring(size(comasinthestring,2)-1)+2:comasinthestring(size(comasinthestring,2))-1));

%%% End VoxPul3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EdPul

EdPulparameters=0;
EdPulw1max=0;
if strcmp(sequence, 'HERMES') || strcmp(sequence, 'MEGAPRESS')

% Determine the duration
method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$EdEnum=<Calculated>');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$EdEnum=<Calculated>');
end

% Save the start line 
start_line = line(index:end-1);

method_fid=fopen([inDir '/method']);
line=fgets(method_fid);
index=findstr(line,'##$EdAmpl=');
while isempty(index)
    line=fgets(method_fid);
    index=findstr(line,'##$EdAmpl=');
end

% Save the end line 
end_line = line(index:end-1);

% Extract data between start_line and end_line
fid=fopen([inDir '/method']);  % open the file
while ~feof(fid)  % loop to the end of file
    l=fgetl(fid);            % read a record of each line
    if strfind(l,start_line)  % if  l, the last line found is equal to start_line, start saving the water data
        EdPulparams=0;
        while (strcmp(l,end_line)==0)  % ... until it reaches the end_line
            l=fgetl(fid);
            if EdPulparams==0
                EdPulparams=l;
            else
                EdPulparams= horzcat(EdPulparams,l);
            end
        end
        
    end
end

% Remove the last characters corresponding to the end_line
EdPulparams = (EdPulparams(1:end-length(end_line)));
comasinthestring = strfind(EdPulparams,',');
fclose(method_fid);
% Determine pulse duration (ms), bandwidth (Hz), FA, time-bandwidth product (Hz*ms), power  
EdPulparameters.Tp = str2num(EdPulparams(strfind(EdPulparams,'(')+1:comasinthestring(1)-1));
EdPulparameters.bw = str2num(EdPulparams(comasinthestring(1)+2:comasinthestring(2)-1));
EdPulparameters.FA = str2num(EdPulparams(comasinthestring(2)+2:comasinthestring(3)-1));
EdPulparameters.bwfac = str2num(EdPulparams(comasinthestring(5)+2:comasinthestring(6)-1));
EdPulparameters.pow = str2num(EdPulparams(comasinthestring(size(comasinthestring,2)-1)+2:comasinthestring(size(comasinthestring,2))-1));
end
%%% End EdPul
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the RF pulses strucures and get peak B1

%  % [RF_struct,w1max]=rf2RFwaveform(rf,type,f0,Tp,w1max,bw);
%   
%   DESCRIPTION:
%   Initialize an RF pulse structure to contain an 
%   RF Pulse waveform as well as its accompanying header
%   information.  This function finds the time-bandwidth
%   product (tbw) and the time-w1 product (tw1) of the pulse,
%   and stores this information in the header fields of the 
%   output RF structure.
%   
%   INPUTS:
%   rf        = rf pulse array. with phase in 1st colume, magnitude in
%               2nd column and time step in 3rd column
%   type      = Excitation ('exc'), Refocusing ('ref') or Inversion ('inv')
%   f0        = centre frequency of the rf pulse [Hz].  Optional. Default=0.
%   Tp        = durarion of pulse [s].  Optional. Default=0.00ß5.
%   w1max     = Default B1 max (typically assuming a 5 ms Tp)
%   bw        = bandwidth

[VoxPul1_struct VoxPul1w1max]= rf2RFwaveform_Bruker(VoxPul1,'exc',0,VoxPul1parameters.Tp/1000,0);
[VoxPul2_struct VoxPul2w1max]= rf2RFwaveform_Bruker(VoxPul2,'ref',0,VoxPul2parameters.Tp/1000,0);
[VoxPul3_struct VoxPul3w1max]= rf2RFwaveform_Bruker(VoxPul3,'ref',0,VoxPul3parameters.Tp/1000,0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exc_RFpulse.waveform = VoxPul1;
rfc_RFpulse.waveform = VoxPul2;

% Use round to obtain phase values of 180 and 0
exc_RFpulse.waveform(:,1) = round(exc_RFpulse.waveform(:,1));
rfc_RFpulse.waveform(:,1) = round(rfc_RFpulse.waveform(:,1));

% % Switch column 1 with column 2 to have phase and magnitude
% exc_RFpulse.waveform(:,[1 2]) = exc_RFpulse.waveform(:,[2 1]);
% rfc_RFpulse.waveform(:,[1 2]) = rfc_RFpulse.waveform(:,[2 1]);

% The numeric values of the phase and magnitude for the exc and rfc
% pulses are a bit weird so the follwing changes must be done:
exc_RFpulse.waveform(1,1) = 180;
exc_RFpulse.waveform(1,2) = 0;
exc_RFpulse.waveform(size(exc_RFpulse.waveform,1),2) = 0;
rfc_RFpulse.waveform(size(rfc_RFpulse.waveform,1),2) = 0;

% Add the third column filled with 1 only (giving the timestep)
exc_RFpulse.waveform(:,3) = ones(length(exc_RFpulse.waveform(:,1)),1);
rfc_RFpulse.waveform(:,3) = ones(length(rfc_RFpulse.waveform(:,1)),1);

% Convert the RF structure into a .pta file
io_writepta(exc_RFpulse,'PRESS_excRF.pta');
io_writepta(rfc_RFpulse,'PRESS_rfcRF.pta');

% To check the code, load the .pta file and plot the magnitude values simulated
% and the RF waveform; they should be identical
excwaveform = 'PRESS_excRF.pta';
excRF = io_loadRFwaveform(excwaveform,'inv',0);

rfcwaveform = 'PRESS_rfcRF.pta';
rfcRF = io_loadRFwaveform(rfcwaveform,'inv',0);

figure, plot(excRF.waveform(:,2))
title('PRESS excitation RF pulse')

saveas(gca,['PRESS_exc_pulse'],'fig');
saveas(gca,['PRESS_exc_pulse'],'jpg');

figure, plot(rfcRF.waveform(:,2))
title('PRESS refocusing RF pulse')
saveas(gca,['PRESS_rfc_pulse'],'fig');
saveas(gca,['PRESS_rfc_pulse'],'jpg');

excmag = excRF.waveform(:,2);
excphase = excRF.waveform(:,1);
excmag(excphase>0) = -excmag(excphase>0);
figure, plot(excmag)
title('PRESS excitation RF pulse (real)')
saveas(gca,['PRESS_exc_pulse(real)'],'fig');
saveas(gca,['PRESS_exc_pulse(real)'],'jpg');

rfcmag = rfcRF.waveform(:,2);
rfcphase = rfcRF.waveform(:,1);
excmag(rfcphase>0) = -rfcmag(rfcphase>0);
figure, plot(rfcmag)
title('PRESS refocusing RF pulse (real)')
saveas(gca,['PRESS_rfc_pulse(real)'],'fig');
saveas(gca,['PRESS_rfc_pulse(real)'],'jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

io_writepta_Bruker(VoxPul1_struct,[homedir 'PRESS_excRF.pta']);
io_writepta_Bruker(VoxPul2_struct,[homedir 'PRESS_rfcRF.pta']);
io_writepta_Bruker(VoxPul3_struct,[homedir 'PRESS_rfcRF.pta']);

if strcmp(sequence, 'HERMES') || strcmp(sequence, 'MEGAPRESS')    
    
[EdPul_struct EdPulw1max]= rf2RFwaveform_Bruker(Ed,'inv',0,EdPulparameters.Tp/1000,0);

% Ed1_struct = rf2RFwaveform(Ed1,'inv',0);
% Ed2_struct = rf2RFwaveform(Ed2,'inv',0);
% Ed3_struct = rf2RFwaveform(Ed3,'inv',0);
% Ed4_struct = rf2RFwaveform(Ed4,'inv',0);

% [mv,sc]=rf_blochSim(RF,tp,fspan,f0,peakB1,ph,M0);
%   
%   DESCRIPTION:
%   Perform a bloch simulation of an RF pulse.  This code simply runs Martyn 
%   Klassen's excellent bloch equation simulator.  For more information, see
%   help file for bes.m.  (~FID-A/rfPulseTools/mklassenTools/bes.m).
%   
%   INPUTS:
%   RF        = RF pulse definition structure
%   tp        = pulse duration in [ms]
%   fspan     = Frequency span in [kHz] (optional.  Default=10kHz)
%   f0        = Centre of frequnecy span [kHz] (optional.  Default=0)
%   peakB1 	= Peak B1 amplitude in [kHz] (optional.  Default=RF.tw1/tp)
%   ph        = Starting phase of the rf pulse [degrees] (optional. Default=0)
%   M0        = Starting magnetization [units of M0] (optional. Default=[0,0,1])
%  
%   OUTPUTS:
%   mv        = Simulated magnetization vector in three columns (x,y,z) as a
%               function of frequency.
%   sc        = Frequency scale (in kHz) corresponding to the simulated mv vectors.
 
% Calculate peak B1 for the editing pulse (single band)
w1max_range = [0.05:0.001:0.2];
for j=1:numel(w1max_range)   
[nv,sc] = rf_blochSim(EdPul_struct, EdPulparameters.Tp,2,0,w1max_range(j));
Mz1(j)=nv(3,size(nv,2)/2);
close all
end

w1max_range = [0.137:0.001:0.141];
for j=1:numel(w1max_range)   
[nv,sc] = rf_blochSim(EdPul_struct, EdPulparameters.Tp,2,0,w1max_range(j));
Mz2(j)=nv(3,size(nv,2)/2);
close all
end

figure,
plot(w1max_range, Mz2)
saveas(gcf, 'SimPulse_w1max_Mz', 'jpg')
ind=find(Mz2==min(Mz2), 1, 'first');

[nv,sc] = rf_blochSim(EdPul_struct, EdPulparameters.Tp,2,0,w1max_range(ind));
s = ((1000*sc)/400)+1.9;
figure,
plot(s,nv(3,:))
ylim([-inf 2])
xlim([-inf inf])
saveas(gcf, 'EdPul_ReImg_GABA', 'jpg')
saveas(gcf, [homedir 'EdPul_ReImg_GABA' sequence], 'jpg')

% Calculate the composite pulse based on the code lines in backbone.o
ind_amp=find(EdPul_struct.waveform(:,2)==max(EdPul_struct.waveform(:,2)), 1, 'first');
EdPulTp=EdPulparameters.Tp/1000;
t=-EdPulTp/2:EdPulTp/200:+EdPulTp/2-EdPulTp/200;
cf=(1.9+4.56)/2;
offset1=(1.9-cf);
offset2=(4.56-cf);
offset1Hz=400.335061568014*offset1;
offset2Hz=400.335061568014*offset2;

% ph1=cos(2.0*pi*offset[0]*t);
% ph2=sin(-2.0*pi*offset[0]*t); //fishy minus sign?
% re =ph1*wavein[2*i]*cos(wavein[2*i+1]*rad)-ph2*wavein[2*i]*sin(wavein[2*i+1]*rad);
% im =ph1*wavein[2*i]*sin(wavein[2*i+1]*rad)+ph2*wavein[2*i]*cos(wavein[2*i+1]*rad);

ph1=cos(2.0*pi*offset1Hz*t);
ph2=sin(-2.0*pi*offset1Hz*t);

ph1=ph1';
ph2=ph2';

% Take magnitude and phase and convert them to real and imaginary
rad=pi/180.0;
re=EdPul_struct.waveform(:,2).*cos(EdPul_struct.waveform(:,1)*rad);
im=EdPul_struct.waveform(:,2).*sin(EdPul_struct.waveform(:,1)*rad);

re_initial= ph1.*re-ph2.*im;
im_initial= ph1.*im+ph2.*re;

figure,
plot(t,re_initial)
hold on
plot(t,im_initial)
saveas(gcf, [homedir 'EdPulComp1_ReImg' sequence], 'jpg')
% ph1=cos(2.0*pi*offset[1]*t);
% ph2=sin(-2.0*pi*offset[1]*t);
% re +=ph1*wavein[2*i]*cos(wavein[2*i+1]*rad)-ph2*wavein[2*i]*sin(wavein[2*i+1]*rad);
% im +=ph1*wavein[2*i]*sin(wavein[2*i+1]*rad)+ph2*wavein[2*i]*cos(wavein[2*i+1]*rad);

ph1=cos(2.0*pi*offset2Hz*t);
ph2=sin(-2.0*pi*offset2Hz*t);

ph1=ph1';
ph2=ph2';

% Take magnitude and phase and convert them to real and imaginary
rad=pi/180.0;
re=EdPul_struct.waveform(:,2).*cos(EdPul_struct.waveform(:,1)*rad);
im=EdPul_struct.waveform(:,2).*sin(EdPul_struct.waveform(:,1)*rad);

re_final= re_initial+ph1.*re-ph2.*im;
im_final= im_initial+ph1.*im+ph2.*re;

figure,
plot(t,re_final)
hold on
plot(t,im_final)
saveas(gcf, [homedir 'EdPulComp2_ReImg' sequence], 'jpg')

% Combine real and imaginary points and get the phase and the magnitude
tmp=re_initial+i*im_initial;
figure
plot(t,abs(tmp));

tmp=re_final+i*im_final;
figure
plot(t,abs(tmp));

magn=abs(tmp)/2;
ph=angle(tmp)/rad;

ph(ph<150)=-ph(ph<150);

% % Only to check RF pulse waveform:
% magn(ph>150)=-magn(ph>150);
% figure,plot(magn)

compRF=zeros(size(magn,1),2);
compRF(:,1)=ph;
compRF(:,2)=magn;

[EdPulComp_struct EdPulCompw1max] = rf2RFwaveform_Bruker(compRF, 'inv',0,EdPulparameters.Tp/1000,2*w1max_range(ind));
[nv,sc] = rf_blochSim(EdPulComp_struct, EdPulparameters.Tp,4,0,2*w1max_range(ind));
saveas(gcf, 'EdPulComp_w1max_Mz', 'jpg')
saveas(gcf, [homedir 'EdPulComp_w1max_Mz_' sequence], 'jpg')

s = ((1000*sc)/400)+(1.9+4.56)/2;
figure,
plot(s,nv(3,:))
ylim([-inf 2])
xlim([-inf inf])
saveas(gcf, 'EdPulComp_Re_Im', 'jpg')
saveas(gcf, [homedir 'EdPulComp_Re_Im_' sequence], 'jpg')

% For plotting purposes:
% For the composite pulse, do not use the structure, but instead 
% use the re_final variable where ph=+/-180(178) => magn=-magn
% and ph=+/-360 => magn=: magn
end
figure,
subplot(2,3,1)
plot(VoxPul1_struct.waveform(:,2),'Linewidth',2)
legend('VoxPul1')
subplot(2,3,2)
plot(VoxPul2_struct.waveform(:,2),'Linewidth',2)
legend('VoxPul2') 
subplot(2,3,3)
plot(VoxPul3_struct.waveform(:,2),'Linewidth',2)
legend('VoxPul3')

if strcmp(sequence, 'HERMES') || strcmp(sequence, 'MEGAPRESS')

subplot(2,3,4)
plot(EdPul_struct.waveform(:,2),'Linewidth',2)
legend('EdPul')
subplot(2,3,5)
plot(re_final,'Linewidth',2)
legend('EdPulComp')
end
saveas(gcf, 'allpulses', 'jpg')
saveas(gcf, [homedir 'allpulses_' sequence], 'jpg')

if strcmp(sequence, 'HERMES') || strcmp(sequence, 'MEGAPRESS')

io_writepta_Bruker(EdPul_struct,[homedir 'EdPul_singleband' sequence '.pta']);
io_writepta_Bruker(EdPulComp_struct,[homedir 'EdPul_dualband' sequence '.pta']);
end
end


