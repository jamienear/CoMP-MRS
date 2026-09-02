% sim_sLASER_shaped_pc.m
% Dana Goerzen (McGill University, 2019).
% Jamie Near (McGill University, 2020).
% Vincent Boer (DRCMR, 2020)
% 
% USAGE:
%  out = sim_sLASER_shaped_pc(n,sw,Bfield,linewidth,sys,te,RF,tp,dx,dy,Gx,Gy,pcSteps,flipAngle,centreFreq)
% 
% DESCRIPTION:
% This function simulates the semi-LASER experiment as described by Oz et al. (2018).
% The excitation is simulated as an instantaneous rotation, and the two pairs of slice
% selective refocusing adiabatic pulse are simulated as shaped RF pulses.

% This code differs from 'sim_sLASER_shaped.m' in that it enables the 
% choice of the number of phase encode steps per AFP pulse.  If a four step 
% phase cycling scheme chosen (pcSteps=4), then each AFP pulse will undergo 
% a 4 step cycle 0-90-180-270.  The density matrix from each phase cycle 
% step is added together and then divided by four.
% 
% Finally, this code simulates the spectrum at a given point in space (x,y),
% given the values of the slice selection gradients (Gx, and Gy).  The pulse
% waveform is assumed to be the same for both refocusing pulses.  In order
% to fully simulate the sLASER experiment, you have to run this
% simulation many times at various points in space (x,y), and then summate
% and scale together the resulting spectra.  
% 
% Feb 2020 - Jamie Near:  This code now accepts gradient modulated pulses.  
% If the input pulse is gradient modulated (waveform has 4 columns), then 
% the input parameters Gx and Gy are scaling factors for the GM function
% in order to achieve the desired slice thickness.
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% te        = echo time of sLASER experiment (ms)
% RF        = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
% tp        = RF pulse duration in [ms]
% dx        = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% dy        = position offset in y-direction (corresponding to second refocusing pulse) [cm]
% Gx        = If RF is not a gradient modulated pulse, Gx is the gradient strength
%             for first selective refocusing pulse [G/cm].  If RF is a gradient
%             modulated pulse, then Gx is a scaling factor for the GM
%             function to achieve the desired slice thickness in the
%             x-direction.  
% Gy        = If RF is not a gradient modulated pulse, Gy is the gradient strength
%             for first selective refocusing pulse [G/cm].  If RF is a gradient
%             modulated pulse, then Gy is a scaling factor for the GM
%             function to achieve the desired slice thickness in the
%             y-direction.  
% pcSteps   = Number of phase encoding steps per AFP pulse.  The four steps
%             will be evenly distributed about a full 2*pi cycle.
% flipAngle = flip angle of refocusing pulses [degrees] (Optional.  Default = 180 deg)
% centreFreq= centre frequency of the spectrum in [ppm] (Optional.  Default = 2.3)
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using PRESS 
%             sequence.

function out = sim_sLASER_shaped_pc(n,sw,Bfield,linewidth,sys,te,RF,tp,dx,dy,Gx,Gy,pcSteps,flipAngle,centreFreq)

if nargin<15
    centreFreq=2.3;
    if nargin<14
        flipAngle=180;
    end
end

%Check if this is a gradient modulated pulse.  If so, scale the GM functions
% according to Gx and Gy, and then set Gx and Gy both equal to zero:
if RF.isGM
    RFX=rf_scaleGrad(RF,Gx);
    RFY=rf_scaleGrad(RF,Gy);
    Gx=0;
    Gy=0;
else
    RFX=RF;
    RFY=RF;
end

 if (te/4)<(tp/1000)
     error('ERROR: the duration of the refocusing pulse cannot be longer than a quarter of the echo time! ABORTING!!');
 end

%initialize evolution times
 tau1=(te/4-tp)/2;
 tau2=te/4-tp;

%Set water to centre
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%calculte the phase cycling step values in Degrees:
phArray=linspace(0,360-(360/pcSteps),pcSteps);

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%BEGIN sLASER PULSE SEQUENCE************ 
d=sim_excite(d,H,'x');                                          %EXCITE instantaneously
d=sim_evolve(d,H,tau1/1000);                                    %Evolve by tau1

d_pc = cell(size(d));                                           %Create empty cell for phase cycling;
for k=1:length(phArray)                                         %Loop through all phase cycle steps for first AFP pulse.  
    d_temp=sim_shapedRF(d,H,RFX,tp,flipAngle,phArray(k),dx,Gx); %1st shaped 180 degree adiabatic refocusing pulse along X gradient
    d_temp=sim_dMul(d_temp,(-1).^(phArray(k)/90));              %Vincent Boer:  Multiply by -1 for 90 and 270 degree
    d_pc=sim_dAdd(d_pc,d_temp);                                 %Add the density matrices from each phase cycle step
end
d=sim_dDiv(d_pc,length(phArray));                               %Divide by the number of steps

d=sim_evolve(d,H,tau2/1000);                                    %Evolve by tau2

d_pc = cell(size(d));                                           %Create empty cell for phase cycling;
for k=1:length(phArray)                                                 %Loop through all phase cycle steps for second AFP pulse. 
    d_temp=sim_shapedRF(d,H,RFX,tp,flipAngle,phArray(k),dx,Gx); %2nd shaped 180 degree adiabatic refocusing pulse along X gradient
    d_temp=sim_dMul(d_temp,(-1).^(phArray(k)/90));              %Vincent Boer:  Multiply by -1 for 90 and 270 degree
    d_pc=sim_dAdd(d_pc,d_temp);                                 %Add the density matrices from each phase cycle step
end
d=sim_dDiv(d_pc,length(phArray));                               %Divide by the number of steps

d=sim_evolve(d,H,tau2/1000);                                    %Evolve by tau2

d_pc = cell(size(d));                                           %Create empty cell for phase cycling;
for k=1:length(phArray)                                         %Loop through all phase cycle steps for third AFP pulse. 
    d_temp=sim_shapedRF(d,H,RFY,tp,flipAngle,phArray(k),dy,Gy); %3rd shaped 180 degree adiabatic refocusing pulse along Y gradient
    d_temp=sim_dMul(d_temp,(-1).^(phArray(k)/90));              %Vincent Boer:  Multiply by -1 for 90 and 270 degree
    d_pc=sim_dAdd(d_pc,d_temp);                                 %Add the density matrices from each phase cycle step
end
d=sim_dDiv(d_pc,length(phArray));                               %Divide by the number of steps

d=sim_evolve(d,H,tau2/1000);                                    %Evolve by tau2

d_pc = cell(size(d));                                           %Create empty cell for phase cycling;
for k=1:length(phArray)                                         %Loop through all phase cycle steps for fourth AFP pulse. 
    d_temp=sim_shapedRF(d,H,RFY,tp,flipAngle,phArray(k),dy,Gy); %4th shaped 180 degree adiabatic refocusing pulse along Y gradient
    d_temp=sim_dMul(d_temp,(-1).^(phArray(k)/90));              %Vincent Boer:  Multiply by -1 for 90 and 270 degree
    d_pc=sim_dAdd(d_pc,d_temp);                                 %Add the density matrices from each phase cycle step
end
d=sim_dDiv(d_pc,length(phArray));                               %Divide by the number of steps

d=sim_evolve(d,H,tau1/1000);                                    %Evolve by tau1

[out,dout]=sim_readout(d,H,n,sw,linewidth,90);                  %Readout along +y axis (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='semi-LASER';
out.te=te;
out.sim='shaped';

%Additional fields for compatibility with FID-A processing tools.
out.sz=size(out.specs);
out.date=date;
out.dims.t=1;
out.dims.coils=0;
out.dims.averages=0;
out.dims.subSpecs=0;
out.dims.extras=0;
out.averages=1;
out.rawAverages=1;
out.subspecs=1;
out.rawSubspecs=1;
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
out.flags.addedrcvrs=1;
out.flags.subtracted=1;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.isFourSteps=0;


end
