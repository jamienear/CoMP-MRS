% compMRS_simsLASER.m
% Dana Goerzen and Jamie Near, McGill University 2021.
% Fast version by Muhammad G Saleh (Johns Hopkins University School of Medicine, 2019)
% Revised for CoMP-MRS - Jamie Near and Diana Rotaru, 2025/2026.
%
% USAGE:
% out=compMRS_simsLASER(spinSys,simPars);
% 
% DESCRIPTION:
% This script simulates a PRESS experiment with fully shaped refocusing 
% pulses. Coherence order filtering is employed to only simulate desired signals
% This results in a 4x speed up compared to phase cycling (see deprecated run_simSemiLASERShaped_fast_phCyc.m)
% Furthermore, simulations are run at various locations in space to account for the 
% within-voxel spatial variation of the metabolite signal.  Summation 
% across spatial positions is performed. The MATLAB parallel computing toolbox 
% (parfor loop) was used to accelerate the simulations.  Acceleration 
% is currently performed in the direction of the slice selective pulse along
% the x-direction, but this can be changed.  Up to a factor of 12 acceleration
% can be achieved using this approach. To achieve 
% faster perfomance compared to the original 'run_simSemiLASER_shaped.m' function,
% this code uses the method described by Yan Zhang et al. Med Phys 2017;44(8): 
% 4169-78.  Some additional acceleration is currently performed using parfor 
% loops in both x and y directions.  To enable the use of the MATLAB
% parallel computing toolbox, initialize the multiple worked nodes using
% "matlabpool size X" where "X" is the number of available processing
% nodes.  If the parallel processing toolbox is not available, then replace
% the "parfor" loop with a "for" loop.  This tool is based on the FID-A
% function 'run_simSemiLASERShaped_fast.m'.
% 
% INPUTS:
% To run this function, there are two input arguments:
% spinSys           = spin system to simulate 
% simPars           = structure variable containing the necessary
%                     simulation parameters as follows:
%           simPars.refocWaveform:  sLASER refocusing pulse waveform in FID-A RF structure format
%           simPars.refTp:          duration of refocusing pulses[ms]
%           simPars.flipAngle:      flip angle of refocusing pulses [degrees]
%           simPars.Bfield:         magnetic field strength in [T]
%           simPars.Npts:           number of spectral points
%           simPars.sw:             spectral width [Hz]
%           simPars.lw:             linewidth of the output spectrum [Hz]
%           simPars.thkX:           slice thickness of x refocusing pulse [cm]
%           simPars.thkY:           slice thickness of y refocusing pulse [cm]
%           simPars.fovX:           full simulation FOV in the x direction [cm]
%           simPars.fovY:           full simulation FOV in the y direction [cm]
%           simPars.nX:             number of spatial grid points to simulate in x-direction
%           simPars.nY:             number of spatial grid points to simulate in y-direction
%           simPars.te:             Total sLASER echo time (TE) [ms]
%           simPars.centreFreq:     Centre frequency of simulation [ppm]
%           simPars.lineshape:      Lineshape to simulate ('L' (lorentzian) or 'G' (gaussian))
%
%
% OUTPUTS:
% out               = Simulation results, summed over all space.


function out=compMRS_simsLASER(sys,simPars)

% ************ASSIGN INPUT PARAMETERS TO VARIABLES**********************************
refocWaveform = simPars.refocWaveform;  %sLASER refocusing pulse waveform in FID-A RF structure format
refTp         = simPars.refTp;          %duration of refocusing pulses[ms]
flipAngle     = simPars.flipAngle;      %flip angle of refocusing pulses [degrees]
Npts          = simPars.Npts;           %number of spectral points
sw            = simPars.sw;             %spectral width [Hz]
Bfield        = simPars.Bfield;         %magnetic field strength in [T]
lw            = simPars.lw;             %linewidth of the output spectrum [Hz]
thkX          = simPars.thkX;           %slice thickness of x refocusing pulse [cm]
thkY          = simPars.thkY;           %slice thickness of y refocusing pulse [cm]
fovX          = simPars.fovX;           %full simulation FOV in the x direction [cm]
fovY          = simPars.fovY;           %full simulation FOV in the y direction [cm]
nX            = simPars.nX;             %number of spatial grid points to simulate in x-direction
nY            = simPars.nY;             %number of spatial grid points to simulate in y-direction
te            = simPars.te;             %Total sLASER echo time (TE) [ms]
centreFreq    = simPars.centreFreq;     %Centre frequency of simulation [ppm]
lineshape     = simPars.lineshape;      %Lineshape to simulate ('L' (lorentzian) or 'G' (gaussian))
% ************END OF INPUT PARAMETERS**********************************

%set up spatial grid
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY); %y positions to simulate [cm]

gamma=42577000; %gyromagnetic ratio

% %Load the RF pulse:
% rfPulse=io_loadRFwaveform(refocWaveform,'inv');
rfPulse=refocWaveform;

%If length of RFpulse is >200pts, resample to 100 pts to reduce
%computational workload
if size(rfPulse.waveform,1)>200
    rfPulse=rf_resample(rfPulse,100);
end

%sys=sysRef0ppm
if ~rfPulse.isGM
    %Non-gradient modulated pulse - Calculating the x and y gradient 
    %strengths for the desired slice thickness
    Gx=(rfPulse.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
    Gy=(rfPulse.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
else
    %Gradient modulated pulse
    %1.  Calculating the unitless scaling factor for the GM waveform.
    Gx=(rfPulse.tthk/(refTp/1000))/thkX;
    Gy=(rfPulse.tthk/(refTp/1000))/thkY;
end

%Initialize structures:
% out_posxy_rpc=cell(length(x),length(y),length(ph1));
out_posx_rpc =cell(length(x),1);
d=cell(1,1); %Initialize a cell for the dentity matrix, with elements for each phase cycle;

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%if you do not have the parallel computing toolbox, uncomment the first
%for loop and delete "parfor X=1:length(x)"
parfor X=1:length(x)
%  for X=1:length(x)
        disp(['Executing X-position ' num2str(X) '!!' ]);
        out_posx_rpc{X}=sim_sLASER_shaped_Ref1(Bfield,sys,te,rfPulse,refTp,x(X),Gx,flipAngle,centreFreq);
%                            sim_sLASER_shaped_Ref1(Bfield,sys,te,RF,       tp,  dx, Gx,ph1,   ph2,  centreFreq)
end

%calculate the average density matrix (Doing this inside a separate for
%loop because I couldn't figure out how to do this inside the parfor loop):
for X=1:length(x)
        d{1}=sim_dAdd(d{1},out_posx_rpc{X});
end

% %Initialize structures:
out_posy_rpc =cell(length(x),1);
out=struct([]);

%Now loop through y direction (second refoc pulse only);
parfor Y=1:length(y) %Use this if you do have the MATLAB parallel processing toolbox
%for Y=1:length(y) %Use this if you don't have the MATLAB parallel processing toolbox
        disp(['Executing Y-position ' num2str(Y) '!!' ]);
        out_posy_rpc{Y}=sim_sLASER_shaped_Ref2(d{1},Npts,sw,Bfield,lw,sys,te,rfPulse,refTp,y(Y),Gy,flipAngle,centreFreq,lineshape);
%                            sim_sLASER_shaped_Ref2(d,   n,sw,Bfield,linewidth,sys,te,RF,       tp, dy,  Gy,ph3,    ph4,  centreFreq)
end

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:
for Y=1:length(y)
        out=op_addScans(out,out_posy_rpc{Y});
end


%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together.
numSims=(nX*nY);
out=op_ampScale(out,1/numSims);

%2.  Scale by the total size of the simulated region, relative to the size
%    of the voxel.
voxRatio=(thkX*thkY)/(fovX*fovY);
out=op_ampScale(out,1/voxRatio);


end




%%%%%NESTED FUNCTIONS BELOW%%%%%%%%%

%% Simulate in X-direction only
function d = sim_sLASER_shaped_Ref1(Bfield,sys,te,RF,tp,dx,Gx,flipAngle,centreFreq)

if nargin<9
    centreFreq=2.3;
    if nargin<8
        flipAngle=180;
    end
end

%Check if this is a gradient modulated pulse.  If so, set Gx equal to zero:
if RF.isGM
    %Scale the GM waveform by the factor Gx and then set Gx equal to zero:
    RF=rf_scaleGrad(RF,Gx);
    Gx=0;
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

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%BEGIN sLASER PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                  %EXCITE instantaneously
d=sim_COF(H,d,-1);
d=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dx,Gx);          %1st shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_COF(H,d,1);
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dx,Gx);          %2nd shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_COF(H,d,-1);
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2

end

%% Simulate in Y-direction only
function out = sim_sLASER_shaped_Ref2(d,n,sw,Bfield,linewidth,sys,te,RF,tp,dy,Gy,flipAngle,centreFreq,shape)

if nargin<14
    shape='L';
    if nargin<13
        centreFreq=2.3;
        if nargin<12
            flipAngle=180;
        end
    end
end

%Check if this is a gradient modulated pulse.  If so, set Gy equal to zero:
if RF.isGM
    %Scale the GM waveform by the factor Gy and then set Gy equal to zero:
    RF=rf_scaleGrad(RF,Gy);
    Gy=0;
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

%Calculate Hamiltonian matrices and starting density matrix.
[H]=sim_Hamiltonian(sys,Bfield);

%BEGIN sLASER PULSE SEQUENCE************
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dy,Gy);          %3rd shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_COF(H,d,1);
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dy,Gy);          %4th shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_COF(H,d,-1);
d=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1

[out,~]=sim_readout(d,H,n,sw,linewidth,90,shape);      %Readout along +y axis (90 degree phase);
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
out.flags.isISIS=0;
end
