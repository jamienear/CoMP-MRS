% compMRS_simSPECIAL.m
% Jamie Near, McGill University 2014.
% Revised for CoMP-MRS - Jamie Near and Diana Rotaru, 2025/2026.
% 
% USAGE:
% out=compMRS_simSPECIAL(spinSys,simPars);
% 
% DESCRIPTION:
% This script simulates a localized spin-echo experiment with a fully shaped  
% refocusing pulse.  Phase cycling of the refocusing pulse is performed.  
% Furthermore, simulations are run at various locations in space (1-D) to 
% account for the within-voxel spatial variation of the metabolite signal 
% due to J-evolution.  Summation across phase cycles and spatial positions is
% performed.  As a result of the phase cycling and spatially resolved simulations, 
% this code takes a long time to run.  Therefore, the MATLAB parallel computing
% toolbox (parfor loop) was used to accelerate the siumulations.  Accelration 
% is performed in the direction of the slice selective refocusing pulse.  
% Up to a factor of 12 acceleration can be achieved using this approach.  
% To enable the use of the MATLAB parallel computing toolbox, initialize 
% the multiple worked nodes using "matlabpool size X" where "X" is the 
% number of available processing nodes.  If the parallel processing toolbox
% is not available, then replace the "parfor" loop with a "for" loop.  This 
% tool is based on the FID-A function 'run_simSpinEchoShaped.m'.
% 
% INPUTS:
% To run this function, there are two input arguments:
% spinSys           = spin system to simulate 
% simPars           = structure variable containing the necessary
%                     simulation parameters as follows:
%           simPars.refocWaveform:  SPECIAL refocusing pulse waveform in FID-A RF structure format
%           simPars.refTp:          duration of the refocusing pulse[ms]
%           simPars.flipAngle:      flip angle of the refocusing pulse [degrees]
%           simPars.Bfield:         magnetic field strength in [T]
%           simPars.Npts:           number of spectral points
%           simPars.sw:             spectral width [Hz]
%           simPars.lw:             linewidth of the output spectrum [Hz]
%           simPars.thkX:           slice thickness of x refocusing pulse [cm]
%           simPars.fovX:           full simulation FOV in the x direction [cm]
%           simPars.nX:             number of spatial grid points to simulate in x-direction
%           simPars.tau1:           SPECIAL echo time (TE) in [ms]
%           simPars.centreFreq:     Centre frequency of simulation [ppm]
%           simPars.lineshape:      Lineshape to simulate ('L' (lorentzian) or 'G' (gaussian))
%
% OUTPUTS:
% out               = Simulation results, summed over all space.

function [out,out_pos]=compMRS_simSPECIAL(spinSys,simPars);

% ************ASSIGN INPUT PARAMETERS TO VARIABLES**********************************
refocWaveform = simPars.refocWaveform;  %SPECIAL refocusing pulse waveform in FID-A RF structure format
refTp         = simPars.refTp;          %duration of the refocusing pulse[ms]
flipAngle     = simPars.flipAngle;      %flip angle of the refocusing pulse [degrees]
Npts          = simPars.Npts;           %number of spectral points
sw            = simPars.sw;             %spectral width [Hz]
Bfield        = simPars.Bfield;         %magnetic field strength in [T]
lw            = simPars.lw;             %linewidth of the output spectrum [Hz]
thkX          = simPars.thkX;           %slice thickness of x refocusing pulse [cm]
fovX          = simPars.fovX;           %full simulation FOV in the x direction [cm]
nX            = simPars.nX;             %number of spatial grid points to simulate in x-direction
tau1          = simPars.tau1;           %SPECIAL echo time (TE) in [ms]
centreFreq    = simPars.centreFreq;     %Centre frequency of simulation [ppm]
lineshape     = simPars.lineshape;      %Lineshape to simulate ('L' (lorentzian) or 'G' (gaussian))
% ************END OF INPUT PARAMETERS**********************************

pos=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]

phCyc=[0,90]; %phase cycling steps for 1st refocusing pulse [degrees]

% %Load RF waveforms
% RF=io_loadRFwaveform(refocWaveform,'ref',0);
RF=refocWaveform;

gamma=42577000; %gyromagnetic ratio

%If length of RF pulse is >200 pts, resample to 100 pts to reduce
%computational workload
if size(RF.waveform,1)>200
    RF=rf_resample(RF,100);
end

Gx=(RF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]

%n=1;
%totalIters=length(x)*length(y)*length(editPhCyc1)*length(editPhCyc2)*length(refPhCyc1)*length(refPhCyc2);

%Initialize structures:
out_pos_rpc=cell(length(pos),length(phCyc));
out_pos=cell(length(pos),1);
out=struct([]);

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%for X=1:length(x);
parfor X=1:length(pos);
                for RP=1:length(phCyc)
                        disp(['Executing X-position ' num2str(X) ' of ' num2str(length(pos)) ', '...
                            'Refoc phase cycle ' num2str(RP) ' of ' num2str(length(phCyc)) '!!!']); 
                        out_pos_rpc{X}{RP}=sim_spinecho_shaped(Npts,sw,Bfield,lw,spinSys,tau1,RF,refTp,Gx,pos(X),phCyc(RP),flipAngle,centreFreq,lineshape);
                        if RP==1
                            out_pos{X}=out_pos_rpc{X}{RP};
                        else
                            out_pos{X}=op_subtractScans(out_pos{X},out_pos_rpc{X}{RP});
                        end
                    end %end of 1st refocusing phase cycle loop
        out=op_addScans(out,out_pos{X});
end %end of spatial loop (parfor) in x direction.
        

%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together.
numSims = nX * length(phCyc);
out=op_ampScale(out,1/numSims);

%2.  Scale by the total size of the simulated region, relative to the size
%    of the voxel.
voxRatio=(thkX)/(fovX);
out=op_ampScale(out,1/voxRatio);





       