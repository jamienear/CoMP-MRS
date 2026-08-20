% op_plotMeanSD.m
% Jamie Near, Sunnybrook Research Institute 2026.
% 
% USAGE:
% out=op_plotMeanSD(in);
% 
% DESCRIPTION:
% Plot a standard mean +/- SD band visualization for a series of spectra.
% Note that the input spectra in the cell array should be aligned and
% should have the same dimensions, otherwise, this function might fail.  
% 
% INPUTS:
% in     = a cell array whose elements are the input spectra in matlab 
%          structure format.
% ppmmin = lower limit of ppm scale to plot (optional.  Default = 0.2 ppm).
% ppmmax = upper limit of ppm scale to plot (optional.  Default = 5.2 ppm).
% xlab   = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
% ylab   = label for the y-axis (optional.  Default = '');
% tit    = label for the title of the plot (optional.  Default = '');
% fact   = scaling factor for SD
%
% OUTPUTS:
% out    = Figure handle.

function [out, mu, sigma] =op_plotMeanSD(in,ppmmin,ppmmax,xlab,ylab,tit,fact)

if nargin<7
    fact = 1;
    if nargin<6
        tit='';
        if nargin<5
            ylab='';
            if nargin<4
                xlab='Frequency (ppm)';
                if nargin<3
                    ppmmax=5.2;
                    if nargin<2
                        ppmmin=0.2;
                        if nargin<1
                            error('ERROR: no input spectrum specified.  Aborting!!');
                        end
                    end
                end
            end
        end
    end
end

%First concatenate the elements of the input into a single FID-A structure
%that treats the spectra as "averages"
avgs=in{1};
avgs.flags.averaged=0;
for n=2:length(in)
    avgs=op_concatAverages(avgs,in{n});
end

%extract the and 'specs' field so that we can calculate mean and SD:
specs=avgs.specs;
out.specs = specs;

%also extract the 'ppm' field:
ppm=avgs.ppm';
out.ppm = ppm;

%Calculate the mean:
mu = mean(real(specs),2);

%Calculate the SD:
sigma = std(real(specs),0,2);

% %Make the plot:
% figure; hold on;
% %Shaded area
% fill([ppm; flipud(ppm)],...
%     [(mu+fact*sigma); flipud(mu-fact*sigma)], ...
%     [0.8 0.8 1],...    % light blue color
%     'EdgeColor','none',...
%     'FaceAlpha',0.8);   % Transparency
% 
% %Mean curve
% plot(ppm,mu,'b','LineWidth',2);
% 
% set(gca,'XDir','reverse'); %Flip x-axis for MRS spectrum convention
% xlim([ppmmin ppmmax]); %Limit the ppm display range
% set(gcf,'Color','w') %Figure background white
% set(gca,'ycolor','w');  %Make the y-axis invisible
% xlabel('Frequency (ppm)');  %Label x-axis




        
    
    

