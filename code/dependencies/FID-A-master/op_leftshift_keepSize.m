% op_leftshift_keepSize.m
% Jamie Near, McGill University 2014.
% Thanh Phong Lê, 2025
% 
% USAGE:
% out=op_leftshift_keepSize(in,ls);
% 
% DESCRIPTION:
% Remove leading datapoints from the fid to get rid of 1st order phase
% errors, and fill zeros at the end such that the size does not change.
% 
% INPUTS:
% in     = input data in matlab strucure format.
% ls     = number of points to remove from the beginning of the fid.
%
% OUTPUTS:
% out    = Output following left shifting.  

function out=op_leftshift_keepSize(in, ls)

if in.flags.leftshifted
    cont=input('WARNING:  Left shifting has already been performed!  Continue anyway?  (y or n)','s');
    if cont=='y'
        %continue;
    else
        error('STOPPING');
    end
end

fids = in.fids;
fids(1:end-ls,:) = fids(1+ls:end,:);
fids(end-ls:end,:) = 0;

specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);
    
%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.leftshifted=1;
