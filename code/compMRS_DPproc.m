%compMRS_DPproc.m
%Jamie Near, Sunnybrook Research Institute, 2025
%Diana Rotaru, Columbia University, 2025
%
% USAGE:
% [out,outw]=compMRS_DPproc.m


function [out,outw]=compMRS_DPproc(DPid)
try
    [check]=compMRS_DPcheck(DPid);
    if check.allSame
        switch char(string(check.vendor(1)))
            case 'BRUKER'
                [out, outw] = compMRS_DPload(DPid);

                for m = 1:check.nSubj
                    for n = 1:check.nSes(m)
                        [out_proc{m,n}]=compMRS_procBrukerRaw(out{m,n},outw{m,n});
                    end
                end
            case 'VARIAN'
                disp('Varian')
        end

    end
catch
    disp([DPid ' error'])
end

end


