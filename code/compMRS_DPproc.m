%compMRS_DPproc.m
%Jamie Near, Sunnybrook Research Institute, 2025
%Diana Rotaru, Columbia University, 2025
%
% USAGE:
% [out,outw]=compMRS_DPproc.m


function [out,outw]=compMRS_DPproc(DPid)

[check]=compMRS_DPcheck(DPid);
if check.allSame
    switch char(string(check.vendor(1)))
        case 'BRUKER'
            [out, outw] = compMRS_DPload(DPid);
            for i=1:length(out)
                outLB = op_filter(out{i,1},1);
                outAlign = op_alignAverages(outLB);
                outAverage = op_averaging(outAlign);
                op_plotspec(outAverage);
            end
        case 'VARIAN'
            
            end
    end
end

