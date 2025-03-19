%





function basis = compMRS_makeBasis(DPid)


%First check the vendor using DPcheck:
check = compMRS_DPcheck(DPid);

%Get the method file from the first subject and first session.  Assume all method files after that are the same.  
if strcmp(check.vendor(1),'BRUKER')
    subjs=dir([DPid filesep 'sub*']);
    sess=dir([DPid filesep subjs(1).name filesep 'ses*']);
    %Find the MRS data path and the REF data path:
    svspath = dir([DPid filesep subjs(1).name filesep sess(1).name filesep 'mrs' filesep '*svs']);
    %read the method file:
    method = compMRS_loadMethod([svspath filesep 'method']);

    %Now extract the relevant params from the method file:
    params.seq=method.Method;
    params.TE1 = method.TE1;
    params.TE2 = method.TE2;
    params.refocWaveform = ; %TBD
    params.refTp = extractBetween(method.VoxPulse2, '(', ',');
    
