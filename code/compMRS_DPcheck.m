%compMRS_DPcheck.m
%Jamie Near, Sunnybrook Research Institute, 2025
%
% USAGE:
% [check]=compMRS_DPcheck(DPid);
%
% DESCRIPTION:
% Simple script to check some useful information about a data packet in
% CoMP-MRS.  
%
% The function returns a scructure called 'check', which contains the
% following fields:
%
% INPUTS:
% DPid:     Data Packet ID (i.e. DP01, DP02, etc.).
% 
% OUTPUTS:
% check:    Structure with the following fields:.
%               -nSubj:  Number of subjects in the data packet
%               -nSes:  Number of sessions per subject
%               -vendor:  What vendor did this DP come from 
%               -version: What vendor software version was the DP generated with
%               -seq:  What pulse sequence was the DP generated with
%               -allSame:  Boolian (True/False) to check if all subjects and sessions in the Data 
%               Packet are using the same MRI vendor and software version.  


function [check]=compMRS_DPcheck(DPid)

%Check the number of subjects
subs=dir([DPid '/sub-*']);
check.nSubj=length(subs);

%Check the number of sessions per subject
for n=1:check.nSubj;
    ses=dir([DPid filesep subs(n).name filesep 'ses-*']);
    check.nSes(n)=length(ses);
    for m=1:check.nSes(n)
        svsDir=dir([DPid filesep subs(n).name filesep ses(m).name filesep 'mrs' filesep '*_svs']);
        if exist([DPid filesep subs(n).name filesep ses(m).name filesep 'mrs' filesep svsDir(length(svsDir)).name filesep 'fid']) && exist([DPid filesep subs(n).name filesep ses(m).name filesep 'mrs' filesep svsDir(length(svsDir)).name filesep 'procpar'])
            %Fill in the vendor and version fields (VARIAN DATA):
            check.vendor{n,m}='VARIAN';
            par=readprocpar([DPid filesep subs(n).name filesep ses(m).name filesep 'mrs' filesep svsDir(length(svsDir)).name filesep 'procpar']);
            check.version{n,m}=par.parversion;
            check.seq{n,m}=par.seqfil;

        elseif exist([DPid filesep subs(n).name filesep ses(m).name filesep 'mrs' filesep svsDir(length(svsDir)).name filesep 'method']) && exist([DPid filesep subs(n).name filesep ses(m).name filesep 'mrs' filesep svsDir(length(svsDir)).name filesep 'acqp'])
            %Fill in the vendor and version fields (BRUKER DATA):
            check.vendor{n,m}='BRUKER';
            headerAcqp = parseBrukerFormat([DPid filesep subs(n).name filesep ses(m).name filesep 'mrs' filesep svsDir(length(svsDir)).name filesep 'acqp']);
            check.version{n,m}=headerAcqp.ACQ_sw_version;
            check.seq{n,m}=headerAcqp.PULPROG;
            

        end
        
    end
end

%Check that all scans within the DP have the same vendor and version:
check.allSame = isequal(check.vendor{:}) && isequal(check.version{:}) && isequal(check.seq{:});

end









function header = parseBrukerFormat(inputFile)
% This subroutine uses regular expressions and case differentiations to
% extract all relevant information from a Bruker-formatted header file
% (acqp, method, etc.)

% Open file
fid = fopen(inputFile);

% Get first line
tline = fgets(fid);

% Loop over subsequent lines
while ~feof(fid)

    % First, get the parameters without a $
    [tokens, matches] = regexp(tline,'##([\w\[\].]*)\s*=\s*([-\(\w\s.\"\\:\.,\)]*)','tokens','match');

    % When a matching string is found, parse the results into a struct
    if length(tokens) == 1

        fieldname = regexprep(tokens{1}{1}, '\[|\]',''); % delete invalid characters

        % Convert numbers to doubles, leave strings & empty lines alone
        if ~isnan(str2double(tokens{1}{2}))
            value = str2double(tokens{1}{2});
        else
            value = strtrim(tokens{1}{2});
        end

        % Convert char to string
        if ischar(value)
            value = string(value);
        end

        % Store
        header.(fieldname) = value;

        % Get next line
        tline = fgets(fid);
        continue

    else

        % If not a match, get the parameters with a $
        [tokens, ~] = regexp(tline,'##\$([\w\[\].]*)\s*=\s*([-\(\w\s.\"\\:\.,\)]*)','tokens','match');


        % When a matching string is found, parse the results into a struct
        if length(tokens) == 1

            fieldname = regexprep(tokens{1}{1}, '\[|\]',''); % delete invalid characters

            % Determine if the value indexes an array (signaled by a number
            % inside a double bracket, e.g. ##$PULPROG=( 32 )), or a single
            % value (signaled by just a string, e.g. ##$ACQ_user_filter_mode=Special)
            [tokensValue, ~] = regexp(tokens{1}{2},'\( (.*) \)','tokens','match');

            % If there's a match, we need to parse the subsequent lines
            % which contain the array
            if length(tokensValue) == 1

                % Arrays can span more than one line, unfortunately, so we
                % need to do some clever pattern matching - basically, we
                % want to extract lines until they lead with ## or $$
                % again:
                endOfBlock = 0;
                multiLine  = '';
                while endOfBlock ~=1

                    % Get next line
                    tline = fgets(fid);

                    if contains(string(tline), ["$$", "##"])
                        endOfBlock = 1;
                    else
                        multiLine = [multiLine, tline];
                    end

                end

                % If the line is bracketed by <>, store that contents as
                % one
                contents = {};
                [tokensBrackets, ~] = regexp(multiLine,'<(.*)>\n','tokens','match');
                if length(tokensBrackets) == 1
                    contents{1} = tokensBrackets{1}{1};

                    % Convert numbers to doubles, leave strings & empty lines alone
                    if ~isnan(str2double(contents{1}))
                        value = str2double(contents{1});
                    else
                        value = strtrim(contents{1});
                    end

                    % Convert char to string
                    if ischar(value)
                        value = string(value);
                    end

                else
                    % If not, it's an array.
                    % Sometimes this array can even contain vectors, for example TPQQ
                    % In this case, let's look for recurring parentheses
                    % again:
                    multiLine = erase(multiLine, newline); % remove new line characters
                    [tokensParentheses, ~] = regexp(multiLine,'\(([^\)]+)\)','tokens','match');

                    if ~isempty(tokensParentheses)
                        for rr = 1:length(tokensParentheses)
                            test = textscan(tokensParentheses{1,rr}{1}, '%s', 'Delimiter', ',');
                            contents{rr} = test{1};
                        end
                    else
                        % use textscan to convert space-delimited vectors to cell array
                        test = textscan(multiLine, '%s');
                        contents = test{1};
                    end

                    % Convert numbers to doubles, leave strings & empty lines alone
                    if ~isnan(str2double(contents))
                        value = str2double(contents);
                    else
                        value = strtrim(contents);
                    end
                    
                    % Convert char to string
                    if ischar(value)
                        value = string(value);
                    end

                    if iscell(value) && length(value) == 1
                        value = value{1};
                    end

                end

                % Store
                header.(fieldname) = value;

                continue

            else

                % Convert numbers to doubles, leave strings & empty lines alone
                if ~isnan(str2double(tokens{1}{2}))
                    value = str2double(tokens{1}{2});
                else
                    value = strtrim(tokens{1}{2});
                end

                % Convert char to string
                if ischar(value)
                    value = string(value);
                end

                if iscell(value) && length(value) == 1
                    value = value{1};
                end

                % Store
                header.(fieldname) = value;

            end

            % Get next line
            tline = fgets(fid);
            continue

        else

            % Get next line
            tline = fgets(fid);
            continue

        end

    end

end

fclose(fid);

end







