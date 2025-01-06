% Test suite for new FID-A Bruker loaders
% This function sets up a series of tests to see if the function
% io_loadspec_bruk_new fails on a series of test data.
%
% Georg Oeltzschner 2025

classdef io_loadspec_bruk_newTest < matlab.unittest.TestCase

    % Load the list of files to be tested
    properties (TestParameter)
        testDataFiles = getTestFiles();
    end

    % List of test methods
    methods(Test)

        % Test 01: Bruker-combined data
        function testCombined(testCase, testDataFiles)
            [testOut, testRef] = testLoadCombined(testDataFiles);

            % Verify that the test output is not zero
            testCase.verifyNotEqual(testOut, 0);

            % Verify that the test ref output is not zero
            testCase.verifyNotEqual(testRef, 0);
        end

        % Test 02: Raw (uncombined) data
        function testRaw(testCase, testDataFiles)
            [testOut, testRef] = testLoadRaw(testDataFiles);

            % Verify that the test output is not zero
            testCase.verifyNotEqual(testOut, 0);

            % Verify that the test ref output is not zero
            testCase.verifyNotEqual(testRef, 0);

        end

    end
end

%% Helper functions below

function testDataFiles = getTestFiles()
% Function to create the list of files to be tested

% Get the path of this test function
testFunctionPath = fileparts(mfilename('fullpath'));
% Divide path into parts
pathparts = strsplit(testFunctionPath, filesep);
% Assemble new path
mainPath = fullfile(pathparts{1:end-2});

testDataFiles = {...
    fullfile(mainPath, 'data', 'DP01', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-press_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP03', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-press_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP04', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-press_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP05', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-press_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP06', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-slaser_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP07', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-slaser_voi-stri_svs'), ...
    fullfile(mainPath, 'data', 'DP08', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-slaser_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP09', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-press_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP10', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_R-Hipppo_svs'), ...
    fullfile(mainPath, 'data', 'DP14', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-steam_voi-stri_svs'), ...
    fullfile(mainPath, 'data', 'DP15', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-semiLASER_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP16', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-semiLASER_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP17', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-press_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP18', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-laser_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP19', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-press_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP20', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-steam_voi-stri_svs'), ...
    fullfile(mainPath, 'data', 'DP21', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-press_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP22', 'sub-01', 'ses-01', 'mrs', 'met'), ...
    fullfile(mainPath, 'data', 'DP23', 'sub-01', 'ses-01', 'mrs', 'met'), ...
    fullfile(mainPath, 'data', 'DP24', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-steam_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP25', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-press_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP26', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses-01_acq-press_voi-hipp_svs'), ...
    fullfile(mainPath, 'data', 'DP27', 'sub-01', 'ses-01', 'mrs', 'sub-01_ses_01_acq_press_voi-hipp_svs'), ...
    };

end



function [out, ref] = testLoadCombined(testDataFile)
% Function to test the Bruker-combined data
try
    [out, ref] = io_loadspec_bruk_new(testDataFile, 'n');
catch
    % If load fails, return zero flag
    out = 0;
    ref = 0;
end

end

function [out, ref] = testLoadRaw(testDataFile)
% Function to test the Bruker-combined data
try
    [out, ref] = io_loadspec_bruk_new(testDataFile, 'y');
catch
    % If load fails, return zero flag
    out = 0;
    ref = 0;
end

end