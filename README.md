# CoMP-MRS

Repository of processing tools and procedures for a multi-site preclinical MRS project

## Contents

`code\notes` contains helpful documentation to understand the various Bruker file formats.

`code\notes\loaders` contains various bits and pieces of code used to load data from different ParaVision versions (courtesy of Brayan Alves, Jessie Mosso, Diana Rotaru).

`code\notes\io_loadspec_bruk.m` is the latest version of the FID-A Bruker loader. We will develop this into a new 'one loads all' function that works for all DPs.

### Test suite

`code\tests\io_loadspec_brukTest.m` contains a test suite that test-runs the `io_loadspec_bruk.m` function on one dataset per DP. Tests currently include two function calls per dataset: 1) for coil-combined data; 2) for coil-uncombined data.

#### Test data

Testing requires the complete contents of the CoMP-MRS test data package stored in a folder `data` at the level of the main repository. All DP folders and sub-folders *must* follow BIDS terminology for the test to function properly; e.g., DP06 has another sub-folder DP06 whose contents must be moved to the parent DP06 directory. At the time of writing, this tests runs for DP01 to DP19.

NOTE:  Today I removed all instances of nested DPXX folders (i.e. DP06/DP06/...) from both the test data package and the complete CoMP-MRS data repository. %JN - 30 Dec 2024 

```md
.git
code
data
└── DP01
    └── sub-01
        └── ses-01
            └── mrs
                └── sub-01_ses-01_acq-press_voi-hipp_svs
└── DP02
└── DP03
└── ...
.gitignore
makeDirs.sh
README.md [this file]
```

#### Running the test suite

The test suite can be run and viewed from the MATLAB command line with the following command:

```matlab
result = runtests('io_loadspec_brukTest')
table(result)
```
