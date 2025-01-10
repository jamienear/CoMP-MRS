# CoMP-MRS

Repository of processing tools and procedures for a multi-site preclinical MRS project

## Update Log (for our internal development)

### Jan 10 2025

Major update for raw data. Main innovations:

- Added (and tested) functionality to load raw (coil-uncombined) data and store the multi-dimensional arrays correctly. Note that it appears from the available test data that PV5 does not store separate channel data.
- Navigator scans (PV5, if available) are now returned as a separate `nav` struct, not `ref`.
- Data with multiple repetitions (which are an outer loop around averages) are now correctly loaded (at least in 'combined' mode), with the `subSpecs` dimension being used to store the repetition dimension.
- Parse the sequence method string from the `PULPROG` variable in `ACQP` file.
- Refactored the FID data load function.

We still need to figure out where to get the proper coil combination coefficients from:

- Phases are in `PVM_PhasedArray`
- Some datasets have values in `PVM_EncChanScaling` (but not sure they are what I think they are)
- Some datasets (definitely not all) seem to have a water signal stored in `PVM_RefScan` with separate coil elements. This might be able to be used.

Notes on scaling - getting that right is going to be very important:

- It seems that (at least for PV5), combined data is simply the **sum** of the individual transients, **not the average**.
- Receiver gain is not loaded at all yet (Jessie does it in her code), but it can differ between water reference and metabolite scans.
  
### Jan 5 2025

First alpha version of the new reader. Main innovations compared to the 'old' reader are:

- TX freqs for PV6 onwards are explicitly stored separately for metabolite and ref scans. headerACQP.BF1 is identical to headerMethod.PVM_FrqWork (and is slightly different from headerMethod.PVM_FrqRef). For PV6 onwards, I am now reading txfrq_ref and txfrq from these two fields and shift the time-domain signals accordingly (but see open questions below)
- Center frequency is *explicitly* stored in the headers from PV6 onwards so I am reading it out and using it to calculate the ppm axis - mostly it's been 4.7 so we'll assume that for PV5 too. For compatibility with the other FID-A functions, we should probably store this centerfreq value in the FID-A header? (Note: I don't think we should end up with different ppm axes for ref and metabolite scans, so I'm just using the frequency shifts from above)

I have tested the 'combined data loading' on Bruker sets between DP01-DP32 (with updated test suite). Most test datasets I've tried load without error (except DP22 and DP23 when trying to load in 'combined' mode, because don't have a `fid` file - I suppose this is because it's a non-standard sSPECIAL sequence? We can either let this throw an error or re-direct the function to open the `ser` file).

Visual inspection of the combined data has appeared reasonable, but some doubts remain:

- DP17 doesn't seem to be shifted right (residual water at 6.4 ppm?), so I am not totally certain I have the frequency/shift referencing nailed perfectly.
- Some datasets (DP15, DP16, DP27) have NAA appear to the left of 2.01 ppm (e.g. 1.88, 1.95 ppm), but the residual water is where it should be. I was tearing my head out whether I'm calculating the ppm axis wrong, but I now think it's a temperature effect? Unfortunately the spreadsheet does not have the temperature for these two DPs.
- I was also under the impression that PV6 onwards *always* store a combined reference scan, but that does not appear to be the case (DP05, DP08, DP09, DP10, DP15, DP16, DP17, DP18)

I have, so far, not visually inspected the output of the 'uncombined' (raw) data loading.

### Jan 3, 2025 (GO)

I have created a new 'parallel' development version `io_loadspec_bruk_new.m` alongside a new test suite `io_loadspec_bruk_newTest`. The new reader has a completely rewritten routine to extract header information from acqp, acqus, and method files, allowing these to be quite flexibly plucked from a struct. The new reader successfully runs the current test suite (15 datasets between DP01 and DP19), i.e., it loads the datasets without error (both with and without the 'rawdata' flag).

### Open Questions

- (Jan 5 GO) Are the `fid` (combined and processed) data already eddy-current-corrected? There appears to be a parameter `EDC_OnOff` indicating this but it does not appear in all data.
- (Jan 5 GO) For now, I'm applying the frequency shift (corresponding to the difference between FrqRef and FrqWork) to the *water-suppressed* data (PV-360 only) or to the ref data (if not PV-360), but I'm not confident that this is correct.
- (Jan 5 GO) Is cutting off the number of points indicated in GRPDLY really sufficient? When I look at the FIDs, I frequently find that the top of the echo appears slightly after that chop-off point

### TO DO

- [ ] (Jan 7 GO) If combined data are requested but not fid file is found, look for ser file (and return that)
- [x] (Jan 5 GO) Test the extended dataset (DP20-DP32)
- [x] (Jan 5 GO) Test raw data loading
- [ ] (Jan 3 GO) Design tests to check that data are actually loaded *correctly*, i.e., do some sanity checks, visual inspection, etc.
- [ ] (Jan 3 GO) Design tests for reference data (in the `mrsref` directories).
- [ ] (Jan 3 GO) I have found some interesting code in Jessie's loader that does some scaling according to the receiver gain - need to look into that
- [ ] (Jan 3 GO) Need to figure out where to find amplitude and phase coefficients for coil combination

## Contents

`code\notes` contains helpful documentation to understand the various Bruker file formats.

`code\notes\loaders` contains various bits and pieces of code used to load data from different ParaVision versions (courtesy of Brayan Alves, Jessie Mosso, Diana Rotaru).

`code\notes\io_loadspec_bruk.m` is the latest version of the FID-A Bruker loader. We will develop this into a new 'one loads all' function that works for all DPs.

### Test suite

`code\tests\io_loadspec_brukTest.m` contains a test suite that test-runs the `io_loadspec_bruk.m` function on one dataset per DP. Tests currently include two function calls per dataset: 1) for coil-combined data; 2) for coil-uncombined data.

#### Test data

Testing requires the complete contents of the CoMP-MRS test data package stored in a folder `data` at the level of the main repository. All DP folders and sub-folders *must* follow BIDS terminology for the test to function properly; e.g., DP06 has another sub-folder DP06 whose contents must be moved to the parent DP06 directory. At the time of writing, this tests runs for DP01 to DP19.

NOTE:  Today I removed all instances of nested DPXX folders (i.e. DP06/DP06/...) from both the test data package and the complete CoMP-MRS data repository.  I also added new 
datasets to the test data package, now up to DP32. %JN - 30 Dec 2024 

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
