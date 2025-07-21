# CoMP-MRS

Repository of processing tools and procedures for a multi-site preclinical MRS project

## Update Log (for our internal development)

### Apr 21 2025

Some progress these days: except for one dataset (DP14), all frequency shifts and phase differences between online and offline processing seem to be [resolved](/code/notes/graphics/20250421_Comparison_Processing.pdf).

Here are the main observations:
- Multi-channel data require ECC to be applied before coil combination (as brought up previously by Georg), using the **uncombined** reference scan (stored as parameter `PVM_RefScan` in the method file). Following that, the channels can be combined using their respective sensitivity factor, **without** any further phase adjustment. This improves DP04, DP18, DP24, DP25.
- The differences previously observed in DP04 and DP09 were because the pdata folder was empty. The code thought that the scanner-processed data were not ECCed, and applied it again, causing the discrepancy.
- The online reconstruction of DP1 (PV5.1) and DP20 (PV6.0.1) was performed with retro frequency lock, so it must be applied to our processing as well for the purpose of the comparison (as brought up previously by Georg). It is now implemented using the FID-A functions `op_freqAlignAverages` and `op_makeFreqDrift` (note that this feature still needs to be tested with PV7 and PV360 data).

I updated `io_loadspec_bruk.m` and `processBrukerRaw.m` accordingly.

Other observations:
- A perfect match is always obtained only when scanner-processed-data were not ECCed.
- DP01, DP05, DP10, DP15, DP16 have no ECC reference (the reference scan is performed in another acquisition (the `_mrsref` folder), so they are not phased for now.
- There is a slight difference between the ECC algorithm used here and Bruker’s one, which is causing the noise-level differences. Bruker seems to use a filter, with a parameter `EdcFilterWidthHz` (in the methreco file). Description from the manual: "`EdcFilterWidthHz` - Specifies the width of the Gaussian filter for eddy current compensation.”
- For DP14, there is a remaining hump around 2.9-3.6ppm, which I am not able to correct.


### Feb 07 2025

Big progress today:

- I figured out the mystery scaling incoherences, thanks to the very useful notes from Thanh!
  1) I removed the receiver gain scaling (it should *NOT* be applied anywhere);
  2) used the coil weighting in `encChanScaling`;
  3) did either coil averaging or summation (depending on the software version);
  4) *averaged* the raw data across the averages/transients dimension, independent of software version. 
- I also added Klose eddy-current correction to the `processBrukerRaw` script (which will eventually end up as an Osprey `osp_loadBruker` wrapper).
- **For many datasets, reconstruction from the raw data agrees with the scanner-processed data either *perfectly* ([DP06](/code/notes/graphics/DP06_perfect_020725.png),07,08,10) or to within noise (DP04,05,15,16,19,21,[26](/code/notes/graphics/DP26_perfect_020725.png),27).**
- There are a few remaining datasets with phase and frequency shift differences resulting in structural differences that I haven’t quite figured out (DP01,09,14,20,24,25).
  - Do we, maybe, need to apply eddy-current correction to each channel *before* combining the coils? (see e.g. DP14,18,24,25 - all of these have `Edc_OnOff = Yes`)
  - Do we, maybe, need to apply the interleaved frequency correction update `RetroFrequencyLock_OnOff`? Some datasets have `PVM_RefScan` and `PVM_RefScanPC` (PC = phase correction?), others `PVM_DriftCompYesNo`
- Finally, some datasets (even scanner-processed ones) have a, sometimes huge, zero-order phase shift. It’s nothing we can’t fix in post-processing and modeling, but maybe we’re overlooking one final phase correction parameter somewhere (it may be related to the other phase issues).

#### Notes per dataset

- DP01: Processed and scanner spectra agree to within noise except for minuscule frequency shift; both have ~45 deg phase; no ECC reference. 
  - If I load `fid.orig` as the scanner-processed file instead, the outcome is *identical*!
- DP04: Processed and scanner spectra agree to within noise but are approximately 170 degrees out of phase 
- DP05: Processed and scanner spectra agree to within noise; both have ~160 deg phase; no ECC reference
- DP06: Processed and scanner spectra agree are perfectly *identical* and both are perfectly phased after ECC
- DP07: Processed and scanner spectra agree are perfectly *identical* and both are perfectly phased after ECC
- DP08: Processed and scanner spectra agree are perfectly *identical* and both are perfectly phased; no ECC reference
- DP09: Processed and scanner spectra agree to within noise except for a small frequency/phase shift; both appear overall well-phased after ECC
- DP10: Processed and scanner spectra agree are perfectly *identical* (only 1 coil present though, so no coil combination); both have ~160 deg phase; no ECC reference
- DP14: Processed and scanner spectra agree to within 10%; somewhat of a difference around the water base; both are perfectly phased after ECC
- DP15: Processed and scanner spectra agree to within noise; both have ~180 deg phase; no ECC reference
- DP16: Processed and scanner spectra agree to within noise; both have ~180 deg phase; no ECC reference
- DP17: Super noisy dataset; even the ref is noisy; I basically discard this as useless
- DP18: Processed and scanner spectra agree to within 10%; phase difference of about maybe 10%?; both are well phased; no ECC reference ALTHOUGH Edc_OnOff is set (MYSTERY, this shouldn’t be able to happen)
- DP19: Processed and scanner spectra agree to within noise; both are perfectly phased after ECC
- DP20: Processed and scanner spectra agree to within noise except for minuscule frequency shift; both are perfectly phased after ECC
- DP21: Processed and scanner spectra agree to within noise; both are perfectly phased after ECC
- DP22: SPECIAL data (not yet supported, working on it)
- DP23: SPECIAL data (not yet supported, working on it)
- DP24: Processed and scanner spectra agree to within 10%; but approximately 30 Hz frequency and -75 deg phase offset between them (scanner spectrum perfectly phased); no ECC reference ALTHOUGH Edc_OnOff is set (MYSTERY, this shouldn’t be able to happen)
- DP25: Processed and scanner spectra agree to within 10%; but approximately 4 Hz frequency and 152 deg phase offset between them (scanner spectrum perfectly phased); no ECC reference ALTHOUGH Edc_OnOff is set (MYSTERY, this shouldn’t be able to happen)
- DP26: Processed and scanner spectra agree to within noise; both are perfectly phased after ECC
- DP27: Processed and scanner spectra agree to within noise; both have approximately 20 deg phase after ECC

### Jan 23 2025

Update for the raw data loader about eddy-current correction. Details:

- The loader returns an additional flag `isECCed` that indicates whether the **processed data** have been eddy-current-corrected (this is determined by the checkbox on the Bruker console).
- The governing parameter `EDC_OnOff` is stored in different places (apparently in method file in PV5 or in methodreco file in other versions) - if not found anywhere, it is assumed to be set to "Off"
- Oddly enough, some datasets have this parameter set to "On", but do not have a reference dataset (according to the manual, the checkbox can only be set if a reference scan is acquired)
- In any case, this extra flag can be used in a subsequent script to dynamically decide whether ECC needs to be applied or not.
So far, a few findings:
- If the dataset has reference scan and `EDC_OnOff` set to `Off`, applying `op_eccKlose` seems to perfectly phase the data (ex: DP19, DP20)
- If the dataset has `EDC_OnOff` set to `On`, the metabolite data are ECC-ed and perfectly phased, but the **reference data are not** (ex: DP14, DP17, DP18, DP21, DP24, DP25, DP26, DP27)
- If the dataset has **no** reference scan and `EDC_OnOff` set to `Off`, things are a little more difficult:
  - DP08 still seems perfectly phased
  - DP04, DP09, DP16, DP16 are not
  - In some cases (DP15), applying the *last element* of `PVM_ArrayPhase` appears to fix the overall phase, but that's definitely not the case consistently.

### Jan 17 2025

Updates for the raw data loader. Details:

- For multi-channel data, the phases are read from `PVM_PhasedArray` and returned in a field `ph` of a new output struct `coilcombos` which can be fed into the FID-A coil combination function `op_addrcrvs.m`. 
- Coil combination *amplitudes* are still a mystery.  My hunch is that they just add the (phased) coil signals up (or average them), but I can't get the scaling between raw and processed data consistent.
- Receiver gain for the water-suppressed data is read from `ACQ_jobs` in the ACQP file (according to the Bruker manual, the `RG` value in the same file is *not* correct) and applied to the raw data.

By applying the coil combination phases to the raw data (and then simply summing the averages up), one can get remarkably close to the Bruker-processed data, but only to a scaling factor that I cannot wrap my head around.
Between the coil combination amplitudes and the receiver gain, we must still be missing something. I think that, at this point, we should e-mail Bruker about the details of the on-scanner recon.

Other additions include a script `processBrukerRaw.m` in `code\tests` that can be used to compare the Bruker-processed spectrum with one that we process with FID-A functions from the raw data.
This function requires a couple of FID-A functions (`op_addrcvrs`, `op_ampScale`, `op_averaging` and `phase`) that I have added for convenience. I've also added my own FID-A-style function for a classic Klose eddy current correction (`op_eccKlose`).

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

- (Jan 10 GO) How does Bruker do the on-scanner processing? What are the coil combination amplitudes? How is the receiver gain applied? Any other scalings that we're missing?
- (Jan 5 GO) Are the `fid` (combined and processed) data already eddy-current-corrected? There appears to be a parameter `EDC_OnOff` indicating this but it does not appear in all data.
- (Jan 5 GO) For now, I'm applying the frequency shift (corresponding to the difference between FrqRef and FrqWork) to the *water-suppressed* data (PV-360 only) or to the ref data (if not PV-360), but I'm not confident that this is correct.
- (Jan 5 GO) Is cutting off the number of points indicated in GRPDLY really sufficient? When I look at the FIDs, I frequently find that the top of the echo appears slightly after that chop-off point

### TO DO

- [ ] (Jan 7 GO) If combined data are requested but not fid file is found, look for ser file (and return that)
- [x] (Jan 5 GO) Test the extended dataset (DP20-DP32)
- [x] (Jan 5 GO) Test raw data loading
- [ ] (Jan 3 GO) Design tests to check that data are actually loaded *correctly*, i.e., do some sanity checks, visual inspection, etc.
- [ ] (Jan 3 GO) Design tests for reference data (in the `mrsref` directories).
- [x] (Jan 3 GO) I have found some interesting code in Jessie's loader that does some scaling according to the receiver gain - need to look into that
- [x] (Jan 3 GO) Need to figure out where to find amplitude and phase coefficients for coil combination

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
