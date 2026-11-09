# Protocol: native Matlab + Letswave 6

This is an independent **diagnose-first** branch.  It stops after ICA fitting and diagnostics. It does not yet produce reviewed/final cleaned epochs.

## Cluster commands

Copy/extract this package into:

`/YOUR_PATH`

Then on the login node:

```bash
cd ~/
mkdir -p ~/logs
bash setup_letswave.sh
sbatch run_one.sbatch
```

`setup_letswave.sh` downloads the official Letswave repository into a **separate** `letswave6_moaz_pin` folder and checks out commit first. Internet access and git are needed for this one-time setup. It refuses a modified existing checkout. No original data are downloaded or modified.

The job loads `matlab` and requires its Signal Processing Toolbox (`butter`, `filtfilt`). Bundled Letswave import, coordinates, ICA and plotting dependencies come from the pinned checkout. The Python `rsa` environment is not used.

If your available Matlab module has a versioned name, first check `module spider matlab`, then submit with that name, for example:

```bash
MATLAB_MODULE=matlab/<your_available_version> sbatch run_one.sbatch
```

Do not copy the literal placeholder. Environment overrides available: `MOAZ_SCRIPT_DIR`, `MOAZ_LW_ROOT`, `MOAZ_BDF`, `MOAZ_OUTPUT`, `MATLAB_MODULE`. Defaults are set in `run_one.sbatch`. The pilot is intentionally restricted to subj04 S1. Do not point it at a different person's BDF while retaining this label.

Output by default:

`~/preprocessed_m_JOBID/subjxx/S*`

Every output directory must be new; an existing directory causes an error, including an incomplete previous run. An 8-hour request is a wall-time limit, not a predicted runtime. The job exits as soon as it finishes. Retaining intermediate native datasets requires several GB of disk space.

## Executed steps and recorded assumptions

1. Native `RLW_import_BDF`; require 512 Hz and one continuous recording.
2. Native international 10/20 template coordinates from `Standard-10-20-Cap81.locs`; remove EX1–8 / EXG1–8 and Status. Require 64 unique scalp labels and coordinates for all. Confirm the BDF physical dimension is uV; no guessing or silent unit conversion.
3. Remove only codes 251, 253, 254, 255 and 80. Retain 252. Save original/deleted-event ledgers. No automatic ghost-event repair; the pilot manual deletion list is empty.
4. Native continuous Butterworth bandpass 0.1–40 Hz. Requested order 4 is the pinned Letswave GUI default, **not an order specified by Moaz**. This version internally uses a second-order bandpass prototype with `filtfilt`.
5. For each 202–217 trigger, segment at −0.02 s for 291.1 s. This duration is fixed from Moaz's text, not recomputed from adjacent trigger intervals. Missing/duplicate block markers cause an error.
6. Within each block, segment ordinary stimuli 1–60 at −0.1 s for 1 s. The native inclusive endpoint gives **513 samples**. Keep the native event propagation. The pilot entrypoint explicitly selects `drop_incomplete_with_audit`: trials whose complete windows exceed a fixed block are excluded and logged in `qc/boundary_excluded_trials.csv`. The pinned native segmenter packs valid trials first but leaves preallocated zero tail epochs; this wrapper trims only that verified zero tail, updates the count, and never passes artificial zero trials to ICA. This is a disclosed compatibility repair, not an assertion of historical Moaz behavior. Direct configs without this policy still stop on incomplete trials.
7. Native whole-trial DC removal, **without** linear detrending.
8. Native amplitude rejection across all retained scalp channels and the complete 1-second trial: reject if **any absolute value > 700 uV after DC removal**. Exactly ±700 is retained. This is not a 700-uV peak-to-peak limit. Record every trial, identity and stimulus race before/after rejection.
9. Optional visual rejection and neighbouring-channel interpolation are **not applied**, because no reviewed list was supplied. No global bad channels are silently inferred. This follows the automatic amplitude-criterion alternative in the supplied procedure.
10. Merge the retained trial datasets across all 16 blocks **within S1 only**; call native `LW_ICA_compute_merged` with its default square `runica` configuration. No extra 1-Hz ICA filter, no notch, no resampling, no boundary guard from v2.2, no reference/baseline before ICA.
11. Save the matrix and native per-block ICA datasets; plot all components. **Exclude nothing and stop for review.**

The GUI normally saves intermediate datasets as single precision. This wrapper uses native `CLW_save` and reloads each stage so these storage round trips and native header sanitation are retained. Native sample/latency conventions are not corrected to MNE conventions. History for the ICA matrix is saved by the native merged-ICA wrapper. The complete scripted protocol is recorded in `run.json`, this README and `config_used.json`; intermediate RLW calls do not create full GUI menu-history entries.

## Limits on the word “strict”

This reproduces the supplied operation order using pinned **native Letswave functions**, rather than substituting MNE algorithms. It cannot establish numerical equivalence to Moaz's historical run until his Letswave version and actual settings are known.

Two source-version details are deliberately disclosed:

- 291.1 s gives **149044 samples** in this native version (`round(duration * 512) + 1`), while Moaz also wrote 149043. The pilot follows the written duration. To choose the sample-count interpretation, an explicit protocol revision is needed; it is not silently substituted.
- The pinned `RLW_ICA_compute` initializes a full zero block and then appends the real epochs. Thus its ICA input includes a zero prefix equal in length to the real retained data. We call the native code unchanged and record this behavior; we do not assert Moaz's historical release had it. Native `runica` also reseeds from the clock, so a second fit need not have identical component numbers/signs. Review and reuse the specific saved matrix.

Sources: [official Letswave repository](https://github.com/NOCIONS/letswave6/tree/7932ccfd570dc55ba1c5ae8e062147afce615496), [filter](https://github.com/NOCIONS/letswave6/blob/7932ccfd570dc55ba1c5ae8e062147afce615496/RLW/RLW_butterworth_filter.m), [segmentation](https://github.com/NOCIONS/letswave6/blob/7932ccfd570dc55ba1c5ae8e062147afce615496/RLW/RLW_segmentation.m), [ICA](https://github.com/NOCIONS/letswave6/blob/7932ccfd570dc55ba1c5ae8e062147afce615496/RLW/RLW_ICA_compute.m).

## What to inspect after the job

`run.json` should say `awaiting_ica_review`. This means diagnosis finished, not that cleaned data were approved.

- `qc/block_counts.csv`: normal-face count before/after rejection per block; oddball count within the block.
- `qc/boundary_excluded_trials.csv`: identity and native stimulus latency of each incomplete trial excluded at a fixed block boundary. Block summaries distinguish boundary from amplitude rejection.
- `qc/identity_block_counts.csv`: each identity's retained repetitions; no forcing counts back to four.
- `qc/trial_ledger.csv`: every ordinary trial, original native latency, race, maximum after DC and acceptance.
- `qc/events_original.csv`, `events_deleted.csv`, `block_triggers.csv`: event audit.
- `ica_diagnostics/topographies_*.png`: all component maps, 16 per page.
- `ica_diagnostics/ICA001.png` … `ICA064.png`: all component maps, trial images, per-trial spectra and trial means.
- `ica_diagnostics/frontal_correlation_ranking.csv`: Fp1/Fp2 correlation ranking **only**. These scalp channels are not confirmed EOG, and ranking does not authorize deletion.
- `ica_matrix.mat`: saved matrix, fit duration, RNG states and native configuration.
- `pre_ica_cache.mat`: retained per-block trials for a later reviewed application, avoiding a refit.
- `datasets/09_ica_*.lw6` with paired `.mat`: native files containing the matrix history, usable for Letswave inspection.
- `ica_review_template.json`: intentionally `review_confirmed: false` and no exclusions. **Numbers are native ONE-BASED**: ICA1 is component 1. Do not transfer a v2/MNE index [0] without reviewing this new decomposition.

After review, the remaining Moaz steps are ICA application using the saved matrix, DC + linear detrend, all-scalp average reference, −0.1…0 baseline, and the final identity/race segmentation/ERP operations. They are deliberately not executed in this pilot. A separate reviewed application step is needed; editing the template alone does not launch cleaning.

## Verification

The package has passed local Matlab static syntax/lint and Bash syntax checks. The development machine has no Matlab or cluster BDF, so no real-data execution is claimed. The batch job first runs `moaz_selftest.m` on small synthetic data to check native filter-order semantics, inclusive segmentation, DC, the exact 700-uV criterion and single-precision storage, then proceeds to the BDF. Source paths and SHA-256 values for critical routines are checked against `source_lock.json`; a mismatch stops the job.
