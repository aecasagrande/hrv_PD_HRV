# Changelog

## 2.6.1 — 2026-10-05 (code review + end-to-end simulation)

The whole of `main` (3,935 lines) was reviewed for logic errors, edge cases and crashes, with
particular attention to array indexing, NaN / empty-data handling (trials too short, every
interval flagged, < 5 beats) and the NumPy/Pandas logic. The detection and reconstruction core
(Kubios detector port, [D13]–[D28] rules, run planning, GP fill, beat breakdown) held up: its
index arithmetic is consistent, short and degenerate inputs return empty (NaN) results rather
than crashing, and per-trial errors are caught so one trial cannot stop a batch.

Problems were in the data-loading and spectral code. **No time-domain value changes on
normal input** (HR, SDNN, RMSSD, SDSD, pNN50 and every Data_Quality / Reconstruction_QC count
are identical to 2.6 on the synthetic recordings). `PIPELINE_VERSION` is now `2.6.1`, so
Run_Config shows which version wrote a workbook.

### Fixed

1. **EDF files with different sampling rates were silently mixed** (`load_recording`).
   The warning said "using per-file rate", but every RR interval is measured in samples and
   converted with the first file's rate, so the other file's intervals were wrong by the
   ratio of the two rates (250 vs 256 Hz: 2.4 %). Such a file is now resampled to the
   first file's rate.
   *Test:* a 75-bpm metronome rhythm split over a 256-Hz and a 250-Hz file. 2.6 gave
   HR 73.0 bpm and SDNN 156 ms for a trial reading both files; 2.6.1 gives 75.0 bpm and
   1.6 ms, the same as a trial in one file.

2. **Overlapping EDF files were interleaved sample by sample** (`load_recording`).
   The files were concatenated and sorted by time stamp, and only *exactly equal* time stamps
   were removed. Where two files overlap with sample grids that don't line up (the
   Equivital's ~256.0164-Hz clock, files starting on different whole seconds, or different
   rates), the samples of both files were merged. That doubles the sample density and
   corrupts the ECG and every sample-based RR interval there. Now the later file's samples
   inside an already-covered span are dropped, with a console warning. Nothing changes when
   files don't overlap. *Test:* the same pair of files overlapping by 10 s. This, together
   with item 1, caused the error above.

3. **Sync-log times after midnight** (`parse_log_file_new`). Only the time of day is
   read from the log, and the date comes from the EDF. A trial that crosses midnight
   stopped before it started, and trials after midnight landed ~24 h *before* the
   recording. Both were skipped as "window outside recording". Times more than 12 h before
   the recording start now go on the next day, and a stop earlier than its start is moved
   forward a day. `process_trial` now also gives a clear message for a stop time that is not
   after the start. *Test:* a session 23:50–00:10 with three trials after midnight. 2.6
   processed 1 of 4 trials, 2.6.1 processes 4 of 4.

4. **Duplicate trial names overwrote each other** (`parse_log_file_new`). Processed trials
   are stored in a dict keyed by name, so a trial logged twice (e.g. a repeat) silently
   replaced the earlier one in every sheet. Later repeats are now named `<trial>_dup2`,
   `_dup3`, … with a console warning. *Test:* `OFF_30s_T1` logged twice gives both trials.

5. **Spectral metrics that cannot be computed were written as 0** (`calculate_lomb_scargle_psd`,
   `calculate_welch_psd`, `calculate_hrv_metrics`, `export_hrv_to_excel`). This is the same
   flaw [D17] fixed for the time domain: every 30-s trial (shorter than
   `MIN_SEGMENT_DURATION`), and every trial with too many exclusions, had VLF = LF = HF =
   LF/HF = 0 in Freq_Lomb_Scargle / Freq_Welch, which reads as a real value. LF/HF, n.u. and
   relative powers were also 0 when the denominator was 0. These cells are now **empty
   (NaN)**, and the `Quality_*` column still says why. The R script reads only Time_Domain
   and Data_Quality (per the comment above `QC_COLUMNS`), so it is unaffected.
   The dashboard shows "—" for missing numbers instead of `nan`. Its benchmark table showed
   a missing Kubios value as `0.00` and crashed on a `None`; both are fixed.

6. **Lomb–Scargle detrending used the sample index instead of time** (`calculate_lomb_scargle_psd`).
   `scipy.signal.detrend` fits its line against the array index. The Lomb–Scargle input is
   deliberately unevenly sampled (excluded beats, gaps; see the Methods statement), so the
   index is not proportional to time and the removed trend was slightly wrong. The line is now
   fitted against the beat times. Effect on the synthetic 5-min trials: ≤ 0.04 % on LF, HF and
   LF/HF. Welch is unchanged, since it detrends on its uniform 4-Hz grid where index and time
   agree.

7. **Spectra spliced recording holes shut** (`calculate_hrv_metrics`, `process_trial`).
   Beat times for the spectra were sample numbers / fs. Across a hole in the recording [D26]
   the samples on either side are adjacent, so the hole had zero duration and the series
   after it was shifted back in time. That contradicts [D6] ("gaps stay gaps in time") for
   exactly the case [D26] targets (5-min quiet trials of P-018 and P-028 have such holes).
   `calculate_hrv_metrics` now takes the segment's sample time stamps (`sample_index`) and
   uses clock times. Time-domain metrics are unchanged. *Test:* 300 s of beats with a 20-s
   hole: the spectral duration goes 299.2 s → 319.2 s (correct), and RMSSD is identical.

8. **EDF without a start time** (`load_recording`). A file without `meas_date` used to be
   placed at *the current time*, which can never match the sync log. It is now skipped with
   a clear message.

### Reviewed, not changed (worth knowing)

- **[D28] is referenced but not documented.** `ALTERNATION_HANDLING = 'auto'` (the default),
  `_alternation_sustained`, `ALTERNATION_CONTEXT_BEATS` / `_SUSTAINED_FRAC` and the
  early-beat counting in `plan_correction_runs` all cite [D28]. The DECISIONS list stops at
  [D27], and [D21] still describes only `'keep'` / `'off'`. The rationale and validation for
  [D28] should be added before the methods section is written. I did not write it, because
  the validation behind the choice is not in the repository.
- `kubios_find_artifacts` advances `i` once per long/short *candidate*, so it can step two
  beats at once. This is inherited from NeuroKit2's implementation, and it is left as is so
  the port stays faithful to the benchmark reference.
- 'Raw' metrics still use the detected series as it is, including an interval across a
  recording hole (documented in [D26]).
- Welch interpolates across excluded stretches. This is documented in the Methods statement,
  and Welch is the secondary estimate.
- A beat within 0.2 s of the very start or end of the *loaded recording* (not the trial)
  has no full morphology epoch and is flagged "abnormal shape". This is conservative, and
  with `PAD_SECONDS` of padding it only affects trials at the edge of a file.
- `remove_small_peaks` assumes the R peaks are positive after `INVERT_ECG`. If a file is
  recorded with the opposite polarity, its amplitude threshold is meaningless. Check the
  dashboard for any new participant.

### Added

- `generate_mock_data.py` uses numpy / scipy / mne (with edfio as MNE's EDF writer).
  It writes synthetic EDF recordings, sync logs and trial run sheets under `./mock_data` for
  **P-001** and **P-015**, laid out exactly as their `PARTICIPANTS` entries expect, so the
  registry's special cases are exercised:
  - The rhythm is sinus with respiratory sinus arrhythmia (0.25 Hz), Mayer waves (0.1 Hz)
    and random variability. The signal has baseline wander (respiration plus slow drift),
    60-Hz mains and broadband noise, and is written with inverted polarity like the
    Equivital.
  - Events are injected at known times: PVCs with compensatory pauses, sinus-resetting
    APCs, missed beats (tiny R wave), a QRS-band noise burst, and a 60-s bigeminy run
    (P-015 ON standing).
  - P-001 OFF has a 10-s recording hole inside `OFF_30s_T2` [D26]. P-001 ON has its 5-min
    sitting in a separate `PD_P001_ON_State_Rest(1).edf` [D25].
  - P-015 OFF has a corrupted ECG 1 (motion noise and dropouts) during `OFF_5min_Sitting`,
    which the registry reads from ECG 2 [D27]. Its OFF file also covers an `ON_` trial
    [D20]. P-015 ON has a 2-s placeholder trial.
  - Every beat's true time and type is saved as `*_truth.csv`. The generator is seeded, so
    re-running it gives byte-identical EDFs.
- `validate_mock_results.py` compares the workbooks with the ground truth: true RMSSD/SDNN
  from adjacent sinus–sinus intervals versus the pipeline's raw and corrected values. It also
  checks the [D23] identity `Unique_Bad_Beats = Dropped + Corrected + Kept` and, optionally,
  compares against another run.
- `run_simulation.py` runs the full batch on `./mock_data` → `./output` without editing
  `main`, then validates.
- `output/HRV_Analysis_P-00{1,15}_{OFF,ON}.xlsx` are the four workbooks from the simulation.

### Simulation results (`BASE = ./mock_data`, `OUTPUT_DIR = ./output`, `BATCH_MODE = True`, `DRY_RUN = False`)

The batch ran over all 45 registry entries in ~30 s. It wrote 4 workbooks with 20 trials and
9 sheets each. The 41 entries without mock data were skipped as "no trials", with no errors.
Each special case behaved as intended:

| Case | Expected | Observed |
|---|---|---|
| [D26] 10-s hole in P-001 `OFF_30s_T2` | interval across the hole excluded, nothing inserted | `recording gap` run excluded, `Unrecorded_Seconds` = 10 |
| [D25] P-001 ON sitting in the Rest file | trial found and analysed | 353 beats, RMSSD 36.8 ms (true 36.9) |
| [D27] P-015 `OFF_5min_Sitting` from ECG 2 | clean lead used | `ECG_Channel` = ECG 2, 3.7 % bad beats |
| [D20] `ON_30s_T9` inside the OFF file | skipped for OFF | skipped, with the D20 note |
| [D21/D28] 60-s bigeminy in P-015 ON standing | sustained alternation excluded | 20 % alternating, 84 beats excluded, RMSSD 31.5 ms (true 33.2) |
| placeholder `ON_IART_Baseline` (2 s) | skipped | listed by the preflight, skipped |

Accuracy against the ground truth over the 20 trials:

| | raw | corrected |
|---|---|---|
| mean \|RMSSD error\| | 210 % | **4.5 %** (median 2.7 %) |
| mean \|SDNN error\| | 130 % | **3.7 %** (median 2.2 %) |

The worst trial is P-015 `OFF_30s_T3`, with RMSSD −26 %. It is a 30-s trial whose PVC is
replaced by smooth GP values, which is the known small-trial underestimate described in
[D2] and VALIDATION 9.

The dashboard callback was also run for every trial with both PSD methods, and the DRY_RUN
preflight on the mock data. The preflight lists the hole, the placeholder and the trial
outside the files.

### How to reproduce

```bash
pip install mne neurokit2 scikit-learn dash plotly openpyxl scipy numpy pandas edfio
python generate_mock_data.py
python run_simulation.py          # writes ./output and prints the validation table
```

For the simulation, `BASE`, `OUTPUT_DIR` and `DRY_RUN` in `main` were temporarily set to
`./mock_data`, `./output` and `False`. They are restored to the lab-drive paths and
`DRY_RUN = True` in the committed script. `run_simulation.py` applies the same settings
without editing the file.
