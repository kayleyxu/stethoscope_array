# Stethoscope Array: Acoustic Source Localization and Differential Beamforming for Body Sounds

MATLAB code and data for research on **multi-microphone (MEMS) stethoscope arrays**. A single stethoscope records one point on the body. An array of microphones adds spatial information, which lets you:

1. **Localize** a sound source inside the body in 3-D, such as a lung crackle. This uses steered-response-power delay-and-sum (SRP-DAS) beamforming over a 3-D search grid. Sparse array geometries (Costas arrays, Golomb rulers) are compared with regular layouts, and a hybrid beamforming + machine-learning estimator is also studied.
2. **Separate** sounds that overlap at the chest wall. A two-mic end-fire differential beamformer places a spatial null on the interferer, such as heart sound masking lung sound or maternal heartbeat masking fetal heartbeat. The change in signal-to-interference ratio (ΔSIR) is measured.
3. **Characterize the hardware and medium.** The repo compares the noise spectra of the custom MEMS microphones with a commercial digital stethoscope, and measures the speed of sound in tissue-mimicking gelatin phantoms.

Each topic is worked through in simulation first, then validated **in vitro** with recordings from a gelatin foam/gel phantom.

## Repository layout

```
src/    MATLAB scripts and functions (run from inside src/; paths are relative, e.g. ..\data\...)
data/   Raw multichannel recordings (.raw, int16 little-endian), metadata, generated results
```

| Data folder | Contents |
|---|---|
| `data/2024_11_26_23_51_25_single_crackle/` | 16-mic recording of a single crackle in air (48 kHz). Also used as the crackle template in simulations |
| `data/localization_test1` … `localization_test8/` | 5-mic Costas-array in-vitro localization recordings (64 kHz) with `metadata.mat` |
| `data/localization_analysis_results/` | Output of `invitro_loc_experiments.m`: per-set error tables, paper table, figures |
| `data/sectionB_invitro_recording/` | In-vitro differential-beamforming recordings: 5 heart-position conditions × 3 repeats (heart/lung/noise WAVs, manifest, SIR verification CSVs) |
| `data/sound_speed_3/` | Two-mic recording used to estimate sound speed in the phantom |
| `data/stethoscope_rec/` | White/pink-noise chamber recordings from the commercial stethoscope and the MEMS mics, for PSD comparison |
| `data/outputs/sound/` | Before/after audio demos of the fetal-vs-maternal heartbeat beamforming |

> Some scripts may need the `folder_data` / file-path lines at the top changed to point to your data folder.

## Requirements

- MATLAB with the **Signal Processing Toolbox**, the **Phased Array System Toolbox** (`phased.ULA`, `diffbfweights`, `delayseq`), the **Statistics and Machine Learning Toolbox** (`fitrensemble`, `fitrgp`, `fitrsvm`, `cvpartition`) and the **Communications Toolbox** (`awgn`).
- Most scripts are experiment notebooks. You switch scenarios by commenting or uncommenting the parameter blocks near the top of the file.

## Code overview

### 1. 3-D source localization: real recordings

- **[`mic_preprocess.m`](src/mic_preprocess.m)**: Preprocessing function for the 16-mic array (48 kHz, int16 `.raw`).
  - Loads the channels, with optional coarse cropping.
  - Finds the most active common 20 ms window.
  - Applies the calibration delay compensation.
  - Computes each mic's relative latency (TDOA) by cross-correlation and plots waveforms and STFTs.
- **[`localization.m`](src/localization.m)**: 16-mic localization in air.
  - Calibration: a source at `[0, 0, 0.22]` m gives each mic's relative latency, measured by cross-correlation.
  - SRP-DAS localization with the calibration applied as delay compensation.
  - Source locations: `[±0.08, 0, 0.22]` and `[0, ±0.08, 0.22]` m.
  - Calls `mic_preprocess.m`.
- **[`invitro_mic_preprocess.m`](src/invitro_mic_preprocess.m)**: Same steps as `mic_preprocess.m`, set up for the 5-mic in-vitro array at 64 kHz.
- **[`invitro_localization.m`](src/invitro_localization.m)**: Like `localization.m`, but uses the 5-mic Costas array on a gelatin phantom. Calls `invitro_mic_preprocess.m`. MVDR and MUSIC variants are included but commented out.
- **[`invitro_loc_experiments.m`](src/invitro_loc_experiments.m)**: Standalone batch analysis of the eight recordings `localization_test1..8` (5 mics, Fs = 64 kHz, c = 30 m/s).
  - Runs a 3-D SRP-DAS grid search against the measured reference positions.
  - Reports x/y/z and total errors.
  - Writes CSV/MAT summaries, a paper-style table and figures to `data/localization_analysis_results/`.

### 2. Array geometry: simulation

- **[`costas_property.m`](src/costas_property.m)**: Shows the 5×5 Costas array geometry and its auto-correlation and displacement properties. These properties are why the layout gives unambiguous delay patterns.
- **[`golomb_circle.m`](src/golomb_circle.m)**: Monte-Carlo comparison of a 5-mic circular array with Golomb-ruler angular spacing against one with even angular spacing.
  - Uses grid-search localization of a Gaussian-modulated pulse in noise.
  - Reports mean, median and variance of the error.
- **[`costas_mic_array_src_trial_with_new_1_offset.m`](src/costas_mic_array_src_trial_with_new_1_offset.m)**: Robustness study of the array layout.
  - 21 six-mic layouts. Each shares one fixed side mic: 1 Costas layout, 10 with random `{-0.5, 0, 0.5}` mic shifts and 10 with random `{-1, 0, 1}` mic shifts.
  - Each layout gets 10 random sources inside the Costas contour × SNR `[-3, 0, 3]` dB × `[150, 250, 400]` Hz × 2 runs = 180 trials.
  - Uses SRP-DAS with an exhaustive grid search for the maximum-energy point.
  - Metrics: RMSE (which penalizes large errors) and variance.
  - Switching between 5 mics (no side mic) and 6 mics shows how much the side mic adds to depth estimation.
  - Calls `costas_mic_array_src_trial_display.m`.
- **[`costas_mic_array_src_trial_display.m`](src/costas_mic_array_src_trial_display.m)**: Error statistics tables by layout, SNR and frequency, plus plots of the mic layouts and predictions.

### 3. Beamforming-assisted machine learning

- **[`costas_mic_array_dataset_gen.m`](src/costas_mic_array_dataset_gen.m)**: Builds simulated datasets for the Costas array: a small one (180 samples) and a large one (3200 samples).
  - Each row of `results_for_ml` holds `source_position`, `snr`, `freq`, `c`, `predicted_position_srp`, `adjusted_delays` (cross-correlation TDOAs), `position_sets` and `raw_mic_signals`.
  - Saved to `src/simulation_results_ml/`.
- **[`ml_location.m`](src/ml_location.m)**: Compares three localization approaches on an 80/20 train/test split:
  - **SRP-DAS only**
  - **ML only**, trained on raw mic signals
  - **BF-assisted ML**, trained on SRP-DAS estimates and TDOA features

  Regressors: LSBoost ensemble, Gaussian process and SVM.

### 4. Differential beamforming for sound separation

- **[`lung_heart_differential_stft.m`](src/lung_heart_differential_stft.m)**: Two-mic end-fire differential array (5 mm spacing) in simulation.
  - Synthetic breathing with a recorded crackle, plus a synthetic heartbeat. The two sources sit at +90° and −90°, so both mics hear a mixture.
  - The STFT splits the wideband signal into narrowband bins. Differential beam weights are computed per bin, applied and summed, then the signal is rebuilt with the ISTFT.
  - Two weight sets are computed, one with the null at +90° and one at −90°. The weight set chosen decides whether the lung sound or the heart sound is cancelled.
- **[`invitro_lung_heart_differential_sfft.m`](src/invitro_lung_heart_differential_sfft.m)**: The same method applied to in-vitro two-mic recordings from the gelatin phantom.
- **[`fetal_maternal_heart.m`](src/fetal_maternal_heart.m)**: Uses differential beamforming to suppress the maternal heartbeat when monitoring the fetal heartbeat. Writes before/after audio to `data/outputs/sound/`.
- **[`DiffBF_sim_and_invitro.m`](src/DiffBF_sim_and_invitro.m)**: Full SIR study, with the results formatted for a paper.
  - **Section A (simulation):** sweeps over azimuth, depth and common-noise level.
  - **Section B (in vitro):** noise-corrected SIR before and after beamforming, and ΔSIR, for 5 heart positions (−120°/−90°/−60° at 5 cm, and −90° at 3 cm and 7 cm), each with N = 3 repeats.
  - Matched simulation-vs-in-vitro tables.
  - Exports PNG/PDF/FIG figures, CSV and LaTeX tables, and a MAT file.

### 5. Hardware and medium characterization

- **[`welch.m`](src/welch.m)**: Welch power spectral density comparison in a noise chamber. Compares a commercial digital stethoscope (Littmann CORE; the recordings are Eko app exports) with the custom MEMS microphones (MP and CMM variants).
- **[`sound_speed.m`](src/sound_speed.m)**: Estimates the speed of sound in the gelatin foam/gel phantom from the time delay between two mics.

### Bundled audio

`src/heartbeat_signal_mom.wav`, `src/heartbeat_signal_fetus.wav`, `src/mic1_signals.wav` and `src/mic2_signals.wav` are source and mixture signals used by the beamforming demos.
